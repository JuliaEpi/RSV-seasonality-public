using ArgParse, DrWatson
using AdvancedHMC.ProgressMeter
using StatsPlots
using RSVSeasonality
using CSV, DataFrames, StatsPlots, DifferentialEquations, Turing, DiffEqSensitivity
using Distributions, StatsBase, LabelledArrays, MCMCChains, Serialization
using Random
using Dates
using DynamicPPL
using ModelingToolkit


include(scriptsdir("utils.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "runs"
    help = "runs to process"
    action = :store_arg
    nargs = '*'
    "--prefix"
    help = "Prefix that name of run should match, e.g. '2021-08-19T14' for a runs started on 19th of August 2021 between 2pm and 3pm. Either this or runs should be provided."
    "--list-runs"
    help = "If specified, available runs will be listed."
    action = :store_true
    "--verbose"
    help = "If specified, more logging information will be provided."
    action = :store_true
    "--ignore-commit"
    help = "If specified, checkout of the correct commit won't be forced."
    action = :store_true
end

include(scriptsdir("parse_args.jl"))
include(scriptsdir("trace_merge.jl"))

const verbose = args["verbose"]

# Get the available runs. Useful for either `--list-runs` or if `--prefix` is given.
available_runs = filter(
    isdir ∘ Base.Fix1(projectdir, "output"),
    readdir(projectdir("output"))
)

# List runs and exit, if specified.
if args["list-runs"]
    println("The following runs are available:")
    println()
    for run in available_runs
        println("  $(run)")
    end
    exit(0)
end

# Get the absolute paths of the runs.
runs = if length(args["runs"]) > 0
    map(abspath, args["runs"])
elseif !isnothing(args["prefix"])
    prefix = args["prefix"]
    map(Base.Fix1(projectdir, "output"), filter(startswith(prefix), available_runs))
else
    error("either runs or prefix has to be provided")
end
verbose && @info "Using runs $(runs)"
@assert length(runs) > 0 "no runs specified"

# Load model, trace and model args
verbose && @info "(♻) Loading model and chain(s)..."
model, trace, run_args = trace_merge(runs, args["ignore-commit"])
verbose && @info "(✓) Done!"

# Predict 
verbose && @info "(♻) Predicting..."
# `decondition` model so that `observations` are now treated as random variables.
predictive_model = if (model.name == :turingmodel_noar)
     DynamicPPL.decondition(model, :obs_ts);
else
    DynamicPPL.decondition(model, :obs_ts, :obs_ar1, :obs_ar2, :obs_ar3);
end

# We also want to simulate for _all_ timesteps (also time points with missing data), not just those that we have data for.
times_predict = let foo = model.args.times 
    # Use the smallest observed timestep as the "default" timestep.
    timestep = minimum(diff(foo))
    collect(first(foo):timestep:last(foo))
end
# Update the indices to be used in the incidence vector
times_ind_predict = convert.(Int, (times_predict ./ 7) .+ 1)

# Update `predictive_model` to use `times_predict`.
predictive_model = RSVSeasonality.replace_args(predictive_model; times = times_predict, times_ind = times_ind_predict);

# Extract the parameters from `trace`.
parameters = DataFrame(MCMCChains.get_sections(trace, :parameters))

function generate_output(model, retval)

    sol = DataFrame(retval.sol)
    rename!(sol,  [(:time); collect(keys(model.args.problem.u0))])
    # Subset: We're only interested in those at `time ≥ 0`.
    sol = subset(sol, :time => ByRow(≥(0.0)))

    incidence = DataFrame(retval.inc_w, :auto)
    rename!(incidence, 1:4 .=> ["inc"] .* string.(0:3))
    incidence[!, :time] = times_predict

    inc_symp = DataFrame(retval.inc_symp, :auto)
    rename!(inc_symp, 1:4 .=> ["inc_symp"] .* string.(0:3))
    inc_symp[!, :time] = times_predict

    # Combine with solution
    traj_sim = outerjoin(sol, incidence, inc_symp, on=:time)

    # Calculate seasonal beta and reporting rate
    betafunc = model.args.problem.f.f.betafunc
    rhofunc = model.args.rho_func
    p = retval.p
    @. traj_sim[!, :betaseason] = betafunc(traj_sim.time, p.φ, p.η, p.β₀, p.k)
    @. traj_sim[!, :rho] = rhofunc(p.ρ₀, p.q * traj_sim.time)

    # Calculate observable and reported incidence
    @. traj_sim[!, :inc_obs] = traj_sim.inc_symp0 + traj_sim.inc_symp1 + traj_sim.inc_symp2
    @. traj_sim[!, :inc_rep] = traj_sim.inc_obs * traj_sim.rho
        
    # Set observations
    traj_sim[!, :simobserror] .= retval.obs_ts

    if (model.name == :turingmodel_noar)
        return traj_sim
    else
        # Return also sampled AR numerators
        ar_sim = DataFrame(hcat(retval.obs_ar1, retval.obs_ar2, retval.obs_ar3), [:ar1, :ar2, :ar3])
        return traj_sim, ar_sim
    end

end

# Simulate from posterior
traj = DataFrame()
if !(model.name == :turingmodel_noar) ar = DataFrame() end
parameter_iterator = Tables.namedtupleiterator(parameters)
@showprogress "Generating..." for (i, θ) in enumerate(parameter_iterator)
    # Condition `model` on the parameters `θ``
    cmodel = predictive_model | θ
    # We wrap this in a `try-catch` because the call to `rand` for generating the
    # observations might fail if we have a particularly badly behaved `θ`.
    try
        # Sample from predictive.
        retval, issuccess = cmodel()
        if !issuccess
            @warn "iteration $i does not result in a successful evaluation of the model; skipped!"
        else
            if !(model.name == :turingmodel_noar)
                traj_sim, ar_sim = generate_output(cmodel, retval)
                traj_sim[!, :replicate] .= i
                traj_sim[!, :chain] .= parameters.chain[i]
                append!(traj, traj_sim)
                ar_sim[!, :replicate] .= i
                append!(ar, ar_sim)
            else
                traj_sim = generate_output(cmodel, retval)
                traj_sim[!, :replicate] .= i
                traj_sim[!, :chain] .= parameters.chain[i]
                append!(traj, traj_sim)
            end
        end
    catch e
        @warn "iteration $i resulted in an error when evaluating the model; skipped!"
    end
end

traj[!, :rhopct] .= traj.rho .* 100
verbose && @info "(✓) Predictions done!"

# Replace subscripts for saving as CSV
rename!(traj, [if !isascii(x) replace(x, x[2] =>  x[2] - '₀') else x end for x in names(traj)])

# Trace plots and diagnostics
verbose && @info "(♻) Generating plots..."
traceplot = plot(trace)
lpplot = plot(trace[:lp])

# Get diagnostics (ESS / RHAT)
ess_rhat_df = DataFrame(ess_rhat(trace))

#=internals = filter(
    !=(:max_hamiltonian_energy_error),
    trace.name_map.internals
)
diagplot = plot(MCMCChains.get_sections(trace, :internals)[internals])


p_diagnostics = plot(diagplot, traceplot, layout=(1,2), size=(1600,2000))
=#

verbose && @info "(✓) Plotting done!"

# Save output
verbose && @info "(♻) Saving output..."
suffix = if (model.name == :turingmodel_noar)
    string(run_args["betafunc"], "_", run_args["startyear"], "_", run_args["endyear"], "_", "noar")
else 
    string(run_args["betafunc"], "_", run_args["startyear"], "_", run_args["endyear"]) 
end 
save(string("output/figS3_", suffix, ".png"), traceplot)
save(string("output/figS4_", suffix, ".png"), lpplot)
CSV.write(string("output/ess_rhat_", suffix, ".txt"), ess_rhat_df)
CSV.write(string("output/traj_sim_", suffix, ".csv"), traj)
if !(model.name == :turingmodel_noar) CSV.write(string("output/ar_sim_", suffix, ".csv"), ar) end
CSV.write(string("output/posterior_quantiles_", suffix, ".csv"), quantile(trace))
CSV.write(string("output/trace_final_", suffix, ".csv"), DataFrame(trace))
verbose && @info "(✓) Finito!"
