using ArgParse, DrWatson, ProgressBars, Pipe, Random
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
    "--seed"
    help = "random seed to use"
    arg_type = Int
    default = 42
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
    "--simyears"
    help = "Duration of simulated vaccination in years"
    arg_type = Int64
    default = 10
    "--duration"
    help = "Vector of simulated duration of pulse in days (must be ∈ (0.0, 363.0])"
    arg_type = String
    default = "[210.0]"
    "--coverage"
    help = "Vector of simulated vaccination coverage(s) (must be ∈ (0.0, 1.0])"
    arg_type = String
    default = "[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]" # this cannot be parsed: "[0.1:0.1:1.0;]"

end

include(scriptsdir("parse_args.jl"))
include(scriptsdir("trace_merge.jl"))

const verbose = args["verbose"]
Random.seed!(args["seed"]);

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
model, trace, run_args = trace_merge(runs)
verbose && @info "(✓) Done!"

# Define simulation settings
simyears = args["simyears"] 
duration = Float64.(Meta.parse(args["duration"]).args) # duration of pulses
coverage = Float64.(Meta.parse(args["coverage"]).args) # vaccination coverage
vaccstart = [0.0:60.0:330.0;] .+ model.args.problem.tspan[1]

# Test arguments for consistency
(minimum(duration) > 0.0 && maximum(duration) <= 363.0) ? nothing : error("duration must be ∈ (0.0, 363.0]")
(minimum(coverage) > 0.0 && maximum(coverage) <= 1.0) ? nothing : error("coverage must be ∈ (0.0, 1.0]")

# Define model settings
statenames = keys(model.args.problem.u0)
estparnames = Tuple(trace.name_map.parameters)
tspan = (model.args.problem.tspan[1], simyears*365.0 + model.args.problem.tspan[1]) # start time has to be the same as for fitting to match the phase shift phi
betafunc = model.args.problem.f.f.betafunc

solvsettings = (
    abstol = 1e-6,
    reltol = 1e-6,
    solver = AutoTsit5(Rosenbrock23()),
    saveat = 1.0
)

# Get generated quantities: u0 at stable periodic orbit
quantities = generated_quantities(model, MCMCChains.get_sections(trace, :parameters));
u0 = map(quantities) do qty
    qty[1].u0
end;
u0 = vec(u0);
u0 = transpose(reduce(hcat, u0))

# Simulate vaccination scenarios --------------------------------------------------

verbose && @info "(♻) Sampling from posterior..."
posteriors = Array(trace, [:parameters])
nsamples = 500
draws = sort(sample(1:size(posteriors,1), nsamples, replace=false))

posterior_sample = posteriors[draws,:] 
u0_sample = u0[draws,:]

@info "(✓) Done!"


# Base scenario (no vaccination) --------------------------------------------------------
verbose && @info "(♻) Simulating base scenario (no vaccination)..."

# Define problem
problem_cont = ODEProblem(SEIRRS4Vacc!(; betafunc=betafunc))
theta_fix_cont = RSVSeasonality.fixedparameters(problem_cont, model.args.prior())


# Simulate from posterior draws
sim_base = DataFrame()
for i in ProgressBar(1:nsamples)
    # Get initial states and posterior estimates for given draw
    u0_new = problem_cont.u0
    u0_new[1:end-2] = u0_sample[i,:]
    theta_est_sample = @LArray posterior_sample[i,:] estparnames 
    p_new = vcat(theta_est_sample, theta_fix_cont)
    # simulate trajectory
    out = ODEwrap_vaccsim(problem_cont, p_new, u0_new, tspan, nothing, solvsettings)
    out[!,:replicate] .= i
    append!(sim_base, out)
end

# Calculate median and 95% CI from all sampled trajectories
quantiles_base = combine(groupby(sim_base, [:time]), 
        [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95CI", cols*"_median", cols*"_up95CI"]
            for cols in deleteat!(names(sim_base), sort(indexin(["time"; "vaccswitch"; "replicate"], names(sim_base))))])

quantiles_base[!,:strategy] .= "base"
quantiles_base[!,:vaccswitch] .= missing
quantiles_base[!,:coverage] .= missing
quantiles_base[!,:vaccstart] .= missing
quantiles_base[!,:duration] .= missing

CSV.write(string("output/vaccsim_base.csv"), quantiles_base)
@info "(✓) Done!"

# Continuous strategy -------------------------------------------------------------
verbose && @info "(♻) Simulating continuous vaccination..."

# Define callbacks
tstarts_cont = 0.0
tstops_cont = [tspan[2]]
cbstart_cont = PresetTimeCallback(tstarts_cont, integrator -> integrator.u.switch = 1.0, save_positions=falses(length(tstarts_cont)))
cbstop_cont = PresetTimeCallback(tstops_cont, integrator -> integrator.u.switch = 0.0, save_positions=falses(length(tstops_cont)))
cb_cont = CallbackSet(cbstart_cont, cbstop_cont)

# Sample from posterior and forward simulate vaccination at different coverages
sim_cont = DataFrame()
for i in ProgressBar(1:nsamples)
    for j in coverage
        # Get initial states and posterior estimates for given draw
        u0_new = problem_cont.u0
        u0_new[1:end-2] = u0_sample[i,:]
        theta_est_sample = @LArray posterior_sample[i,:] estparnames 
        p_new = vcat(theta_est_sample, theta_fix_cont)
        p_new.v = j  #update vaccination coverage
        # simulate
        out = ODEwrap_vaccsim(problem_cont, p_new, u0_new, tspan, cb_cont, solvsettings)
        out[!,:replicate] .= i
        out[!,:coverage] .= j
        append!(sim_cont, out)
    end
end

# Calculate median and 95% CI from all sampled trajectories
quantiles_cont = combine(groupby(sim_cont, [:time, :coverage]), 
        [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95CI", cols*"_median", cols*"_up95CI"]
            for cols in deleteat!(names(sim_cont), sort(indexin(["time"; "vaccswitch"; "coverage"; "replicate"], names(sim_cont))))])

quantiles_cont[!,:strategy] .= "continuous"
quantiles_cont[!,:vaccstart] .= vaccstart[1]
quantiles_cont[!,:duration] .= missing

# Add other vars back to dataset
quantiles_cont = @pipe sim_cont |> 
            select(_,[:time, :coverage, :vaccswitch]) |> 
            unique(_) |> 
            outerjoin(_, quantiles_cont; on = [:time, :coverage])

@info "(✓) Done!"


# Pulse/seasonal vaccination strategy --------------------------------------------------------
verbose && @info "(♻) Simulating seasonal (pulse) vaccination..."

# Define problem
problem_pulse = ODEProblem(SEIRRS4Vacc_pulse!(; betafunc=betafunc))
theta_fix_pulse = RSVSeasonality.fixedparameters(problem_pulse, model.args.prior())

# We cannot simulate all seasonal/pulse scenarios in one loop due to OOM error. 
# Instead, we create a function and call it on each element of the the vaccstart vector
function simulate_pulse(vaccstart)

    sim_pulse = DataFrame()

    for i in ProgressBar(1:nsamples)
        # Simulate all durations
        for l in duration
            # Define callbacks
            theta_fix_pulse.d = l # update duration
            tstarts_pulse = collect(range(vaccstart, step=365.0, length=simyears))
            tstops_pulse = collect(range(vaccstart + theta_fix_pulse.d, step=365.0, length=simyears))
            cbstart_pulse = PresetTimeCallback(tstarts_pulse, integrator -> integrator.u.switch = 1.0, save_positions=falses(length(tstarts_pulse)))
            cbstop_pulse = PresetTimeCallback(tstops_pulse, integrator -> integrator.u.switch = 0.0, save_positions=falses(length(tstops_pulse)))
            cb_pulse = CallbackSet(cbstart_pulse, cbstop_pulse)
            # Simulate all coverages
            for j in coverage
                # Get initial states and posterior estimates for given draw
                u0_new = problem_pulse.u0
                u0_new[1:end-2] = u0_sample[i,:]
                theta_est_sample = @LArray posterior_sample[i,:] estparnames 
                p_new = vcat(theta_est_sample, theta_fix_pulse)
                #update vaccination coverage
                # equation for pulse vacc is (log(1-v)/-d), therefore v cannot be 1. replace with 0.99999999. I know it's a bit hacky
                j==1.0 ? p_new.v = 0.99999999 : p_new.v = j
                # simulate
                out = ODEwrap_vaccsim(problem_pulse, p_new, u0_new, tspan, cb_pulse, solvsettings)
                out[!,:replicate] .= i
                out[!,:coverage] .= j
                out[!,:vaccstart] .= vaccstart
                out[!,:duration] .= l
                append!(sim_pulse, out)
            end
        end
    end
    
    # Calculate median and 95% CI from all sampled trajectories
    quantiles_pulse = combine(groupby(sim_pulse, [:time, :coverage, :vaccstart, :duration]), 
            [cols => (x -> (quantile(x, [0.025, 0.5, 0.975]))') => [cols*"_low95CI", cols*"_median", cols*"_up95CI"]
                for cols in deleteat!(names(sim_pulse), sort(indexin(["time"; "vaccswitch"; "coverage"; "replicate"; "vaccstart"; "duration"], names(sim_pulse))))])
    
    quantiles_pulse[!,:strategy] .= "seasonal"
    
    # Add vaccswitch back to dataset
    quantiles_pulse = @pipe sim_pulse |> 
                select(_,[:time, :coverage, :vaccstart, :duration, :vaccswitch]) |> 
                unique(_) |> 
                outerjoin(_, quantiles_pulse; on = [:time, :coverage, :vaccstart, :duration])

    return quantiles_pulse

end

# Simulate quantiles for pulse/seasonal scenario for each vaccstart in a loop and combine
quantiles_pulse = DataFrame()
for i in ProgressBar(vaccstart)
    out = simulate_pulse(i)
    append!(quantiles_pulse, out)
end

@info "(✓) Done!"


# Combine
verbose && @info "(♻) Combining data..."
sim = vcat(quantiles_pulse, quantiles_cont, cols=:union)
CSV.write("output/vaccsim.csv", sim)
@info "(✓) Done!"


