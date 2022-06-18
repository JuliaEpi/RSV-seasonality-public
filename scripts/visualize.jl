using ArgParse, DrWatson
using Setfield
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
    help = "runs to visualize"
    action = :store_arg
    nargs = '*'
    "--prefix"
    help = "Prefix that name of run should match, e.g. '2021-08-19T14' for a runs started on 19th of August 2021 between 2pm and 3pm. Either this or runs should be provided."
    # Plot/figure related options.
    "--filename"
    help = "File name of resulting visualization."
    default = "summary.pdf"
    "--theme"
    help = "Plotting theme to use."
    default = :default
    arg_type = Symbol
    # Alternatives.
    "--list-runs"
    help = "If specified, available runs will be listed."
    action = :store_true
    # General.
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
model, trace, run_args = trace_merge(runs)
verbose && @info "(✓) Done!"

# Plotting
pyplot()
const _theme = args["theme"]
verbose && @info "Using plotting theme $(_theme)."
theme(_theme)

const filename = args["filename"]

# Quantiles plotting function
plot_quantiles(args...; kwargs...) = plot_quantiles!(plot(), args...; kwargs...)
function plot_quantiles!(p::Plots.Plot, xs, q=0.95; kwargs...)
    Δ = (1 - q) / 2
    quantiles = mapreduce(hcat, xs) do x
        quantile(x, [Δ, 0.5, 1 - Δ])
    end

    plot!(
        p,
        quantiles[2, :];
        ribbon=(quantiles[2, :] - quantiles[1, :], quantiles[3, :] - quantiles[2, :]),
        label="$(q * 100)% credible intervals",
        kwargs...
    )
    return p
end

# Display the summary of the trace.
show(stdout, MIME"text/plain"(), trace)


# Load the data
verbose && @info "(♻) Loading data..."
version = run_args["version"]
const data_ts = RSVSeasonality.load_data(; version=version)[3];
                                                                #startdate=Dates.Date(run_args["startyear"],1,1),
                                                                #enddate=Dates.Date(run_args["endyear"],12,31));
#const calcmobility = RSVSeasonality.load_mobility(; version=version);
verbose && @info "(✓) Data loaded!"

# Trace plots.
p1 = plot(trace)

internals = filter(
    !=(:max_hamiltonian_energy_error),
    trace.name_map.internals
)
p_internals = plot(MCMCChains.get_sections(trace, :internals)[internals])

# Posterior predictive
verbose && @info "(♻) Predicting..."

# `decondition` model so that `observations` are now treated as random variables.
predictive_model = DynamicPPL.decondition(model, :obs_ts)

# Predict
num_tries = 1
max_tries = 20
predictions = nothing
while isnothing(predictions) && num_tries ≤ max_tries
    try
        global predictions
        predictions = predict(predictive_model, trace)
    catch e
        global num_tries
        num_tries +=1
    end
end
predictions_per_chain = [Array(predictions[:, :, i]) for i = 1:size(predictions, 3)]
verbose && @info "(✓) Predictions done!"

p2 = plot(; legend=:topleft)
for preds in predictions_per_chain
    plot_quantiles!(p2, eachcol(preds), 0.95)
end
# plot!(p2, model.args.observations, color=:black, label="observations")
#observations = DynamicPPL.conditioned(model).obs_ts
plot!(p2, data_ts, color=:black, label="observations")
title!(p2, "Posterior predictive")


# Seasonal beta.
verbose && @info "(♻) Computing posterior of seasonal β..."
betafunc = model.args.problem.f.f.betafunc

φs = vec(trace[:φ])
ηs = vec(trace[:η])
β₀s = vec(trace[:β₀])
ks = if :k in keys(trace)
    vec(trace[:k])
else
    model.args.problem.p.k
end

seasonal_betas_weekly = map(1:52) do week
    betafunc.(7 * week, φs, ηs, β₀s, ks)
end
p3 = plot_quantiles(seasonal_betas_weekly)
title!(p3, "Seasonal β")

verbose && @info "(✓) Seasonal β computed!"

# Generated quantities.
verbose && @info "(♻) Computing generated quantities..."
quantities = generated_quantities(model, MCMCChains.get_sections(trace, :parameters));

# Reporting rate.
rhos = map(quantities) do qty
    qty[1].rho
end;
rhos = vec(rhos); # combine chains into a single vector
rhos = transpose(reduce(hcat, rhos))
p4 = plot_quantiles(eachcol(rhos))
title!(p4, "Reporting proportion ρ")

verbose && @info "(✓) Reporting proportion computed!"

p_full = plot(
    p_internals, p1, plot(p2, p3, p4, layout=(3, 1)),
    layout=(1, 3),
    size=(3000, 1600)
)


@info "Saving plots to $(plotsdir(filename))."

save(plotsdir(filename), p_full)

# Additional information.
@info "(♻) Let's try to also plot some diagnostics using ArviZ.jl..."
using ArviZ
using ArviZ: PyCall
# HACK: Necessary to save the plots from ArviZ.jl
plt = PyCall.pyimport("matplotlib.pyplot")

idata = from_mcmcchains(
    trace;
    library="Turing"
);

plot_energy(idata)
base_filename, ext = Base.Filesystem.splitext(filename)
plt.savefig(plotsdir("$(base_filename)_energy.$(ext)"))
@info "(✓) Finito!"