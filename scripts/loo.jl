using ArgParse, DrWatson
using StatsPlots
using RSVSeasonality
using DynamicPPL
using Turing, ArviZ
using Distributions, StatsBase, LabelledArrays, MCMCChains, Serialization
using Random
using Dates
using DynamicPPL
using ModelingToolkit

#using LinearAlgebra
#using Setfield

include(scriptsdir("utils.jl"))

s = ArgParseSettings()
@add_arg_table! s begin
    "runs"
    help = "runs to process"
    action = :store_arg
    nargs = '*'
    "--prefix"
    help = "Prefix that name of run should match, e.g. '2021-08-19T14' for a runs started on 19th of August 2021 between 2pm and 3pm. Either this or runs should be provided."
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

# Load the data
@info "(♻) Loading data..."
const times, dates, data_ts, data_ar = RSVSeasonality.load_data();
const N_ar = data_ar[:,2]
const n_seasons = convert(Int, ceil((maximum(times)-minimum(times))/365-2))
const times_ind = convert.(Int, (times ./ 7) .+ 1) # indices of incidence vector for which we have data (some missing data points)
const calcmobility = RSVSeasonality.load_mobility()
@info "(✓) Data loaded!"

# Get posterior predictive
@info "(♻) Calculating posterior predictive..."
model_predict = DynamicPPL.decondition(model, :obs_ts);
posterior_predictive = predict(model_predict, trace);
@info "(✓) Done!"

# Get pointwise loglikelihoods
@info "(♻) Calculating pointwise log likelihoods..."
ll = Turing.pointwise_loglikelihoods(model, MCMCChains.get_sections(trace, :parameters));
# reshape to array
ll_array = permutedims(cat(values(ll)...; dims=3), (2,1,3))
@info "(✓) Done!"

# Generate inferenceData object 
@info "(♻) Make idata objects..."
idata = from_mcmcchains(
                trace;
                posterior_predictive=posterior_predictive,
                log_likelihood=Dict("observations" => ll_array),
                observed_data=Dict("observations" => data_ts),
                library="Turing",
                )
@info "(✓) Done!"

# Compute LOO
@show ArviZ.loo(idata)
