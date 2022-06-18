module RSVSeasonality

using Interpolations: Interpolations
using DifferentialEquations, DiffEqSensitivity, DiffEqCallbacks
using LabelledArrays
using Turing
using UnPack
using DataFrames, CSV
using StatsFuns
using DocStringExtensions
using Dates


# This defines what should be made available in the scope where
# we do `using RSVSeasonality`, for example since `SEIRS4!`
# is exported, after `using RSVSeasonality` we can access it by just
# doing `SEIRS4!` rather than `RSVSeasonality.SEIRS4!`.
export
# need to rename the models before publishing the repo
    SEIRRS4!,
    SEIRRS4init!,
    SEIRRS4Vacc!,
    SEIRRS4Vacc_pulse!,
    turingmodel,
    turingmodel_noar,
    mises,
    cosine,
    calcmobility,
    gompertz,
    scaled_exp,
    NegativeBinomial2,
    ODEwrap_vaccsim,
    replace_args,
    setmissings

include("utils.jl")
include("foi.jl")
include("reporting_functions.jl")
include("distributions.jl")
include("data.jl")
include("mobility.jl")
include("differential_equations.jl")
include("models.jl")
include("priors.jl")
include("ODEwrap_vaccsim.jl")
include("callbacks.jl")
include("serial_tempering.jl")
include("mixture_sampler.jl")

end
