"""
    timespan(ode)

Return default timespan for `ode`.
"""
timespan(ode::ODEProblem) = ode.timespan

"""
    parameternames(ode)

Return the names of the parameters for `ode`.
"""
parameternames(f) = parameternames(typeof(f))
function parameternames(ode::ODEProblem)
    # This just forwards the implementation to function defining the dynamics.
    return parameternames(ode.f.f)
end

"""
    defaultparameters(ode)

Return the default parameters for `ode`.
"""
function defaultparameters(ode::ODEProblem)
    return ode.p
end

"""
    statenames(ode)

Return the names of the states for `ode`.
"""
statenames(f) = statenames(typeof(f))
function statenames(ode::ODEProblem)
    # This just forwards the implementation to function defining the dynamics.
    return statenames(ode.f.f)
end

"""
    initialstates(ode)

Return the default initial state for `ode`.
"""
function initialstates(ode::ODEProblem)
    return ode.u0
end

"""
    solvesettings(ode)

Return the default settings used for `solve` for this `ode`.
"""
function solvesettings(ode::ODEProblem)
    # This just forwards the implementation to function defining the dynamics.
    return solvesettings(ode.f.f)
end

include("differential_equations/seirrs4.jl")
include("differential_equations/seirrs4_init.jl")
include("differential_equations/seirrs4_vacc.jl")
include("differential_equations/seirrs4_vacc_pulse.jl")
