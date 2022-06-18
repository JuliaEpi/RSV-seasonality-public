using Distributions: Distributions
using Distributions: StatsBase

using AdvancedHMC: Random
using AdvancedHMC: AdvancedHMC
using AdvancedHMC: AbstractMCMC
using AdvancedHMC.AbstractMCMC: BangBang
using AdvancedHMC.AbstractMCMC.BangBang.SetfieldImpl: Setfield

using Turing.Inference: AdvancedMH

"""
    MixtureSampler <: AbstractMCMC.AbstractSampler

Represents a mixture of samplers.

$(FIELDS)
"""
struct MixtureSampler{Cs, Ws} <: AbstractMCMC.AbstractSampler
    "component/samplers"
    components::Cs
    "weights of different samplers/components"
    weights::Ws
end

function MixtureSampler(components)
    n = length(components)
    # Using `map` means that we get the same container type as `components`,
    # e.g. if `components isa Tuple` then weights is also a `Tuple`.
    return MixtureSampler(components, fill(1/n, n))
end
Base.getindex(sampler::MixtureSampler, I) = sampler.components[I]
Distributions.components(sampler::MixtureSampler) = sampler.components
Distributions.weights(sampler::MixtureSampler) = sampler.weights

"""
    ManyModels <: AbstractMCMC.AbstractModel

Wraps multiple `models` each which can be extracted using [`getmodels`](@ref).

$(FIELDS)
"""
struct ManyModels{Ms} <: AbstractMCMC.AbstractModel
    models::Ms
end
ManyModels(models::AbstractMCMC.AbstractModel...) = ManyModels(models)

Base.getindex(models::ManyModels, I...) = Base.getindex(models.models, I...)

"""
    getmodel(models, I...)

Return the model corresponding to `models[i]` if `models isa `ManyModels`,
otherwise return `models`.
"""
getmodel(models::ManyModels, I...) = models[I...]
getmodel(model, I...) = model

"""
    MixtureState{T,S}

Represents the state of a [`MixtureSampler`](@ref).

$(FIELDS)
"""
struct MixtureState{T,S}
    "current transition"
    transition::T
    "all the previously used states for the different samplers"
    states::S
    "index of the current state"
    current_index::Int
end
Base.getindex(state::MixtureState, I) = state.states[I]
function BangBang.setindex!!(state::MixtureState, value, indices...)
    states = BangBang.setindex!!(state.states, value, indices...)
    return MixtureState(states, state.current_index)
end

"""
    current(state::MixtureState)

Return the state in `state` which is considered the current state.
"""
current(state::MixtureState) = state[state.current_index]

"""
    setcurrent!!(state::MixtureState, index)

Return new instance of `MixtureState` but now with the current state
represented by the state at index `index`.
"""
function setcurrent!!(state::MixtureState, index)
    return MixtureState(state.states, index)
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    models,
    mixture_sampler::MixtureSampler;
    kwargs...
)
    # TODO: Should we also allow specification of per-sampler kwargs?
    samplers = Distributions.components(mixture_sampler)
    # TODO: Also allow `samplers` to NOT be a `Tuple`.
    initial_transitions_states = ntuple(length(samplers)) do i
        AbstractMCMC.step(rng, getmodel(models, i), samplers[i]; kwargs...)
    end
    initial_transitions = map(first, initial_transitions_states)
    initial_states = map(Base.Fix2(getindex, 2), initial_transitions_states)

    # Sample the first sampler we use.
    sampler_index = rand(rng, Categorical(Distributions.weights(mixture_sampler)))
    initial_transition = initial_transitions[sampler_index]
    initial_state = MixtureState(initial_transition, initial_states, sampler_index)

    return initial_transition, initial_state
end

function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    models,
    mixture_sampler::MixtureSampler,
    mixture_state::MixtureState;
    kwargs...
)
    # Sample the sampler to use in this iteration.
    current_index = rand(rng, Categorical(Distributions.weights(mixture_sampler)))
    current_sampler = mixture_sampler[current_index]

    # Transfer information from previous state to the state to be used
    # in the current iteration.
    previous_state = current(mixture_state)
    current_state = state_from_transition(
        previous_state, mixture_state[current_index]
    )

    # Take a step for the chosen sampler, using the corresponding model.
    transition, state = AbstractMCMC.step(
        rng,
        getmodel(models, current_index),
        current_sampler,
        current_state
    )

    # Update the state.
    mixture_state = BangBang.setindex!!(mixture_state, state, current_index)
    mixture_state = setcurrent!!(mixture_state, current_index)

    return transition, mixture_state
end

"""
    params(state)

Returns the parameters in `state`.
"""
StatsBase.params

"""
    setparams!(state, params)

Return a new instance of `state` with parameters set to `params`.
"""
setparams!!

"""
    state_from_transition(state, transition)

Return a new instance of `state` using necessary information from `transition`.

Default implementation extracts the parameters from `transition` using
[`StatsBase.params`](@ref) and then updates the parameters used in `state` using
[`setparams!!`](@ref), if necessary.
"""
function state_from_transition(state, transition)
    θ = StatsBase.params(transition)
    return setparams!!(state, θ)
end

StatsBase.params(transition::AdvancedMH.Transition) = state.params
function setparams!!(state::AdvancedMH.Transition, θ)
    Setfield.@set state.params = θ
end

StatsBase.params(transition::AdvancedHMC.Transition) = transition.z.θ
function setparams!!(state::AdvancedHMC.HMCState, θ)
    # TODO: Check if we need to update the logp to,
    # for example by re-evaluation.
    Setfield.@set state.transition.z.θ = θ
end
