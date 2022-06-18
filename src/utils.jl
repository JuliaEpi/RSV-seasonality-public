"""
    replace_args(model::DynamicPPL.Model; kwargs...)

Return new `Model` with arguments specified in `kwargs` replaced in `model`.
"""
replace_args(model::DynamicPPL.Model; kwargs...) = replace_args(model, NamedTuple(kwargs))
function replace_args(model::DynamicPPL.Model, args)
    return DynamicPPL.Model{DynamicPPL.getmissings(model)}(
        model.name, model.f, merge(model.args, args), model.defaults
    )
end

"""
    setmissings(model::DynamicPPL.Model, syms::Symbol...)

Return new `Model` with `syms` set as `missing`.
"""
function setmissings(model::DynamicPPL.Model, syms::Symbol...)
    return DynamicPPL.Model{syms}(model.name, model.f, model.args, model.defaults)
end


"""
    fixedparameters(problem::ODEProblem, prior::Model)

Return `LArray` of parameters from `problem` considered fixed by `prior`.

It is assumed that
- `problem.p` holds an example of all the parameters, and
- `prior()` returns an iterable supporting `keys(...)`, e.g. `NamedTuple`,
  representing the variables which are sampled.

# Examples
```julia
julia> using RSVSeasonality, DifferentialEquations

julia> # Set up the problem.
       problem = ODEProblem(SEIRRS4!());

julia> # Let's see what parameters are present in `problem`.
       keys(problem.p)
       (:ρ₀, :ψ, :q, :β₀, :η, :ω, :φ, :k, :δ₁, :δ₂, :δ₃, :γ₀, :γ₁, :γ₂, :γ₃, :μ, :σ, :p₀, :p₁, :p₂, :p₃)

julia> # And what is sampled from the prior
       keys(RSVSeasonality.prior()())
(:ρ₀, :ψ, :q, :β₀, :η, :ω, :φ, :k, :δ₁, :δ₂, :δ₃)

julia> # Hence `fixedparameters` returns those present in `problem` but _not_ in `prior`:
       RSVSeasonality.fixedparameters(problem, RSVSeasonality.prior())

```
"""

function fixedparameters(problem::ODEProblem, prior::Turing.Model)
    θ_free = prior()
    θ_fixed_labels = ((k for k in keys(problem.p) if k ∉ keys(θ_free))..., )
    return @LArray problem.p[collect(θ_fixed_labels)] θ_fixed_labels
end

