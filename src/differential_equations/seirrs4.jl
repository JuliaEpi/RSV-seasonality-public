# ODE model
Base.@kwdef struct SEIRRS4!{C,B}
    calcmobility::C=nothing #calcmobility
    betafunc::B=mises
end

function DifferentialEquations.ODEProblem(
    f::SEIRRS4!;
    u0=initialstates(f),
    tspan=timespan(f),
    p=defaultparameters(f),
    kwargs...
)
    return DifferentialEquations.ODEProblem(f, u0, tspan, p; kwargs...)
end

timespan(::SEIRRS4!) = (-7.0, 4109.0)

parameternames(::Type{<:SEIRRS4!}) = (:ρ₀, :ψ, :q, :β₀, :η, :ω, :φ, :k, :δ₁, :δ₂, :δ₃, :γ₀, :γ₁, :γ₂, :γ₃, :μ, :σ, :p₀, :p₁, :p₂, :p₃)

function defaultparameters(f::SEIRRS4!)
    return @LArray [
        0.004, 0.01, 0.000857155, 3.5, 0.5, 250.0, 180.0, 1.0, 
        0.89, 0.81, 0.33, 1.0/9.0, 1.0/9.0, 1.0/3.9, 1.0/1.6, 1.0/(80.0*365.0), 1.0/4.0, 
        0.091, 0.173, 0.173, 0.778
    ] parameternames(f)
end

statenames(::Type{<:SEIRRS4!}) = (:S₀, :E₀, :I₀, :R₀¹, :R₀², :S₁, :E₁, :I₁, :R₁¹, :R₁², :S₂, :E₂, :I₂, :R₂¹, :R₂², :S₃, :E₃, :I₃, :R₃¹, :R₃², :C₀, :C₁, :C₂, :C₃)
function initialstates(f::SEIRRS4!)
    return @LArray [
        8163900.0, 100.0, 0.0, 0.0, 0.0, # Level 0
        0.0, 0.0, 0.0, 0.0, 0.0, # Level 1
        0.0, 0.0, 0.0, 0.0, 0.0, # Level 2
        0.0, 0.0, 0.0, 0.0, 0.0, # Level 3
        0.0, 0.0, 0.0, 0.0 # Cumulative incidence
    ] statenames(f)
end

function solvesettings(::SEIRRS4!)
    return (
        abstol = 1.0e-8,
        reltol = 1.0e-8,
        isoutofdomain = (u,p,t)->any(<(0),u),
        #isoutofdomain = (u,p,t)->any(x->x<(0),u),
        saveat = 7.0
    )
end

function (wrapper::SEIRRS4!)(du, u, p, t)

    # states
    S₀, E₀, I₀, R₀¹, R₀², S₁, E₁, I₁, R₁¹, R₁², S₂, E₂, I₂, R₂¹, R₂², S₃, E₃, I₃, R₃¹, R₃² = u[1:20]
    N = S₀ +  E₀ +  I₀ +  R₀¹ +  R₀² +  S₁ +  E₁ +  I₁ +  R₁¹ +  R₁² +  S₂ +  E₂ +  I₂ +  R₂¹ +  R₂² +  S₃ +  E₃ +  I₃ +  R₃¹ +  R₃²

    # params
    @unpack β₀, η, φ, δ₁, δ₂, δ₃, γ₀, γ₁, γ₂, γ₃, μ, σ, k = p
    ω = 1.0 / p.ω

    # google mobility reduction
    i = interpolation_evaluator(wrapper.calcmobility, t)  

    # FOI
    βeff = wrapper.betafunc(t, φ, η, β₀, k) * i
    λ₀ = βeff*(I₀ + I₁ + I₂ + I₃)/N
    λ₁ = λ₀*δ₁
    λ₂ = λ₀*δ₁*δ₂
    λ₃ = λ₀*δ₁*δ₂*δ₃

    # change in states
    @inbounds begin
        du[1] = μ*N - λ₀*S₀ - μ*S₀
        du[2] = λ₀*S₀ - σ*E₀ - μ*E₀
        du[3] = σ*E₀ - γ₀*I₀ - μ*I₀
        du[4] = γ₀*I₀ - ω*2*R₀¹ - μ*R₀¹
        du[5] = ω*2*R₀¹ - ω*2*R₀² - μ*R₀²
    
        du[6] = ω*2*R₀² - λ₁*S₁ - μ*S₁
        du[7] = λ₁*S₁ - σ*E₁ - μ*E₁
        du[8] = σ*E₁ - γ₁*I₁ - μ*I₁
        du[9] = γ₁*I₁ - ω*2*R₁¹ - μ*R₁¹
        du[10] = ω*2*R₁¹ - ω*2*R₁² - μ*R₁²
    
        du[11] = ω*2*R₁² - λ₂*S₂ - μ*S₂
        du[12] = λ₂*S₂ - σ*E₂ - μ*E₂
        du[13] = σ*E₂ - γ₂*I₂ - μ*I₂
        du[14] = γ₂*I₂ - ω*2*R₂¹ - μ*R₂¹
        du[15] = ω*2*R₂¹ - ω*2*R₂² - μ*R₂²
    
        du[16] = ω*2*R₂² + ω*2*R₃² - λ₃*S₃ - μ*S₃
        du[17] = λ₃*S₃ - σ*E₃ - μ*E₃
        du[18] = σ*E₃ - γ₃*I₃ - μ*I₃
        du[19] = γ₃*I₃ - ω*2*R₃¹ - μ*R₃¹
        du[20] = ω*2*R₃¹ - ω*2*R₃² - μ*R₃²
    
        du[21] = σ*E₀ # cumulative incidence level 0
        du[22] = σ*E₁ # cumulative incidence level 1
        du[23] = σ*E₂ # cumulative incidence level 2
        du[24] = σ*E₃ # cumulative incidence level 3
    end

    return nothing
end
