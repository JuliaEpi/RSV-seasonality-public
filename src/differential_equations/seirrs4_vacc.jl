Base.@kwdef struct SEIRRS4Vacc!{B, Int}
    betafunc::B=mises
    simyears::Int=10
end

function DifferentialEquations.ODEProblem(f::SEIRRS4Vacc!; kwargs...)
    u0=initialstates(f)
    tspan=timespan(f)
    p = defaultparameters(f)
    return DifferentialEquations.ODEProblem(f, u0, tspan, p; kwargs...)

end


timespan(f::SEIRRS4Vacc!) = (-7.0, f.simyears*365.0)

function defaultparameters(f::SEIRRS4Vacc!)
    # Create a dummy so we can get the default initial conditions.
    dummy = SEIRRS4!(betafunc=f.betafunc)
    p = defaultparameters(dummy)
    return vcat(p, LVector(v=0.0))
end

function initialstates(f::SEIRRS4Vacc!)
    # Create a dummy so we can get the default initial conditions.
    dummy = SEIRRS4!(betafunc=f.betafunc)
    u0 = initialstates(dummy)
    T = eltype(u0)
    return vcat(u0, LVector(V=zero(T), switch=zero(T))) 
end


function solvesettings(::SEIRRS4Vacc!)
    return (
        abstol = 1.0e-8,
        reltol = 1.0e-8,
        saveat = 1.0,
        save_idxs = [21, 22, 23, 24, 25, 26],
        isoutofdomain = (u,p,t)->any(<(0),u)
    )
end

function (wrapper::SEIRRS4Vacc!)(du, u, p, t)
    # states
    S₀, E₀, I₀, R₀¹, R₀², S₁, E₁, I₁, R₁¹, R₁², S₂, E₂, I₂, R₂¹, R₂², S₃, E₃, I₃, R₃¹, R₃² = u[1:20]
    N = S₀ +  E₀ +  I₀ +  R₀¹ +  R₀² +  S₁ +  E₁ +  I₁ +  R₁¹ +  R₁² +  S₂ +  E₂ +  I₂ +  R₂¹ +  R₂² +  S₃ +  E₃ +  I₃ +  R₃¹ +  R₃²

    # params
    @unpack β₀, η, φ, δ₁, δ₂, δ₃, γ₀, γ₁, γ₂, γ₃, μ, σ, k, v = p
    ω = 1.0 / p.ω

    # FOI
    βeff = wrapper.betafunc(t, φ, η, β₀, k)
    λ₀ = βeff*(I₀ + I₁ + I₂ + I₃)/N
    λ₁ = λ₀*δ₁
    λ₂ = λ₀*δ₁*δ₂
    λ₃ = λ₀*δ₁*δ₂*δ₃

    # vaccination
    switch = u[26] # Switch for vaccination (1=on, 0=off)
    switch == 1.0 ? v=v : v=0.0

    # change in states
    du[1] = (1-v)*μ*N - λ₀*S₀ - μ*S₀ 
    du[2] = λ₀*S₀ - σ*E₀ - μ*E₀
    du[3] = σ*E₀ - γ₀*I₀ - μ*I₀
    du[4] = v*μ*N + γ₀*I₀ - ω*2*R₀¹ - μ*R₀¹
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

    du[25] = v*μ*N # cumulative vaccinated in level 0
    du[26] = 0.0 # Switch for vaccination (1=on, 0=off)

    return nothing
end
