Base.@kwdef struct SEIRRS4Vacc_pulse!{B, Int}
    betafunc::B=mises
    simyears::Int=10
end

function DifferentialEquations.ODEProblem(f::SEIRRS4Vacc_pulse!; kwargs...)
    u0=initialstates(f)
    tspan=timespan(f)
    p = defaultparameters(f)
    return DifferentialEquations.ODEProblem(f, u0, tspan, p; kwargs...)

end

timespan(f::SEIRRS4Vacc_pulse!) = (-7.0, f.simyears*365.0)

function defaultparameters(f::SEIRRS4Vacc_pulse!)
    # Create a dummy so we can get the default initial conditions.
    dummy = SEIRRS4!(betafunc=f.betafunc)
    p = defaultparameters(dummy)
    return vcat(p, LVector(v=0.0, d=210.0))
end

function initialstates(f::SEIRRS4Vacc_pulse!)
    return initialstates(SEIRRS4Vacc!(betafunc=f.betafunc))
end

function solvesettings(f::SEIRRS4Vacc_pulse!)
    return solvesettings(SEIRRS4Vacc!(betafunc=f.betafunc))
end


function (wrapper::SEIRRS4Vacc_pulse!)(du, u, p, t)
    # states
    S₀, E₀, I₀, R₀¹, R₀², S₁, E₁, I₁, R₁¹, R₁², S₂, E₂, I₂, R₂¹, R₂², S₃, E₃, I₃, R₃¹, R₃² = u[1:20]
    N = S₀ +  E₀ +  I₀ +  R₀¹ +  R₀² +  S₁ +  E₁ +  I₁ +  R₁¹ +  R₁² +  S₂ +  E₂ +  I₂ +  R₂¹ +  R₂² +  S₃ +  E₃ +  I₃ +  R₃¹ +  R₃²

    # params
    @unpack β₀, η, φ, δ₁, δ₂, δ₃, γ₀, γ₁, γ₂, γ₃, μ, σ, k, v, d = p
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
    # convert coverage and duration to vaccination rate
    ε = (log(1-v)/-d)*switch

    # change in states
    du[1] = (1-v)*μ*N - λ₀*S₀ - μ*S₀ - ε*S₀
    du[2] = λ₀*S₀ - σ*E₀ - μ*E₀ - ε*E₀
    du[3] = σ*E₀ - γ₀*I₀ - μ*I₀ - ε*I₀
    du[4] = v*μ*N + γ₀*I₀ + ε*(S₀ + E₀ + I₀) - ω*2*R₀¹ - μ*R₀¹
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

    du[25] = v*μ*N + ε*(S₀ + E₀ + I₀) # cumulative vaccinated in level 0
    du[26] = 0.0 # Switch for vaccination (1=on, 0=off)

    return nothing
end
