# Modified Von Mises forcing
function mises(t, φ, η, β₀, k)
    current = exp((cos(2.0*π*(t-φ)/365.0)-1)/k^2)
    min = exp((cos(π)-1)/k^2)
    beta_eff = η * ((current - min)/(1 - min)) + β₀
    return beta_eff
end

# Cosine forcing
function cosine(t, φ, η, β₀, k) # k has no function in the cosine, but we leave it here so we don't have to change the model function
    k = 1.0
    beta_eff = β₀ * (1.0 + η * cos(2.0 *π*(t-φ)/365.0))
    return beta_eff
end
