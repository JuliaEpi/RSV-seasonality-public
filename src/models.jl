"""
    turingmodel(times, times_ind, prior, theta_fix, problem_init, problem, solvsettings; kwargs...)

# Returns
Execution of the model returns a tuple of length 2 where
1. is a `NamedTuple` with keys `sol, p = p, u0 = u0, inc_w, rho, inc_symp, inc_rep, obs_ts, obs_ar1, obs_ar2, obs_ar3,`.
2. is a flag `issuccess` specifying whether any of the variables went outside of their bounds/
   did not satisfy their constraints.

## `issuccess`
The `issuccess` flag might seem a bit surprising, but it can be useful for cases which are very sensitive to
initial conditions. In the case where we force a variable to satisfy some constraints by
usage of `max`/`min`/`clamp`, the gradient information will be lost, and this is a big issue for
samplers which rely on gradient information, e.g. Hamiltonian Monte Carlo (HMC).
_But_ we can still follow the gradient wrt. the rest of the parameter space during adaptation/warmup;
As long as we do not step outside of the bounds during inference, it's a-okay.
Therefore allowing full execution of the model even when constraints are broken can be
used to move from a badly-behaved region of the space to a more well-behaved region,
even if the gradient wrt. some of the parameters is 0.

Once we have a chain from this model, we can use `Turing.generated_quantities` to recover the
values returned by the model, e.g. `issuccess`, and thus check whether the bounds were actually
broken during sampling. If not, we're good! If they were, well then we need to address by changing the prior.

# Arguments
- `times`: the times corresponding to `observations`, i.e. `length(times) == length(observations)`.
- `times_ind`: the index of all time points for which there is data
- `prior`: the prior `Turing.Model`. Expected to return a `LArray` with the parameters.
- `theta_fix`: the fixed part of the parameters, i.e. not to be sampled.
- `problem_init`: the `ODEProblem` to solve for the run-in phase to reach table periodic orbit
- `problem`: the `ODEProblem` to solve for the period where we have data (in quasi equilibrium)
- `solvsettings`: solver settings.

# Keyword arguments
- `cb_terminate`: whether to use a callback function to derminate the init model at equilibrium (defaults to nothing)
- `saveat_init`: either a scalar -7 to denote at which time step to save the init model (without the callback), or a 1D array [730] at which time step to save the init model (with callback)
- `solver`: the solver used in the call to `solve`. 
- `N_ar`: the denominator for the attack rate data. 
- `n_seasons`: the number of pre-pandemic seasons to fit the AR data to.
- `rho_func`: the function of the reporting rate. Defaults to gompertz.
- `sensealg=ForwardDiffSensitivity()`: the sensitivity algorithm passed to `solve`.
  See [https://diffeqflux.sciml.ai/stable/ControllingAdjoints/](https://diffeqflux.sciml.ai/stable/ControllingAdjoints/)
  for more info on how to choose this keyword.

# Variables
- `obs_ts ~ arraydist(@. NegativeBinomial2(ψ, increp))`
- obs_ar1 ~ arraydist(@. Binomial(N_ar[1], ar[:,1]))
- obs_ar2 ~ arraydist(@. Binomial(N_ar[2], ar[:,2]))
- obs_ar3 ~ arraydist(@. Binomial(N_ar[3], ar[:,3]))
- Whatever variables are present in the `prior` model.
"""
@model function turingmodel(
    times,
    times_ind,
    prior,
    theta_fix,
    problem_init,
    problem,
    solvsettings;
    cb_terminate,
    saveat_init,
    solver,
    N_ar,
    n_seasons,
    rho_func,
    sensealg=ForwardDiffSensitivity(),
)
    # Used to keep track of whether we're doing something bad.
    issuccess = true

    # Sample prior parameters.
    theta_est = @submodel prior()

    # Update `p`.
    p = vcat(theta_est, theta_fix) # `vcat` on two `LArray`s gives back a `LArray`
    
    # Run model with inits to stable periodic orbit
    problem_init_new = remake(problem_init; p=p) 
    sol_init = solve(problem_init_new, 
                    solver;
                    sensealg=sensealg, 
                    callback=cb_terminate,
                    saveat=saveat_init,
                    reltol=solvsettings.reltol,
                    abstol=solvsettings.abstol,
                    isoutofdomain = (u,p,t)->any(<(0),u));


    # Return early if integration failed
    issuccess &= (sol_init.retcode === :Success || sol_init.retcode === :Terminated)
    if !issuccess
        Turing.@addlogprob! -Inf #convert(eltype(p), -Inf) #-Inf
        return nothing
    end

    # If successful, make array with new inits at stable periodic orbit     
    u0 = sol_init[end] 
                     
    # Remake problem, solve ODE and return incidence.
    problem_new = remake(problem; u0=u0, p=p)
    sol = solve(problem_new,
                solver;
                sensealg=sensealg,
                solvsettings...)

    # Switch `issuccess` if integration failed.
    # This is equvalent to:
    #
    #    issuccesss = issuccesss && (sol.retcode === :Success)
    #
    # i.e. BOTH need to hold `true` for `issuccess` to be `true`
    # after this statement.
    issuccess &= (sol.retcode === :Success)

    # Return early if integration failed
    if !issuccess
        Turing.@addlogprob! -Inf #convert(eltype(p), -Inf) #-Inf
        return (; sol, p = p), false
    end

    # Convert solution to array
    sol_array = Array(sol')
    
    # Compute daily incidence.
    inc_w = diff(sol_array[:,end-3:end], dims=1)

    # Avoid negative incidence due to numerical instability issue.
    inc_w = max.(eltype(inc_w)(1e-9), inc_w)

    # Compute symptomatic incidence.
    inc_symp = inc_w .* [(1-p.p₀) (1-p.p₁) (1-p.p₂) (1-p.p₃)] 

    # Instantaneous reporting rate.
    rho = rho_func.(p.ρ₀, p.q .* times)

    # Switch `issuccess` to `false` if we're not within the wanted bounds. 
    # Note: this is only relevant when using the exponential reporting function
    issuccess &= all(0.0 .≤ rho .≤ 1.0)

    # NOTE: We still continue even if `issuccess` is false because
    # using the gradient information from this imperfect solution,
    # we can move towards more well-behaved regions of the parameter space.
    # But this should all be checked after inference.
    rho = min.(one(eltype(rho)), rho) # reporting rate cannot be larger than 1

    # Subset to timepoints where there is data (for calculation of likelihood).
    inc_obs = inc_symp[times_ind,:]

    # Compute total reported incidence from levels 0-2.
    inc_rep = vec(sum(inc_obs[:,1:3], dims=2)) .* rho

    # Calculate the cumulative incidences in each level at the end of the year
    cumul = inc_w[1:n_seasons*52,:] # select pre-pandemic data (n seasons) 
    cumul = sum(cumul, dims=1) ./ n_seasons # calculate average cumulative cases per year
    # For Munywoki data: sum level 1 and 2
    cumul = hcat(cumul[:,1], sum(cumul[:,2:3], dims=2), cumul[:,4])

    # Calculate the average number of individuals in each level at the end the year
    pop = Array(sol')[(sol.t .== 357.0), 1:end-4]
    pop = reshape(pop, div(size(pop, 2), 4), 4)
    pop = sum(pop, dims=1)
    # For Munywoki data: sum level 1 and 2
    pop = hcat(pop[:,1], sum(pop[:,2:3], dims=2), pop[:,4])

    # Compute yearly attack rates by level.
    # numerator: cumulative cases at the end of the year for level i 
    # denominator: N of individuals at the end of the year for level i
    ar = cumul ./ pop

    # replace values > 1.0 with 1.0 (attach rate cannot be larger than 1)
    ar = min.(one(eltype(ar)), ar)
    
    # repeat for 10 seasons (Loglikelihood is evaluated each year before the pandemic)
    ar = repeat(ar, outer=n_seasons)

    # Likelihood function.
    # HACK(t-tfjelde): makes it so that when we compute pointwise log-likelihood,
    # we treat the observations as a collection of independent observations rather
    # than a single multivariate observation.
    if __context__ isa DynamicPPL.PointwiseLikelihoodContext
        obs_ts = Array{Int}(undef, size(inc_rep))
        for i in eachindex(inc_rep)
            obs_ts[i] ~ NegativeBinomial2(p.ψ, inc_rep[i])
        end

        n = size(ar, 1)
        obs_ar1 = Vector{Int}(undef, n)
        obs_ar2 = Vector{Int}(undef, n)
        obs_ar3 = Vector{Int}(undef, n)
        for i = 1:n
            obs_ar1[i] ~ Binomial(N_ar[1], ar[i,1])
            obs_ar2[i] ~ Binomial(N_ar[2], ar[i,2])
            obs_ar3[i] ~ Binomial(N_ar[3], ar[i,3])
        end
    else
        obs_ts ~ arraydist(@. NegativeBinomial2(p.ψ, inc_rep))
        obs_ar1 ~ arraydist(@. Binomial(N_ar[1], ar[:,1]))
        obs_ar2 ~ arraydist(@. Binomial(N_ar[2], ar[:,2]))
        obs_ar3 ~ arraydist(@. Binomial(N_ar[3], ar[:,3]))
    end
    # Let's indicate in the return-value whether we succeeded or not
    # at running the full model.
    return (; sol, 
            p = p, 
            u0 = u0,
            inc_w, 
            rho,
            inc_symp,
            inc_rep, 
            obs_ts, 
            obs_ar1,
            obs_ar2,            
            obs_ar3,
            ), issuccess

end

@model function turingmodel_noar(
    times,
    times_ind,
    prior,
    theta_fix,
    problem_init,
    problem,
    solvsettings;
    cb_terminate,
    saveat_init,
    solver,
    rho_func,
    sensealg=ForwardDiffSensitivity(),
)
    # Used to keep track of whether we're doing something bad.
    issuccess = true

    # Sample prior parameters.
    theta_est = @submodel prior()

    # Update `p`.
    p = vcat(theta_est, theta_fix) # `vcat` on two `LArray`s gives back a `LArray`
    
    # Run model with inits to stable periodic orbit
    problem_init_new = remake(problem_init; p=p) 
    sol_init = solve(problem_init_new, 
                    solver;
                    sensealg=sensealg, 
                    callback=cb_terminate,
                    saveat=saveat_init,
                    reltol=solvsettings.reltol,
                    abstol=solvsettings.abstol,
                    isoutofdomain = (u,p,t)->any(<(0),u));


    # Return early if integration failed
    issuccess &= (sol_init.retcode === :Success || sol_init.retcode === :Terminated)
    if !issuccess
        Turing.@addlogprob! -Inf #convert(eltype(p), -Inf) #-Inf
        return nothing
    end

    # If successful, make array with new inits at stable periodic orbit     
    u0 = sol_init[end] 
                     
    # Remake problem, solve ODE and return incidence.
    problem_new = remake(problem; u0=u0, p=p)
    sol = solve(problem_new,
                solver;
                sensealg=sensealg,
                solvsettings...)

    # Switch `issuccess` if integration failed.
    # This is equvalent to:
    #
    #    issuccesss = issuccesss && (sol.retcode === :Success)
    #
    # i.e. BOTH need to hold `true` for `issuccess` to be `true`
    # after this statement.
    issuccess &= (sol.retcode === :Success)

    # Return early if integration failed
    if !issuccess
        Turing.@addlogprob! -Inf #convert(eltype(p), -Inf) #-Inf
        return (; sol, p = p), false
    end

    # Convert solution to array
    sol_array = Array(sol')
    
    # Compute daily incidence.
    inc_w = diff(sol_array[:,end-3:end], dims=1)

    # Avoid negative incidence due to numerical instability issue.
    inc_w = max.(eltype(inc_w)(1e-9), inc_w)

    # Compute symptomatic incidence.
    inc_symp = inc_w .* [(1-p.p₀) (1-p.p₁) (1-p.p₂) (1-p.p₃)] 

    # Instantaneous reporting rate.
    rho = rho_func.(p.ρ₀, p.q .* times)

    # Switch `issuccess` to `false` if we're not within the wanted bounds. 
    # Note: this is only relevant when using the exponential reporting function
    issuccess &= all(0.0 .≤ rho .≤ 1.0)

    # NOTE: We still continue even if `issuccess` is false because
    # using the gradient information from this imperfect solution,
    # we can move towards more well-behaved regions of the parameter space.
    # But this should all be checked after inference.
    rho = min.(one(eltype(rho)), rho) # reporting rate cannot be larger than 1

    # Subset to timepoints where there is data (for calculation of likelihood).
    inc_obs = inc_symp[times_ind,:]

    # Compute total reported incidence from levels 0-2.
    inc_rep = vec(sum(inc_obs[:,1:3], dims=2)) .* rho

    # Likelihood function.
    # HACK(t-tfjelde): makes it so that when we compute pointwise log-likelihood,
    # we treat the observations as a collection of independent observations rather
    # than a single multivariate observation.
    if __context__ isa DynamicPPL.PointwiseLikelihoodContext
        obs_ts = Array{Int}(undef, size(inc_rep))
        for i in eachindex(inc_rep)
            obs_ts[i] ~ NegativeBinomial2(p.ψ, inc_rep[i])
        end
    else
        obs_ts ~ arraydist(@. NegativeBinomial2(p.ψ, inc_rep))
    end
    # Let's indicate in the return-value whether we succeeded or not
    # at running the full model.
    return (; sol, 
            p = p, 
            u0 = u0,
            inc_w, 
            rho,
            inc_symp,
            inc_rep, 
            obs_ts, 
            ), issuccess

end
