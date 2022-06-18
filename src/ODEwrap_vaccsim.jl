# Returns a Dataframe of all incident cases (total and each level) for vaccsim
function ODEwrap_vaccsim(problem, p, u0, tspan, cb, solvsettings)
    
    # remake problem and solve
    problem_new = remake(problem; p=p, u0=u0, tspan=tspan)
    sol = solve(problem_new, 
                solvsettings.solver, 
                abstol=solvsettings.abstol, 
                reltol=solvsettings.reltol, 
                isoutofdomain=(u,p,t)->any(<(0),u),
                callback=cb,
                saveat=solvsettings.saveat)    
    
    # calculate incidence
    out = diff(Array(sol')[:,end-5:end-2], dims=1)
    # Calculate total incidence
    out = DataFrame(out, :auto)
    rename!(out, 1:4 .=> ["inc"] .* string.(0:3))
    # calculate new vaccinated cases
    out[!,:nvacc] = diff(Array(sol)[end-1,:])
    # add vaccswitch
    out[!,:vaccswitch] = Array(sol)[end,2:end]
    # add time
    insertcols!(out, 1, :time => sol.t[2:end])
    
    return out

end

