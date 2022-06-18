
function load_mobility(; version="2022-02-06")

    # Get mobility data and define interpolation function
    mobility = DataFrame(CSV.File(joinpath(dirname(pathof(RSVSeasonality)), "..", "data", string("mobility_fit_", version, ".csv"))))
    calcmobility = Interpolations.LinearInterpolation(mobility.time, mobility.scaling, extrapolation_bc=Interpolations.Flat())

    return calcmobility

end

# NOTE: Allows us to do `Symbolics.@register`, thus allowing us
# to create a symbolic version of the ODE.
interpolation_evaluator(itr, t) = itr(t)

