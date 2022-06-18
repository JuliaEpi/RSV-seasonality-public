"""
    load_data([path]; startdate::Dates.DateTime, enddate::Dates.DateTime)

Return a matrix representing time (first column) and the number of cases (second column).

If `path` isn't specified, data included in the package will be used.
"""
function load_data(;

    version="2022-02-06",
    startdate=Dates.Date(0, 1, 1),
    enddate=Dates.Date(9999, 12, 31))

    # Generate path
    path = joinpath(dirname(pathof(RSVSeasonality)), "..", "data", string("data_NSW_RSV_weekly_", version, ".csv"));
    # Get RSV data
    data = DataFrame(CSV.File(path))
    # Drop observations where there are no case data.
    data = dropmissing(data, :npos)
    # Filter based on dates.
    data = data[startdate .≤ data.date .≤ enddate, :]

    # 2. Attack rate point estimate data 
    data_ar = [33 55; 52 82; 73+38+9 165+147+44] # Munywoki 2015, Table 1 (numerator) and table S1 (denominator)

    return  data.time*1.0, data.date, data.npos, data_ar

end
