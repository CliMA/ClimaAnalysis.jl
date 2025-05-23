# Customize these variables to allow other names
LONGITUDE_NAMES = ["long", "lon", "longitude"]
LATITUDE_NAMES = ["lat", "latitude"]
TIME_NAMES = ["t", "time", "valid_time"]
DATE_NAMES = ["date"]
ALTITUDE_NAMES = ["z", "z_reference", "z_physical"]
PRESSURE_NAMES = ["pfull", "pressure_level"]

export times,
    dates,
    longitudes,
    latitudes,
    altitudes,
    pressures,
    time_name,
    date_name,
    longitude_name,
    latitude_name,
    altitude_name,
    pressure_name,
    has_time,
    has_date,
    has_longitude,
    has_latitude,
    has_altitude,
    has_pressure,
    conventional_dim_name,
    find_corresponding_dim_name

"""
    _dim_name(dim_names, allowed_names)

Return the `dim_name` within `dim_names` that is contained in `allowed_names`.
Return `nothing` is not available.

Example
==========

```jldoctest
julia> ClimaAnalysis.Var._dim_name(["z", "lat", "lon"], ["lon", "long"])
1
```
"""
function _dim_name(dim_names, allowed_names)
    return findfirst(possible_name -> possible_name in dim_names, allowed_names)
end

"""
    has_time(var::OutputVar)

Return whether `var` has a `time` dimesnion.
"""
has_time(var::OutputVar) = !isnothing(_dim_name(keys(var.dims), TIME_NAMES))

"""
    has_date(var::OutputVar)

Return whether `var` has a `date` dimesnion.
"""
has_date(var::OutputVar) = !isnothing(_dim_name(keys(var.dims), DATE_NAMES))

"""
    has_longitude(var::OutputVar)

Return whether `var` has a `longitude` dimesnion.
"""
has_longitude(var::OutputVar) =
    !isnothing(_dim_name(keys(var.dims), LONGITUDE_NAMES))

"""
    has_latitude(var::OutputVar)

Return whether `var` has a `latitude` dimesnion.
"""
has_latitude(var::OutputVar) =
    !isnothing(_dim_name(keys(var.dims), LATITUDE_NAMES))

"""
    has_altitude(var::OutputVar)

Return whether `var` has a `altitude` dimesnion.
"""
has_altitude(var::OutputVar) =
    !isnothing(_dim_name(keys(var.dims), ALTITUDE_NAMES))

"""
    has_pressure(var::OutputVar)

Return whether `var` has a `pressure` dimension.
"""
has_pressure(var::OutputVar) =
    !isnothing(_dim_name(keys(var.dims), PRESSURE_NAMES))

"""
    find_dim_name(dim_names::Iterable, allowed_names::Iterable)

Find the first occurrence of a name in `allowed_names` in `dim_names`.

Example
==========

```jldoctest
julia> ClimaAnalysis.Var.find_dim_name(["z", "lat", "lon"], ["lon", "long"])
"lon"
```
"""
function find_dim_name(dim_names, allowed_names)
    dim_name = _dim_name(dim_names, allowed_names)
    isnothing(dim_name) &&
        error("var does not have $(first(allowed_names)) among its dimensions")
    return allowed_names[dim_name]
end

"""
    time_name(var::OutputVar)

Return the name of the `time` dimension in `var`.
"""
time_name(var::OutputVar) = find_dim_name(keys(var.dims), TIME_NAMES)

"""
    times(var::OutputVar)

Return the `time` dimension in `var`.
"""
times(var::OutputVar) = var.dims[time_name(var)]

"""
    date_name(var::OutputVar)

Return the name of the `date` dimension in `var`.
"""
date_name(var::OutputVar) = find_dim_name(keys(var.dims), DATE_NAMES)

"""
    dates(var::OutputVar)

Return the `date` dimension in `var`.

If `dates` is a dimension, return that.

If not, try computing the dates from the `start_date` and the `times`. In this, we assume
`times` is seconds.
"""
function dates(var::OutputVar)
    has_date(var) && return var.dims[date_name(var)]

    # Assuming time is in seconds
    if has_time(var) && haskey(var.attributes, "start_date")
        return time_to_date.(
            Dates.DateTime(var.attributes["start_date"]),
            times(var),
        )
    end

    error(
        "`var` does not have dates among its dimension or start_date among its attributes",
    )
end

"""
    longitude_name(var::OutputVar)

Return the name of the `longitude` dimension in `var`.
"""
longitude_name(var::OutputVar) = find_dim_name(keys(var.dims), LONGITUDE_NAMES)

"""
    longitudes(var::OutputVar)

Return the `longitude` dimension in `var`.
"""
longitudes(var::OutputVar) = var.dims[longitude_name(var)]

"""
    latitude_name(var::OutputVar)

Return the name of the `latitude` dimension in `var`.
"""
latitude_name(var::OutputVar) = find_dim_name(keys(var.dims), LATITUDE_NAMES)

"""
    latitudes(var::OutputVar)

Return the `latitude` dimension in `var`.
"""
latitudes(var::OutputVar) = var.dims[latitude_name(var)]

"""
    altitude_name(var::OutputVar)

Return the name of the `altitude` dimension in `var`.
"""
altitude_name(var::OutputVar) = find_dim_name(keys(var.dims), ALTITUDE_NAMES)

"""
    altitudes(var::OutputVar)

Return the `altitude` dimension in `var`.
"""
altitudes(var::OutputVar) = var.dims[altitude_name(var)]

"""
    pressure_name(var::OutputVar)

Return the name of the `pressure` dimension in `var`.
"""
pressure_name(var::OutputVar) = find_dim_name(keys(var.dims), PRESSURE_NAMES)

"""
    pressures(var::OutputVar)

Return the `pressure` dimension in `var`.
"""
pressures(var::OutputVar) = var.dims[pressure_name(var)]

"""
    conventional_dim_name(dim_name::AbstractString)

Return the type of dimension as a string from longitude, latitude, time, date, altitude, or
pressure if possible or `dim_name` as a string otherwise.
"""
function conventional_dim_name(dim_name::AbstractString)
    dim_name in LONGITUDE_NAMES && return "longitude"
    dim_name in LATITUDE_NAMES && return "latitude"
    dim_name in TIME_NAMES && return "time"
    dim_name in DATE_NAMES && return "date"
    dim_name in ALTITUDE_NAMES && return "altitude"
    dim_name in PRESSURE_NAMES && return "pressure"
    return dim_name
end

"""
    find_corresponding_dim_name(dim_name::AbstractString, dim_names::Iterable)

Find the corresponding dimension name in `dim_names` that matches `dim_name`.

Two names for a dimension match if both correspond to the same conventional name
as determined by `Var.conventional_dim_name`.

Example
==========

```jldoctest
julia> ClimaAnalysis.Var.find_corresponding_dim_name("longitude", ["lat", "lon"])
"lon"

julia> ClimaAnalysis.Var.find_corresponding_dim_name("time", ["lat", "lon", "t"])
"t"

julia> ClimaAnalysis.Var.find_corresponding_dim_name("height", ["lat", "lon", "height"])
"height"
```
"""
function find_corresponding_dim_name(dim_name::AbstractString, dim_names)
    dim_names = collect(dim_names)
    standard_dim_name = conventional_dim_name(dim_name)
    corresponding_dim_name =
        _dim_name([standard_dim_name], conventional_dim_name.(dim_names))
    isnothing(corresponding_dim_name) &&
        error("Cannot find corresponding dimension for $dim_name in $dim_names")
    return dim_names[corresponding_dim_name]
end

"""
    find_corresponding_dim_name(dim_name::AbstractString, var)

Find the corresponding dimension name in `var`'s dimension that matches `dim_name`.

Two names for a dimension match if both correspond to the same conventional name
as determined by `Var.conventional_dim_name`.

If nothing matches, return an error.

This function is useful to identify a dimension from a given name, even if it is
not its conventional name (e.g., `lat` vs `latitude`).

Example
==========

```julia-repl
julia> keys(var.dims)
("lon", "lat", "time", "potatoes")

julia> ClimaAnalysis.Var.find_corresponding_dim_name_in_var("t", var)
"time"

julia> ClimaAnalysis.Var.find_corresponding_dim_name("potatoes", var)
"potatoes"
```
"""
function find_corresponding_dim_name_in_var(dim_name::AbstractString, var)
    dim_name_in_var = try
        find_corresponding_dim_name(dim_name, keys(var.dims))
    catch
        dim_name
    end
    haskey(var.dims, dim_name_in_var) ||
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    return dim_name_in_var
end
