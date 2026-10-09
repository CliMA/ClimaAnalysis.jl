module Grids

import NCDatasets: CFTime
import Dates
import ..Utils: nearest_index, _isequispaced


# TODO: Determine if we should export these or not
export AbstractCoordinate, Dim, AuxCoord, AbstractGrid, Grid
export AbstractSelector, NearestValue, MatchValue, Index
export dim, aux_coord, coordinate, dim_names, aux_coord_names, dim_index
export has_dim, has_aux_coord, has_coord
export is_z_1D, reference_date, calendar
export has_time, has_date, has_longitude, has_latitude, has_altitude
export has_pressure
export time_name, longitude_name, latitude_name, altitude_name, pressure_name
export times, longitudes, latitudes, altitudes, pressures
export spans, units, attributes, isequispaced
export set_units!, set_values!, transform_coords!
export select_indices, get_index, window_indices

include("grid_coordinates.jl")

abstract type AbstractGrid end

struct Grid{D <: Tuple{Vararg{Dim}}, X <: Tuple{Vararg{AuxCoord}}} <:
       AbstractGrid
    dims::D
    aux_coords::X
end

function Grid(dims::Dim...; aux_coords = ())
    aux_coords = Tuple(aux_coords)
    dim_names_of_grid = map(name, dims)
    aux_coord_names_of_grid = map(name, aux_coords)

    # Check coordinate names are all different
    coord_names = (dim_names_of_grid..., aux_coord_names_of_grid...)
    duplicate_names =
        unique(n for n in coord_names if count(==(n), coord_names) > 1)
    isempty(duplicate_names) || error(
        "Coordinate names must be unique, but $(join(duplicate_names, ", ")) appear more than once",
    )

    # Check spans of aux_coords is in the names of dims
    for aux_coord in aux_coords
        unknown_spans = filter(∉(dim_names_of_grid), spans(aux_coord))
        isempty(unknown_spans) || error(
            "Auxiliary coordinate $(name(aux_coord)) spans $(join(unknown_spans, ", ")), which are not dimensions of the grid $dim_names_of_grid",
        )
    end

    # Check size of aux_coords match dimensions
    dim_lengths = Dict(name(dim) => length(coord_values(dim)) for dim in dims)
    for aux_coord in aux_coords
        expected_size = Tuple(dim_lengths[span] for span in spans(aux_coord))
        actual_size = size(coord_values(aux_coord))
        actual_size == expected_size || error(
            "Auxiliary coordinate $(name(aux_coord)) has size $actual_size, but its spans $(spans(aux_coord)) have lengths $expected_size",
        )
    end

    # TODO: Does it make sense for aux_coords to be a tuple
    return Grid(dims, aux_coords)
end

# TODO: This need CFTime support
function calendar(grid::Grid)
    has_time(grid) || return nothing
    return get(attributes(dim(grid, "time")), "calendar", "standard")
end

# TODO: Implement this
function read_grid() end

# Construction
# TODO: Change the name of this since we are just checking the data size match
# the grid size
function check_grid_consistency(grid::Grid, data::AbstractArray)
    return size(grid) == size(data)
end

# Lookup
function dim(grid::Grid, name)
    maybe_dim = _find_coord(grid.dims, name)
    return isnothing(maybe_dim) ?
           error("Cannot find dimension with name $name") : maybe_dim
end

function aux_coord(grid::Grid, name)
    maybe_aux_coord = _find_coord(grid.aux_coords, name)
    return isnothing(maybe_aux_coord) ?
           error("Cannot find auxiliary coordinate with name $name") :
           maybe_aux_coord
end

function coordinate(grid::Grid, name)
    maybe_coord = _find_coord((grid.dims..., grid.aux_coords...), name)
    return isnothing(maybe_coord) ?
           error("Cannot find coordinate with name $name") : maybe_coord
end

function _find_coord(coords, name)
    # Exact match first: z_reference and z_physical share a conventional name
    for coord in coords
        coord.name == name && return coord
    end
    conventional_name = conventional_dim_name(name)
    for coord in coords
        if conventional_name == conventional_dim_name(coord.name)
            return coord
        end
    end
    return nothing
end

function dim_names(grid::Grid)
    return map(name, grid.dims)
end

function aux_coord_names(grid::Grid)
    return map(name, grid.aux_coords)
end

function dim_index(grid::Grid, name)
    (; dims) = grid
    name = conventional_dim_name(name)
    for (i, dim) in enumerate(dims)
        if conventional_dim_name(dim.name) == name
            return i
        end
    end
    # TODO: Better error message when it is an auxillary coordinate
    return error(
        "Cannot find index of the dimension with the name $name in the grid",
    )
end

function has_dim(grid::Grid, name)
    return !isnothing(_find_coord(grid.dims, name))
end

function has_aux_coord(grid::Grid, name)
    return !isnothing(_find_coord(grid.aux_coords, name))
end

function has_coord(grid::Grid, name)
    return !isnothing(_find_coord((grid.dims..., grid.aux_coords...), name))
end

function Base.size(grid::Grid)
    return map(dim -> length(coord_values(dim)), grid.dims)
end

function Base.ndims(grid::Grid)
    return length(grid.dims)
end

function Base.isempty(grid::Grid)
    return isempty(grid.dims) && isempty(grid.aux_coords)
end

function is_z_1D(grid)
    # Find dimension with z coordinates
    # Check for aux coords I think?
    # TODO: Look at the existing netcdf file for this
    error("TODO!")
end

# TODO: Check the dates functions like dates, has_date,
# and date_to_time for conciseness and correctness and type stability

function reference_date(grid::Grid)
    # TODO: Check this function too. It looks weird to me
    has_time(grid) || return nothing
    time_units = units(dim(grid, "time"))
    occursin(" since ", time_units) || return nothing
    cal = calendar(grid)
    ref_date = CFTime.origin(CFTime.timetype(cal, time_units))
    return cal in ("standard", "gregorian", "proleptic_gregorian") ?
           Dates.DateTime(ref_date) : ref_date
end

# TODO: This is wrong since CF conventions said that we identify
# dimensions by their units and not their names
function has_time(grid::Grid)
    return has_dim(grid, "time")
end

function has_date(grid::Grid)
    return has_time(grid) && !isnothing(reference_date(grid))
end

function has_longitude(grid::Grid)
    return has_dim(grid, "longitude")
end

function has_latitude(grid::Grid)
    return has_dim(grid, "latitude")
end

function has_altitude(grid::Grid)
    return has_dim(grid, "altitude")
end

function has_pressure(grid::Grid)
    return has_dim(grid, "pressure")
end

function time_name(grid::Grid)
    return name(dim(grid, "time"))
end

function longitude_name(grid::Grid)
    return name(dim(grid, "longitude"))
end

function latitude_name(grid::Grid)
    return name(dim(grid, "latitude"))
end

function altitude_name(grid::Grid)
    return name(dim(grid, "altitude"))
end

function pressure_name(grid::Grid)
    return name(dim(grid, "pressure"))
end

# Only get dimensions which are 1D

function times(grid::Grid)
    return coord_values(dim(grid, "time"))
end

function dates(grid::Grid)
    time_dim = dim(grid, "time")
    has_date(grid) || error(
        "Cannot compute dates because the units of $(name(time_dim)) ($(units(time_dim))) have no reference date",
    )
    cal = calendar(grid)
    time_dates = CFTime.timedecode(coord_values(time_dim), units(time_dim), cal)
    eltype(time_dates) <: Dates.DateTime && return time_dates
    return map(date -> _calendar_date(cal, date), time_dates)
end

function longitudes(grid::Grid)
    return coord_values(dim(grid, "longitude"))
end

function latitudes(grid::Grid)
    return coord_values(dim(grid, "latitude"))
end

function altitudes(grid::Grid)
    return coord_values(dim(grid, "altitude"))
end

function pressures(grid::Grid)
    return coord_values(dim(grid, "pressure"))
end

function set_units!(grid::Grid, name, u::AbstractString)
    coord = coordinate(grid, name)
    coord.units[] = u
    return nothing
end

function set_values!(grid::Grid, name, new_values; bounds = nothing)
    coord = coordinate(grid, name)
    # TODO: Do we want to support broadcasting here?
    size(new_values) == size(coord.values) || error(
        "New values of $(coord.name) have size $(size(new_values)), but $(coord.name) has size $(size(coord.values))",
    )
    if !isnothing(bounds)
        isnothing(coord.bounds) && error(
            "Cannot set bounds of $(coord.name) in place because it has no bounds",
        )
        size(bounds) == size(coord.bounds) || error(
            "New bounds of $(coord.name) have size $(size(bounds)), but its bounds have size $(size(coord.bounds))",
        )
    end
    coord.values .= new_values
    isnothing(bounds) || (coord.bounds .= bounds)
    return nothing
end

function transform_coords!(f, grid::Grid, name; units = nothing)
    coord = coordinate(grid, name)
    coord_bounds = bounds(coord)
    # Compute everything before mutating, so a failing f leaves the grid unchanged
    new_values = f.(coord_values(coord))
    new_bounds = isnothing(coord_bounds) ? nothing : f.(coord_bounds)
    set_values!(grid, name, new_values; bounds = new_bounds)
    isnothing(units) || set_units!(grid, name, units)
    return nothing
end

function _date_to_time(grid::Grid, date)
    has_date(grid) ||
        error("$date is a date, but the time dimension has no reference date")
    cal = calendar(grid)
    return CFTime.timeencode(
        _to_calendar(date, cal),
        units(dim(grid, "time")),
        cal,
    )
end

# A DateTime is read as the same calendar date in the time dimension's calendar
function _to_calendar(date::Dates.DateTime, cal)
    cal in ("standard", "gregorian", "proleptic_gregorian") && return date
    # TODO: Look into whether we need to do this?
    return _calendar_date(cal, date)
end

function _to_calendar(date, cal)
    return date
end

# CFTime date types encode their units and reference date, so rebuild dates with the
# default type parameters to get one type per calendar
function _calendar_date(cal, date)
    return CFTime.timetype(cal)(
        Dates.year(date),
        Dates.month(date),
        Dates.day(date),
        Dates.hour(date),
        Dates.minute(date),
        Dates.second(date),
        Dates.millisecond(date),
    )
end

function Base.getindex(grid::Grid, inds...)
    length(inds) == ndims(grid) || error(
        "Grid has $(ndims(grid)) dimensions $(dim_names(grid)), but got $(length(inds)) indices",
    )
    ind_of = Dict(zip(dim_names(grid), inds))
    new_dims = Tuple(
        dim[ind_of[name(dim)]] for
        dim in grid.dims if !(ind_of[name(dim)] isa Integer)
    )
    new_aux_coords = Tuple(
        aux_coord[(ind_of[span] for span in spans(aux_coord))...] for
        aux_coord in grid.aux_coords if
        !all(ind_of[span] isa Integer for span in spans(aux_coord))
    )
    return Grid(new_dims...; aux_coords = new_aux_coords)
end

include("grid_selectors.jl")

end
