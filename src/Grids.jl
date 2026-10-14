module Grids

import NCDatasets: CFTime
import NCDatasets
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

function read_grid(path::String, short_name)
    return NCDatasets.NCDataset(path) do nc
        read_grid(nc, short_name)
    end
end

"""
    read_grid(nc::NCDatasets.NCDataset, short_name)

Create a `Grid` from the variable with the name `short_name` in the NetCDF
dataset `nc`.

This function does not close `nc`.
"""
function read_grid(nc::NCDatasets.NCDataset, short_name)
    # We will not try to automatically determine the short name here since this
    # function will be called by the OutputVar constructor
    haskey(nc, short_name) ||
        error("Variable $short_name is not available in the NetCDF file")

    dim_names = NCDatasets.dimnames(nc[short_name])
    dims = map(dim_names) do dim_name
        # TODO: Check the line below since I haven't seen it and maybe it is
        # mentioned in the CF conventions?
        # A dimension without a coordinate variable becomes an index axis
        haskey(nc, dim_name) ||
            return Dim(dim_name, collect(1:nc.dim[dim_name]))
        variable = nc[dim_name]
        dim_attribs = Dict{String, Any}(variable.attrib)
        units = pop!(dim_attribs, "units", "")
        bounds_name = pop!(dim_attribs, "bounds", nothing)
        # The number of vertices for the bounds of Dim is 2.
        bounds = _read_bounds(
            nc,
            dim_name,
            bounds_name,
            (dim_name,);
            nvertices = 2,
        )

        # TODO: Need to handle the time dimension as a special case since it
        # could be some Dates.DateTime

        return Dim(
            dim_name,
            _coord_values(variable);
            units,
            bounds,
            attributes = dim_attribs,
        )
    end

    # TODO: Is this split right?
    aux_coord_names = split(get(nc[short_name].attrib, "coordinates", " "))
    unique!(aux_coord_names)
    aux_coord_names = filter(
        name -> _is_valid_aux_coord(nc, name, dim_names),
        aux_coord_names,
    )
    aux_coords = map(aux_coord_names) do aux_coord_name
        variable = nc[aux_coord_name]
        spans = _coord_dimnames(variable)
        aux_coord_attribs = Dict{String, Any}(variable.attrib)
        units = pop!(aux_coord_attribs, "units", "")
        bounds_name = pop!(aux_coord_attribs, "bounds", nothing)
        # Aux coordinate cells can be polygons, so the vertex count is not fixed
        bounds = _read_bounds(nc, aux_coord_name, bounds_name, spans)
        return AuxCoord(
            aux_coord_name,
            spans,
            _coord_values(variable);
            units,
            bounds,
            attributes = aux_coord_attribs,
        )
    end

    return Grid(dims...; aux_coords)
end

"""
TODO: This entire docstring

Read the bounds named by a coordinate's `bounds` attribute (CF 7.1), or warn and return
nothing if they are missing or malformed.
"""
function _read_bounds(
    nc,
    coord_name,
    bounds_name,
    coord_dims;
    nvertices = nothing,
)
    isnothing(bounds_name) && return nothing
    if !haskey(nc, bounds_name)
        @warn "Bounds $bounds_name of $coord_name are not in the file; skipping them"
        return nothing
    end

    bounds_dims = _coord_dimnames(nc[bounds_name])
    # The number of dimensions for the bounds is one more than the number of
    # dimensions of the coordinates due to the bound dimension. Furthermore,
    # bound dimension should be the first one and not present in the coordinate
    # dimensions.
    is_valid_layout =
        bounds_dims[2:end] == Tuple(coord_dims) && bounds_dims[1] ∉ coord_dims
    if !is_valid_layout
        @warn "Bounds $bounds_name of $coord_name have dimensions $bounds_dims, but a vertex dimension followed by $(Tuple(coord_dims)) is expected; skipping them"
        return nothing
    end

    # Check the size of the bounds match
    file_nvertices = nc.dim[bounds_dims[1]]
    if !isnothing(nvertices) && file_nvertices != nvertices
        @warn "Bounds $bounds_name of $coord_name have $file_nvertices vertices, but $nvertices are expected; skipping them"
        return nothing
    end
    return _coord_values(nc[bounds_name])
end

# TODO: Verify if _coord_dimnames and _coord_values work correctly for strings
# A char array's first dimension holds the characters of each string (CF §2.2)
function _coord_dimnames(variable)
    dims = NCDatasets.dimnames(variable)
    return eltype(variable) <: Char ? dims[2:end] : dims
end

function _coord_values(variable)
    values = Array(variable)
    eltype(values) <: Char || return values
    # CF pads strings with trailing NULs or spaces
    return map(eachslice(values; dims = Tuple(2:ndims(values)))) do chars
        String(rstrip(c -> c == '\0' || c == ' ', String(UInt8.(chars))))
    end
end

# TODO for read_grid, easiest first
# Tidy-ups (no behaviour change)
# 5. Factor the duplicated units/bounds/attributes extraction into one helper
# Small behaviour changes
# 6. The "not available" error lists the variables that do exist
# 7. Accept any NCDatasets.AbstractNCDataset, e.g. multi-file datasets (NCCatalog)
# 8. A variable named like its dim but not 1-D over it (legacy z(z, lat, lon)):
#    warn and use a 1:n axis instead of letting Dim throw (D-7)
# 9. Reject ragged DSG files (sample_dimension/instance_dimension) (IO-11)
# Need a decision first
# 10. Encoding attributes (_FillValue, scale_factor, ...) end up in attributes:
#     drop them or keep them apart (CR-5)
# 11. _FillValue makes values Union{Missing, T}, e.g. CALIPSO lat/lon: convert
#     to NaN or keep missing
# 12. Warn when bounds units disagree with the coordinate's (UN-4, units PR)
# Hard
# 13. Decode time to relative floats, not DateTime (DM-14), including time
#     bounds and date(time): NCDatasets decodes bounds with the parent's
#     units/calendar
# 14. Reference date from time units, then date units, then start_date (IO-13)
# 15. Aux discovery tiers 3-4, e.g. unlinked z_physical, without promoting data
#     variables (IO-3, IO-4)


function _is_valid_aux_coord(nc, aux_coord_name, dim_names)
    aux_coord_name in dim_names && return false
    if !haskey(nc, aux_coord_name)
        @warn "$aux_coord_name is not in the file; skipping it"
        return false
    end
    # We do not support 0D spans right now
    spans = _coord_dimnames(nc[aux_coord_name])
    if isempty(spans)
        @warn "$aux_coord_name is a scalar coordinate; skipping it"
        return false
    end
    # TODO: Mention edge cases here
    if !issubset(spans, dim_names)
        @warn "$aux_coord_name has dimensions $spans not in $dim_names; skipping it"
        return false
    end
    return true
end

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
    # TODO: Better error message when it is an auxiliary coordinate
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
