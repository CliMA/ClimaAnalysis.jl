module Grids

export AbstractCoordinate, Dim, AuxCoord, AbstractGrid, Grid
export AbstractSelector, NearestValue, MatchValue, Index
export dim, aux_coord, coordinate, dim_names, aux_coord_names, dim_index
export hasdim, is_z_1D, reference_date, calendar
export has_time, has_date, has_longitude, has_latitude, has_altitude
export has_pressure
export time_name, date_name, longitude_name, latitude_name, altitude_name
export pressure_name
export times, longitudes, latitudes, altitudes, pressures
export spans, units, attributes, isequispaced
export set_units!, set_values!, transform_coords!
export select_indices, get_index, window_indices, selection_record
export reduce_dims, provenance, replace_dim, transform_coords, rename_dim
export permutation, shift_longitude, shift_longitude_indices, swap_dim
export add_aux_coords, drop_aux_coords, replace_aux_coord, rename_aux_coord
export replace_bounds, guess_bounds
export coordinate_field, latitude_weights, cell_widths
export iscompatible, check_compatible, combine_grids
export read_grid, default_short_name, aux_discovery

abstract type AbstractCoordinate end

# TODO: Maybe rename to DimCoord
struct Dim{A <: AbstractVector, B <: Union{Nothing, AbstractMatrix}} <:
       AbstractCoordinate
    name::String
    values::A
    units::Base.RefValue{String}
    bounds::B
    attributes::Dict{String, Any}
end

struct AuxCoord{
    N,
    A <: AbstractArray{<:Any, N},
    B <: Union{Nothing, AbstractArray},
} <: AbstractCoordinate
    name::String
    spans::NTuple{N, String}
    values::A
    units::Base.RefValue{String}
    bounds::B
    attributes::Dict{String, Any}
end

abstract type AbstractGrid end

struct Grid{D <: Tuple{Vararg{Dim}}, X <: Tuple{Vararg{AuxCoord}}} <:
       AbstractGrid
    dims::D
    aux_coords::X
end

"""
    AbstractSelector

An object that determines which indices are selected.

`AbstractSelector`s have to provide one function, `get_index`

The function has to have the signature
`get_index(var, dim_name, idx_or_val, ::AbstractSelector)` and return a single index. You
can assume that `dim_name` is in `ClimaAnalysis.dim_names(var)`.

The function `get_index` is used by [`slice`](@ref) and [`window`](@ref). For instance, if
you use `ClimaAnalysis.slice(var, time = 2)`, then `dim_name` is `time`, and `idx_or_val` is
`2`.

For example, to implement an `AbstractSelector` that always take the first index of the
dimension regardless of the value or index passed in, you can write the following:

```julia
struct FirstIndex <: ClimaAnalysis.AbstractSelector end

function ClimaAnalysis.get_index(var, dim_name, idx_or_val, ::FirstIndex)
    return 1
end

# Get the first time slice of var. The parameter 10 does not do anything.
ClimaAnalysis.slice(var, time = 10, by = FirstIndex())
```
"""
abstract type AbstractSelector end

"""
    NearestValue

Get the index of the nearest value.

If the dimension is not one dimesional, then an error is thrown.
"""
struct NearestValue <: AbstractSelector end

"""
    MatchValue

Get the index of the approximately matched value.

If the value does not exist, or the dimension is not one-dimensional, then an error is
thrown.
"""
struct MatchValue <: AbstractSelector end

"""
    Index

Select the index in the dimension.

If the index is out of bounds, then an error is thrown.
"""
struct Index <: AbstractSelector end

# Construction
# TODO: Change the name of this since we are just checking the data size match
# the grid size
function check_grid_consistency(grid::Grid, data::AbstractArray) end

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
    conventional_name = conventional_dim_name(name)
    for coord in coords
        if conventional_name == conventional_dim_name(coord.name)
            return coord
        end
    end
    return nothing
end

function dim_names(grid::Grid)
    return Tuple(dim.name for dim in grid.dims)
end

function aux_coord_names(grid::Grid)
    return Tuple(aux_coord.name for aux_coord in grid.aux_coords)
end

function dim_index(grid::Grid, name)
    (; dims) = grid
    name = conventional_dim_name(name)
    for (i, dim) in enumerate(dims)
        if conventional_dim_name(dim.name) == name
            return i
        end
    end
    return error(
        "Cannot find index of the dimension with the name $name in the grid",
    )
end

function hasdim(grid::Grid, name)
    (; dims) = grid
    name = conventional_dim_name(name)
    for dim in dims
        if conventional_dim_name(dim.name) == name
            return true
        end
    end
    return false
end

function Base.getindex(grid::Grid, name::AbstractString) end
function Base.haskey(grid::Grid, name) end
function Base.size(grid::Grid) end
function Base.ndims(grid::Grid) end
function Base.isempty(grid::Grid) end

function is_z_1D(grid)
    # Find dimension with z coordinates
    # Check for aux coords I think?
    # TODO: Look at the existing netcdf file for this
    error("TODO!")
end

function reference_date(grid::Grid)
    # TODO: Use CFTime to parse the units attribute and get the start date from that
    return nothing
end

function has_time(grid::Grid)
    return _has_dimname(grid, "time")
end

function has_date(grid::Grid)
    # TODO: Check units attribute to see if dates make sense
    return has_time(grid)
end

function has_longitude(grid::Grid)
    return _has_dimname(grid, "longitude")
end

function has_latitude(grid::Grid)
    return _has_dimname(grid, "llatitude")
end

function has_altitude(grid::Grid)
    return _has_dimname(grid, "altitude")
end

function has_pressure(grid::Grid)
    return _has_dimname(grid, "pressure")
end

function _has_dimname(grid::Grid, dimname)
    # TODO: This is wrong since CF conventions said that we identify
    # dimensions by their units and not their names
    (; dims) = grid
    for dim in dims
        if conventional_dim_name(dim.name) == dimname
            return true
        end
    end
    return false
end

function time_name(grid::Grid)
    return _find_dimname(grid, "time")
end

function longitude_name(grid::Grid)
    return _find_dimname(grid, "longitude")
end

function latitude_name(grid::Grid)
    return _find_dimname(grid, "longitude")
end

function altitude_name(grid::Grid)
    return _find_dimname(grid, "altitude")
end

function pressure_name(grid::Grid)
    return _find_dimname(grid, "pressure")
end

# TODO: Not sure to support coord too or just dim only
function _find_dimname(grid::Grid, dimname)
    (; dims) = grid
    for dim in dims
        if conventional_dim_name(dim.name) == dimname
            return dim.name
        end
    end
    error("Grid does not have $dimname among its dimensions")

end

function times(grid::Grid) end

function dates(grid::Grid) end

function longitudes(grid::Grid) end

function latitudes(grid::Grid) end

function altitudes(grid::Grid) end

function pressures(grid::Grid) end

# Coordinates
# TODO: Think about consistency between including get or not
function get_name(coord::AbstractCoordinate)
    return coord.name
end

function get_coord_values(coord::AbstractCoordinate)
    return coord.values
end

function spans(dim::Dim)
    return (get_name(dim),)
end

function spans(aux_coord::AuxCoord)
    return aux_coord.spans
end

function units(coord::AbstractCoordinate)
    return coord.units
end

function bounds(coord::AbstractCoordinate)
    return coord.bounds
end

function attributes(coord::AbstractCoordinate)
    return coord.attributes
end

# TODO: Could be aux coords if 1D?
function isequispaced(dim::Dim)
    dim_vals = get_coord_values(dim)
    # TODO: Use is_equispaced from Utils
    return all(diff(dim_vals) .≈ dim_vals[begin + 1] - dim_vals[begin])
end

function set_units!(grid::Grid, name, u::AbstractString)
    coord = coordinate(grid, name)
    coord.units[] = u
    return nothing
end

function set_values!(grid::Grid, name, new_values; bounds = nothing)
    coord = coordinate(grid, name)
    # TODO: Do we want to support broadcasting here?
    coord.values .= new_values
    isnothing(bounds) || (coord.bounds .= bounds)
    return nothing
end

function transform_coords!(grid::Grid, (name, f)::Pair; units = nothing) end

# Selection and structural indexing

function get_index(grid::Grid, dim_name, value, ::NearestValue) end
function get_index(grid::Grid, dim_name, value, ::MatchValue) end
function get_index(grid::Grid, dim_name, value, ::Index) end

function select(grid::Grid; by::AbstractSelector = NearestValue(), inds...) end

function select_indices(
    grid::Grid;
    by::AbstractSelector = NearestValue(),
    inds...,
) end

function window_indices(
    grid::Grid,
    name;
    left = nothing,
    right = nothing,
    by::AbstractSelector = NearestValue(),
) end


function _selection_record(grid::Grid, inds::Tuple) end

# Operation verbs
function reduce_dims(grid::Grid; dims) end
function provenance(old::Grid, new::Grid) end
function replace_dim(grid::Grid, (name, new_dim)::Pair{<:AbstractString, <:Dim}) end
function replace_dim(
    grid::Grid,
    (name, new_values)::Pair{<:AbstractString, <:AbstractVector},
) end
function transform_coords(grid::Grid, (name, f)::Pair; units = nothing) end
function rename_dim(grid::Grid, (old_name, new_name)::Pair) end
function permutation(grid::Grid, perm) end
function shift_longitude(grid::Grid, lower_lon::Real, upper_lon::Real) end
function shift_longitude_indices(grid::Grid, lower_lon::Real, upper_lon::Real) end
function swap_dim(grid::Grid, (dim_name, aux_coord_name)::Pair) end
function Base.dropdims(grid::Grid; dims) end
function Base.reverse(grid::Grid; dims) end
function Base.permutedims(grid::Grid, perm) end
function Base.cat(grid::Grid, grids::Grid...; dims) end
function Base.ones(grid::Grid) end
function Base.ones(::Type{T}, grid::Grid) where {T} end

# Auxiliary coordinate and bounds management
function add_aux_coords(grid::Grid, aux_coords::AuxCoord...) end
function drop_aux_coords(grid::Grid, names::AbstractString...) end
function replace_aux_coord(
    grid::Grid,
    (name, new_aux_coord)::Pair{<:AbstractString, <:AuxCoord},
) end
function rename_aux_coord(grid::Grid, (old_name, new_name)::Pair) end
function replace_bounds(
    grid::Grid,
    (name, new_bounds)::Pair{<:AbstractString, <:Union{Nothing, AbstractArray}},
) end
function guess_bounds(grid::Grid, name) end

# Geometry and interpolation support
function coordinate_field(grid::Grid, name) end
function latitude_weights(grid::Grid) end
function cell_widths(dim::Dim) end
function make_interpolant(grid::Grid, data::AbstractArray) end
function extrapolation_bc(dim::Dim) end

# Compatibility
function iscompatible(a::Grid, b::Grid) end
function check_compatible(a::Grid, b::Grid; dims = dim_names(a)) end
function combine_grids(a::Grid, b::Grid) end
function Base.isapprox(a::Grid, b::Grid; kwargs...) end

# Reading files
function read_grid(ds, short_name::AbstractString) end
function default_short_name(ds) end
function aux_discovery(ds, short_name::AbstractString) end

# Equality, hashing, display
function Base.:(==)(dim1::Dim, dim2::Dim) end
function Base.:(==)(aux_coord1::AuxCoord, aux_coord2::AuxCoord) end
function Base.:(==)(a::Grid, b::Grid) end
function Base.hash(dim::Dim, h::UInt) end
function Base.hash(aux_coord::AuxCoord, h::UInt) end
function Base.hash(grid::Grid, h::UInt) end
function Base.copy(dim::Dim) end
function Base.copy(aux_coord::AuxCoord) end
function Base.copy(grid::Grid) end
function Base.show(io::IO, dim::Dim) end
function Base.show(io::IO, aux_coord::AuxCoord) end
function Base.show(io::IO, grid::Grid) end
function Base.summary(io::IO, grid::Grid) end

end
