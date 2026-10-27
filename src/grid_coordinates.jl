# Customize these variables to allow other names
const LONGITUDE_NAMES = ["long", "lon", "longitude"]
const LATITUDE_NAMES = ["lat", "latitude"]
const TIME_NAMES = ["t", "time", "valid_time"]
const DATE_NAMES = ["date"]
const ALTITUDE_NAMES = ["z", "z_reference", "z_physical", "height"]
const PRESSURE_NAMES = ["pfull", "pressure_level"]

"""
    conventional_dim_name(dim_name::AbstractString)

Return the type of dimension as a string from longitude, latitude, time, date,
altitude, or pressure if possible or `dim_name` as a string otherwise.
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
    AbstractCoordinate

An object representing either a dimension coordinates or auxillary coordinates
with metadata attached to it such as name, units, bounds, and attributes.
"""
abstract type AbstractCoordinate end

"""
    Dim

A dimension of 1D coordinates with metadata attached to it such as name, units,
bounds, and attributes.
"""
struct Dim{A <: AbstractVector, B <: Union{Nothing, AbstractMatrix}} <:
       AbstractCoordinate
    name::String
    values::A
    units::Base.RefValue{String}
    bounds::B
    attributes::Dict{String, Any}
end

"""
    AuxCoord

A dimension of `N`-D auxillary coordinates with metadata attached to it such as
name, spans, units, bounds, and attributes. Each auxillary coordinate must span
`N` `Dim`ensions.
"""
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

"""
    Dim(
        name,
        values::AbstractVector;
        units = "",
        bounds = nothing,
        attributes = Dict{String, Any}(),
    )

Construct a `Dim`ension from its `name` and `values` with optional metadata of
`units`, `bounds`, and `attributes`.

The `values` must be one-dimensional. The `bounds`, if given, have size
`(2, length(values))`, where `bounds[:, i]` are the start and end of the `i`th
cell. An empty `units` means the dimension is dimensionless.

Also, see [`AuxCoord`](@ref).
"""
function Dim(
    name,
    values::AbstractVector;
    units = "",
    bounds = nothing,
    attributes = Dict{String, Any}(),
)
    name = String(name)
    isempty(name) && error("Dimension name cannot be empty")
    _check_bounds(name, values, bounds)
    # For a dimension, the number of vertices for bounds must always be two
    isnothing(bounds) ||
        size(bounds, 1) == 2 ||
        error(
            "Bounds of dimension $name must have 2 vertices; got $(size(bounds, 1))",
        )
    return Dim(
        name,
        values,
        Ref(String(units)),
        bounds,
        Dict{String, Any}(attributes),
    )
end

"""
    AuxCoord(
        name,
        spans,
        values::AbstractArray;
        units = "",
        bounds = nothing,
        attributes = Dict{String, Any}(),
    )

Construct a `AuxCoord` from its `name`, `spans`, and `values` with optional
metadata of `units`, `bounds`, and `attributes`.

The `spans` are the names of the dimensions that the auxiliary coordinate lies
over, given as a single name or a collection of names. There is one name for
each dimension of `values`, in the same order, with no repeats. The `bounds`, if
given, have size `(nv, size(values)...)`, where `nv` is the number of vertices
of each cell, so `bounds[:, I...]` are the vertices of the cell at
`values[I...]`. An empty `units` means the coordinate is dimensionless.
"""
function AuxCoord(
    name,
    spans,
    values::AbstractArray;
    units = "",
    bounds = nothing,
    attributes = Dict{String, Any}(),
)
    name = String(name)
    isempty(name) && error("Auxiliary coordinate name cannot be empty")
    spans = spans isa AbstractString ? (String(spans),) : Tuple(String.(spans))
    isempty(spans) &&
        error("Auxiliary coordinates with no span are not supported")

    length(spans) == ndims(values) || error(
        "Auxiliary coordinate $name spans $(length(spans)) dimensions, but its values are $(ndims(values))-D",
    )
    allunique(spans) || error(
        "Auxiliary coordinate $name repeats a dimension in its spans $spans",
    )
    _check_bounds(name, values, bounds)
    return AuxCoord(
        name,
        spans,
        values,
        Ref(String(units)),
        bounds,
        Dict{String, Any}(attributes),
    )
end

"""
    _check_bounds(name, values, bounds)

If `bounds` is not `nothing`, check that it has size `(nv, size(values)...)`,
where `nv` is the number of vertices of each cell.
"""
function _check_bounds(name, values, bounds)
    isnothing(bounds) && return nothing
    # The first dimension of bounds should represent the vertices
    (
        ndims(bounds) == ndims(values) + 1 &&
        size(bounds)[2:end] == size(values)
    ) || error(
        "Bounds of $name must have size (nv, $(join(size(values), ", "))). Got $(size(bounds))",
    )
    return nothing
end

# Coordinates
"""
    name(coord::AbstractCoordinate)

Return the `name` of `coord`.
"""
function name(coord::AbstractCoordinate)
    return coord.name
end

"""
    coord_values(coord::AbstractCoordinate)

Return the values of `coord`.

This is the stored array in `coord`, not a copy.
"""
function coord_values(coord::AbstractCoordinate)
    return coord.values
end

"""
    spans(dim::Dim)

Return the spans of `dim` which is a 1-tuple of the dimension name.
"""
function spans(dim::Dim)
    return (name(dim),)
end

"""
    spans(aux_coord::AuxCoord)

Return the spans of `aux_coord`, which are the names of the dimensions it lies
over, in the same order as the dimensions of its values.
"""
function spans(aux_coord::AuxCoord)
    return aux_coord.spans
end

"""
    units(coord::AbstractCoordinate)

Return the units of `coord`.
"""
function units(coord::AbstractCoordinate)
    return coord.units[]
end

"""
    bounds(coord::AbstractCoordinate)

Return the bounds of `coord`.

This is the stored array in `coord`, not a copy.
"""
function bounds(coord::AbstractCoordinate)
    return coord.bounds
end

"""
    attributes(coord::AbstractCoordinate)

Return the attributes of `coord`.

This is the stored dictionary in `coord`, not a copy.
"""
function attributes(coord::AbstractCoordinate)
    return coord.attributes
end

"""
    isequispaced(dim::Dim)

Return whether the values of `dim` is equispaced.
"""
function isequispaced(dim::Dim)
    return _isequispaced(coord_values(dim))
end

"""
    Base.getindex(dim::Dim, ind::Union{AbstractVector, Colon})

Return a subset of `dim` as selected by the indices `inds` as a `Dim`.

!!! note "Integer indices"
    Integer indices are not supported because that would remove the dimension.
"""
function Base.getindex(dim::Dim, ind::Union{AbstractVector, Colon})
    dim_bounds = bounds(dim)
    return Dim(
        name(dim),
        coord_values(dim)[ind];
        units = units(dim),
        bounds = isnothing(dim_bounds) ? nothing : dim_bounds[:, ind],
        attributes = deepcopy(attributes(dim)),
    )
end

# Throw an informative error message instead of a MethodError
function Base.getindex(dim::Dim, ::Integer)
    error("Indexing dimension $(name(dim)) with an integer is not allowed")
end

"""
    Base.getindex(
        aux_coord::AuxCoord,
        inds::Union{AbstractVector, Colon, Integer}...,
    )

Return a subset of `aux_coord` as selected by the indices `inds` as a
`AuxCoord`.
"""
function Base.getindex(
    aux_coord::AuxCoord,
    inds::Union{AbstractVector, Colon, Integer}...,
)
    aux_coord_spans = spans(aux_coord)
    length(inds) == length(aux_coord_spans) || error(
        "Auxiliary coordinate $(name(aux_coord)) spans $aux_coord_spans, but got $(length(inds)) indices",
    )
    kept_spans = Tuple(
        span for
        (span, ind) in zip(aux_coord_spans, inds) if !(ind isa Integer)
    )
    isempty(kept_spans) && error(
        "Indexing auxiliary coordinate $(name(aux_coord)) with only integers would leave it with no spans",
    )
    aux_bounds = bounds(aux_coord)
    return AuxCoord(
        name(aux_coord),
        kept_spans,
        coord_values(aux_coord)[inds...];
        units = units(aux_coord),
        bounds = isnothing(aux_bounds) ? nothing : aux_bounds[:, inds...],
        attributes = deepcopy(attributes(aux_coord)),
    )
end

"""
    Base.firstindex(dim::Dim)

Return the first index of `dim`.
"""
function Base.firstindex(dim::Dim)
    return firstindex(coord_values(dim))
end

"""
    Base.lastindex(dim::Dim)

Return the last index of `dim`.
"""
function Base.lastindex(dim::Dim)
    return lastindex(coord_values(dim))
end

"""
    Base.firstindex(aux_coord::AuxCoord, d)

Return the first index along the `d`-th dimension of `aux_coord`.
"""
function Base.firstindex(aux_coord::AuxCoord, d)
    firstindex(coord_values(aux_coord), d)
end

"""
    Base.lastindex(aux_coord::AuxCoord, d)

Return the last index along the `d`-th dimension of `aux_coord`.
"""
function Base.lastindex(aux_coord::AuxCoord, d)
    lastindex(coord_values(aux_coord), d)
end

# TODO: show(io::IO, ::MIME"text/plain", x)
# TODO: Fix this for string intp and vectors
# TODO: 0 elements vs 0 element

"""
    Base.show(io::IO, dim::Dim)

TODO: Write docs and review this.
"""
function Base.show(io::IO, dim::Dim)
    printstyled(io, "Attributes:\n", bold = true, color = :green)
    _print_attributes(io, _shown_attributes(dim), "  ")
    printstyled(io, "Dimension:\n", bold = true, color = :green)
    _print_coord(io, dim, length(name(dim)))
end

"""
    Base.show(io::IO, aux_coord::AuxCoord)

TODO: Write docs and review this.
"""
function Base.show(io::IO, aux_coord::AuxCoord)
    printstyled(io, "Attributes:\n", bold = true, color = :green)
    _print_attributes(io, _shown_attributes(aux_coord), "  ")
    printstyled(io, "Auxiliary coordinate:\n", bold = true, color = :green)
    _print_coord(io, aux_coord, length(name(aux_coord)))
end

# Units are stored apart from the other attributes, so show them first
function _shown_attributes(coord::AbstractCoordinate)
    shown = Pair{String, Any}[]
    isempty(units(coord)) || push!(shown, "units" => units(coord))
    return append!(shown, attributes(coord))
end

function _print_attributes(io, attribs, indent)
    max_length = maximum(length(key) for (key, _) in attribs; init = 0)
    for (key, val) in attribs
        print(io, indent * rpad(key, max_length) * " => ")
        printstyled(io, "$(val)\n", color = :light_cyan)
    end
    return nothing
end

# One line per coordinate, with the name padded to name_width
function _print_coord(io, dim::Dim, name_width)
    printstyled(io, "  " * rpad(name(dim), name_width), color = :light_green)
    print(io, " with ")
    values = coord_values(dim)
    printstyled(io, "$(length(values)) ", color = :light_cyan)
    if length(values) >= 2
        print(io, "elements (")
        if issorted(values) || issorted(values, rev = true)
            printstyled(io, "$(first(values))", color = :light_cyan)
            print(io, " to ")
            printstyled(io, "$(last(values))", color = :light_cyan)
        else
            printstyled(io, "not sorted", color = :red)
        end
        print(io, ")")
    elseif length(values) == 1
        print(io, "element (")
        printstyled(io, "$(first(values))", color = :light_cyan)
        print(io, ")")
    else
        print(io, "elements")
    end
    isnothing(bounds(dim)) || print(io, " and bounds")
    return nothing
end

# Aux coordinates need not be sorted, so show the range of their values
function _print_coord(io, aux_coord::AuxCoord, name_width)
    printstyled(
        io,
        "  " * rpad(name(aux_coord), name_width),
        color = :light_green,
    )
    print(io, " over ($(join(spans(aux_coord), ", "))) with ")
    values = coord_values(aux_coord)
    printstyled(io, join(size(values), "×"), color = :light_cyan)
    print(io, length(values) == 1 ? " element" : " elements")
    numbers = skipmissing(values)
    if nonmissingtype(eltype(values)) <: Real && !isempty(numbers)
        low, high = extrema(numbers)
        print(io, " (")
        printstyled(io, "$low", color = :light_cyan)
        print(io, " to ")
        printstyled(io, "$high", color = :light_cyan)
        print(io, ")")
    end
    isnothing(bounds(aux_coord)) || print(io, " and bounds")
    return nothing
end
