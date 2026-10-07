# Customize these variables to allow other names
const LONGITUDE_NAMES = ["long", "lon", "longitude"]
const LATITUDE_NAMES = ["lat", "latitude"]
const TIME_NAMES = ["t", "time", "valid_time"]
const DATE_NAMES = ["date"]
const ALTITUDE_NAMES = ["z", "z_reference", "z_physical", "height"]
const PRESSURE_NAMES = ["pfull", "pressure_level"]

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

# TODO: CHECK BELOW
# Constructors
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

Dim(name, values::AbstractArray; kwargs...) = error(
    "Dimension $name must be 1-D, but its values are $(ndims(values))-D. Use AuxCoord instead",
)

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
    spans =
        spans isa Union{AbstractString, Symbol} ? (String(spans),) :
        Tuple(String.(spans))
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

function _check_bounds(name, values, bounds)
    # TODO: Add examples of what is mean for bounds to be > 2 dimensional
    isnothing(bounds) && return nothing
    (
        ndims(bounds) == ndims(values) + 1 &&
        size(bounds)[2:end] == size(values)
    ) || error(
        "Bounds of $name must have size (nv, $(join(size(values), ", "))); got $(size(bounds))",
    )
    return nothing
end

# Coordinates
function name(coord::AbstractCoordinate)
    return coord.name
end

function coord_values(coord::AbstractCoordinate)
    return coord.values
end

function spans(dim::Dim)
    return (name(dim),)
end

function spans(aux_coord::AuxCoord)
    return aux_coord.spans
end

function units(coord::AbstractCoordinate)
    return coord.units[]
end

function bounds(coord::AbstractCoordinate)
    return coord.bounds
end

function attributes(coord::AbstractCoordinate)
    return coord.attributes
end

# TODO: Could be aux coords if 1D?
function isequispaced(dim::Dim)
    return _isequispaced(coord_values(dim))
end

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

function Base.getindex(dim::Dim, ind::Integer)
    error(
        "Indexing dimension $(name(dim)) with an integer is not allowed",
    )
end

function Base.getindex(aux_coord::AuxCoord, inds...)
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
