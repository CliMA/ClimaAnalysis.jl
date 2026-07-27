module Interpolation

import OrderedCollections: OrderedDict

"""
    Node

A staggering where the coordinates are the edges of the cells, so the first
and last coordinates are the boundaries of the domain.
"""
struct Node end

"""
    Center

A staggering where the coordinates are the centers of the cells, so the domain
extends half a cell beyond the first and last coordinates.
"""
struct Center end

"""
    AbstractExtrapolation

An object that determines how values outside of the domain are handled when
interpolating.
"""
abstract type AbstractExtrapolation end

"""
    Throw

A boundary condition where an error is thrown if it is outside of the domain.
"""
struct Throw <: AbstractExtrapolation end

"""
    Flat

A boundary condition where values outside of the domain evaluate to the value
at the nearest boundary point.
"""
struct Flat <: AbstractExtrapolation end

"""
    Periodic

A boundary condition where the domain repeats with a physical period of
`period`, so values outside of the domain are wrapped back into it.
"""
struct Periodic{FT} <: AbstractExtrapolation
    period::FT
end

"""
    AbstractInterpolationMethod

An object that determines the interpolation method of the regridder.
"""
abstract type AbstractInterpolationMethod end

"""
    LinearInterpolation

Linear interpolation, where NaNs propagate through the weighted sum.
"""
struct LinearInterpolation <: AbstractInterpolationMethod end

"""
    Regridder

Linearly interpolate gridded data with support for node and centers and throw,
flat, and periodic boundary conditions.
"""
struct Regridder{
    N,
    FT,
    SG <: NTuple{N, AbstractVector},
    A <: AbstractArray,
    E <: NTuple{N, AbstractExtrapolation},
    S <: NTuple{N, Union{Node, Center}},
    M <: AbstractInterpolationMethod,
}

    "Coordinates of the source data"
    src_grid::SG

    "Source data to interpolate from`"
    src_data::A

    "Extrapolation conditions for how queries outside the source coordinates are
    handled"
    extrapolation::E

    "Whether the coordinates of each dimension are cell centers (`Center`) or
    cell edges (`Node`)."
    staggering::S
    "Interpolation method (`LinearInterpolation`)"
    method::M

    "Preallocated per-dimension buffers of 1D stencils
    `(index1, index2, weight)`, one per destination coordinate."
    stencil_cache::NTuple{N, Vector{Tuple{Int, Int, FT}}}
end

"""
    Regridder(
        src_data::AbstractArray,
        src_grid,
        extrapolation,
        staggering,
        method = LinearInterpolation(),
    )

Construct a `Regridder` object from `src_data`, `src_grid`, `extrapolation`,
and `staggering`, interpolating with `method`.

The argument `src_grid` is a tuple of strictly increasing coordinate vectors,
one per dimension of `src_data`. The argument `extrapolation` is an iterable of
`AbstractExtrapolation` (e.g. `(Periodic(360.0), Flat())`) and `staggering` an
iterable of `Node` or `Center`, one entry per dimension.

The number of dimensions of `src_data` and `src_grid` and the length of
`extrapolation` and `staggering` must all be the same.
"""
function Regridder(
    src_data::AbstractArray,
    src_grid,
    extrapolation,
    staggering,
    method::AbstractInterpolationMethod = LinearInterpolation(),
)
    extrapolation = Tuple(extrapolation)
    staggering = Tuple(staggering)

    num_src_dims = ndims(src_data)

    # Check dimensions match
    length(extrapolation) == num_src_dims || error(
        "Number of extrapolation conditions is not the same as the number of dimensions of the source data",
    )
    length(staggering) == num_src_dims || error(
        "Number of staggering is not the same as the number of dimensions of the source data",
    )
    length(src_grid) == num_src_dims || error(
        "Number of dimensions in src_grid is not the same as the number of dimensions of the source data",
    )
    map(length, src_grid) == size(src_data) || error(
        "Lengths of dimensions in src_grid ($(map(length, src_grid))) do not match size(src_data) ($(size(src_data)))",
    )

    for (dim_idx, coords) in enumerate(src_grid)
        length(coords) >= 2 || error(
            "Dimension $dim_idx has less than two coordinates, which is not supported yet",
        )
        # Strictly increasing rules out duplicate coordinates, which would be a
        # division by zero in the weights
        issorted(coords, lt = <=) ||
            error("Dimension $dim_idx is not strictly increasing")
    end

    # Check periodic boundary condition
    for (dim_idx, (coords, extp)) in enumerate(zip(src_grid, extrapolation))
        extp isa Periodic || continue
        extp.period > 0 || error(
            "Period of dimension $dim_idx ($(extp.period)) is not positive",
        )
        # A span larger than the period would silently alias wrapped queries
        # onto the wrong cells
        coords[end] - coords[begin] <= extp.period || error(
            "Span of dimension $dim_idx ($(coords[end] - coords[begin])) is larger than the period ($(extp.period))",
        )
    end

    # Create empty cache for storing results of stencil computation
    # nonmissingtype strips Missing from the eltype, since some datasets are
    # loaded with eltype Union{Missing, FT} even when no value is missing
    FT = float(promote_type(map(nonmissingtype ∘ eltype, src_grid)...))
    stencil_cache = ntuple(_ -> Tuple{Int, Int, FT}[], Val(num_src_dims))

    return Regridder(
        src_grid,
        src_data,
        extrapolation,
        staggering,
        method,
        stencil_cache,
    )
end

"""
    regrid(regridder::Regridder, dest_grid::OrderedDict)

Regrid the data in `regridder` to match `dest_grid`.
"""
function regrid(regridder::Regridder, dest_grid::OrderedDict)
    return regrid(regridder, Tuple(values(dest_grid)))
end

"""
    regrid!(dest_data, regridder::Regridder, dest_grid::OrderedDict)

Regrid `dest_data` in-place using `regridder` to match `dest_grid`.
"""
function regrid!(
    dest_data::AbstractArray,
    regridder::Regridder,
    dest_grid::OrderedDict,
)
    return regrid!(dest_data, regridder, Tuple(values(dest_grid)))
end

"""
    regrid(regridder::Regridder, dest_grid::Tuple)

Regrid the data in `regridder` to match `dest_grid`.
"""
function regrid(regridder::Regridder, dest_grid::Tuple)
    FT = float(nonmissingtype(eltype(regridder.src_data)))
    dest_data = Array{FT}(undef, map(length, dest_grid))
    regrid!(dest_data, regridder, dest_grid)
    return dest_data
end

"""
    regrid!(
        dest_data,
        regridder::Regridder{N, FT},
        dest_grid::Tuple,
    ) where {N, FT}

Regrid `dest_data` in-place using `regridder` to match `dest_grid`.
"""
function regrid!(
    dest_data::AbstractArray,
    regridder::Regridder{N, FT},
    dest_grid::Tuple,
) where {N, FT}
    length(dest_grid) == N || throw(
        ArgumentError(
            "Destination grid has $(length(dest_grid)) dimensions, but the regridder has $N dimensions for the source grid",
        ),
    )
    size(dest_data) == map(length, dest_grid) || throw(
        DimensionMismatch(
            "The number of dimensions of the destination data $(size(dest_data)) does not match the dimension sizes of destination grid ($(map(length, dest_grid)))",
        ),
    )

    # Precompute the cache for each dimension
    for (dim_idx, (src_coords, dest_coords)) in
        enumerate(zip(regridder.src_grid, dest_grid))
        _fill_dim_cache!(
            regridder.stencil_cache[dim_idx],
            dest_coords,
            src_coords,
            regridder.extrapolation[dim_idx],
            regridder.staggering[dim_idx],
        )
    end

    (; src_data, stencil_cache, method) = regridder
    T = float(nonmissingtype(eltype(dest_data)))
    @inbounds for J in CartesianIndices(dest_data)
        stencils = ntuple(d -> stencil_cache[d][J[d]], Val(N))
        dest_data[J] = _apply_stencils(T, src_data, stencils, method)
    end
    return dest_data
end

"""
    _fill_dim_cache!(
        cache,
        dest_coords::AbstractVector,
        src_coords::AbstractVector,
        extrapolation::AbstractExtrapolation,
        staggering::Union{Node, Center},
    )

Fill out the stencil `cache` for a single dimension.
"""
function _fill_dim_cache!(
    cache,
    dest_coords::AbstractVector,
    src_coords::AbstractVector,
    extrapolation::AbstractExtrapolation,
    staggering::Union{Node, Center},
)
    resize!(cache, length(dest_coords))
    for (j, x) in enumerate(dest_coords)
        cache[j] = _stencil_1d(x, src_coords, extrapolation, staggering)
    end
    return cache
end

"""
    _apply_stencils(
        ::Type{T},
        src_data::AbstractArray,
        stencils,
        method::AbstractInterpolationMethod,
    ) where {T}

Compute a single value for linear interpolation for a single point using the
`stencils` and `src_data`.
"""
@inline function _apply_stencils(
    ::Type{T},
    src_data::AbstractArray,
    stencils,
    ::LinearInterpolation,
) where {T}
    # length of a tuple constant-folds, so Val(N) is still a compile-time
    # constant even without annotating stencils
    N = length(stencils)
    # Iterate over all 2^N corners and compute a weighted average
    corners = CartesianIndices(ntuple(_ -> 2, Val(N)))
    acc = zero(T)
    @inbounds for corner in corners
        idx = ntuple(
            d -> corner[d] == 1 ? stencils[d][1] : stencils[d][2],
            Val(N),
        )
        w = one(stencils[1][3])
        for d in 1:N
            w *= corner[d] == 1 ? stencils[d][3] : one(w) - stencils[d][3]
        end
        acc += T(w * src_data[idx...])
    end
    return acc
end

# All stencil functions return (index1, index2, w) where w is the weight used to
# compute val = w * data[index1] + (1 - w) * data[index2]

"""
    _stencil_1d(
        x::Number,
        coords::AbstractVector,
        extrapolation::AbstractExtrapolation,
        staggering::Union{Node, Center},
    )

Compute the stencil `(index1, index2, w)` given `x` and `coords` with
`extrapolation` and `staggering`.
"""
function _stencil_1d(
    x::Number,
    coords::AbstractVector,
    extrapolation::AbstractExtrapolation,
    staggering::Union{Node, Center},
)
    FT = float(promote_type(typeof(x), nonmissingtype(eltype(coords))))
    (i1, i2, w) = if _in_interior(x, coords)
        _interior_stencil(x, coords)
    else
        _boundary_stencil(x, coords, extrapolation, staggering)
    end
    return (i1, i2, FT(w))
end

"""
    _in_interior(x::Number, coords::AbstractVector)

Return whether `x` is in the domain of `coords`.
"""
@inline function _in_interior(x::Number, coords::AbstractVector)
    first_coord = coords[begin]
    last_coord = coords[end]
    return first_coord <= x <= last_coord
end

"""
    _interior_stencil(x::Number, coords::AbstractVector)

Return the stencil `(index1, index2, w)` for `x` and `coords` when `x` is in
the domain of `coords`.
"""
@inline function _interior_stencil(x::Number, coords::AbstractVector)
    i = searchsortedfirst(coords, x)
    first_index = firstindex(coords)
    (i1, i2) = i != first_index ? (i - 1, i) : (first_index, first_index + 1)
    w = (coords[i2] - x) / (coords[i2] - coords[i1])
    return (i1, i2, w)
end

"""
    _boundary_stencil(
        x::Number,
        coords::AbstractVector,
        ::Throw,
        ::Node,
)

Throw an error when `x` is outside of `coords`.
"""
@inline function _boundary_stencil(
    x::Number,
    coords::AbstractVector,
    ::Throw,
    ::Node,
)
    first_coord = coords[begin]
    last_coord = coords[end]
    return throw(
        DomainError(
            x,
            "Cannot interpolate because $x is not in the domain [$first_coord, $last_coord]",
        ),
    )
end

"""
    _boundary_stencil(
        x::Number,
        coords::AbstractVector,
        ::Throw,
        ::Node,
)

Return the stencil `(index1, index2, w)` for `x` and `coords` when `x` is in the
center of the leftmost or rightmost nodes and throw an error otherwise.
"""
@inline function _boundary_stencil(
    x::Number,
    coords::AbstractVector,
    ::Throw,
    ::Center,
)
    left = firstindex(coords)
    right = lastindex(coords)
    half_width_left = (coords[left + 1] - coords[left]) / 2
    half_width_right = (coords[right] - coords[right - 1]) / 2
    if coords[left] - half_width_left <= x < coords[left]
        w = (coords[left + 1] - x) / (coords[left + 1] - coords[left])
        return (left, left + 1, w)
    elseif coords[right] < x <= coords[right] + half_width_right
        w = (coords[right] - x) / (coords[right] - coords[right - 1])
        return (right - 1, right, w)
    else
        throw(
            DomainError(
                x,
                "x is outside the domain [$(coords[left] - half_width_left), $(coords[right] + half_width_right)]",
            ),
        )
    end
end

"""
    _boundary_stencil(
        x::Number,
        coords::AbstractVector,
        ::Flat,
        ::Union{Center, Node},
    )

Return the stencil `(index1, index2, w)` for `x` and `coords` which compute the
same value at the leftmost or rightmost center or node.

Note that linear interpolation is not done for `Center`.
"""
function _boundary_stencil(
    x::Number,
    coords::AbstractVector,
    ::Flat,
    ::Union{Center, Node},
)
    left = firstindex(coords)
    right = lastindex(coords)
    w = one(float(nonmissingtype(eltype(coords))))
    return x < coords[left] ? (left, left, w) : (right, right, w)
end

"""
    _boundary_stencil(
        x::Number,
        coords::AbstractVector,
        extrapolation::Periodic,
        ::Union{Center, Node},
    )

Return the stencil `(index1, index2, w)` for `x` and `coords` which accounts
for `coords` being `Periodic`.
"""
function _boundary_stencil(
    x::Number,
    coords::AbstractVector,
    extrapolation::Periodic,
    ::Union{Center, Node},
)
    (; period) = extrapolation
    left = firstindex(coords)
    right = lastindex(coords)
    x_wrapped = coords[left] + mod(x - coords[left], period)
    if x_wrapped <= coords[right]
        return _interior_stencil(x_wrapped, coords)
    else
        # x_wrapped is in the wrap gap (coords[end], coords[begin] + period),
        # so interpolate between the last and first points
        wrap_coord = coords[left] + period
        w = (wrap_coord - x_wrapped) / (wrap_coord - coords[right])
        return (right, left, w)
    end
end

"""
    interpolate(
        regridder::Regridder{N, FT},
        point::Tuple{Vararg{Number, N}},
    ) where {N, FT}

Interpolate a single `point` using `regridder`.
"""
function interpolate(
    regridder::Regridder{N, FT},
    point::Tuple{Vararg{Number, N}},
) where {N, FT}
    # Since it is a single point, we do not need to cache the results
    stencils = map(
        (x, coords, extp, stag) -> begin
            (i1, i2, w) = _stencil_1d(x, coords, extp, stag)
            (i1, i2, FT(w))
        end,
        point,
        regridder.src_grid,
        regridder.extrapolation,
        regridder.staggering,
    )
    T = float(nonmissingtype(eltype(regridder.src_data)))
    return _apply_stencils(T, regridder.src_data, stencils, regridder.method)
end

"""
    interpolate(regridder::Regridder, points::Vector)

Interpolate `points` using `regridder`.
"""
function interpolate(regridder::Regridder, points::Vector)
    return map(point -> interpolate(regridder, point), points)
end

end
