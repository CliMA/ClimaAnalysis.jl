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
    NaNLinearInterpolation(threshold = 0.5)

Linear interpolation that is aware of NaNs in the source data. The `threshold`
must be between zero and one.

For each interpolated point, the interpolation weights of the NaN values are
summed. If this sum is greater than the `threshold`, the result is `NaN`.
Otherwise, the result is the weighted sum of the non-NaN values renormalized by
the sum of their weights.

A `threshold` of zero means any `NaN` with a nonzero weight in the interpolation
produces `NaN`, and a `threshold` of one means the result is `NaN` only when all
the values with nonzero weights in the interpolation are `NaN`.

A `missing` value anywhere in the interpolation produces `missing`, regardless
of its weight or the `threshold`.
"""
struct NaNLinearInterpolation{FT <: AbstractFloat} <:
       AbstractInterpolationMethod
    threshold::FT
    function NaNLinearInterpolation(threshold::Real = 0.5)
        0 <= threshold <= 1 ||
            error("Threshold ($threshold) must be between 0 and 1")
        threshold = float(threshold)
        return new{typeof(threshold)}(threshold)
    end
end

"""
    Regridder

Linearly interpolate gridded data with support for node and centers and throw,
flat, and periodic boundary conditions.
"""
struct Regridder{
    N,
    A <: AbstractArray,
    SG <: NTuple{N, AbstractVector},
    E <: NTuple{N, AbstractExtrapolation},
    S <: NTuple{N, Union{Node, Center}},
    M <: AbstractInterpolationMethod,
}

    "Source data to interpolate from"
    src_data::A

    "Coordinates of the source data"
    src_grid::SG

    "Extrapolation conditions for how queries outside the source coordinates are
    handled"
    extrapolation::E

    "Whether the coordinates of each dimension are cell centers (`Center`) or
    cell edges (`Node`)."
    staggering::S

    "Interpolation method (`LinearInterpolation` or `NaNLinearInterpolation`)"
    method::M
end

"""
    Regridder(
        src_data::AbstractArray,
        src_grid,
        extrapolation,
        staggering;
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

The keyword `method` (`LinearInterpolation` or `NaNLinearInterpolation`)
determines how NaNs in `src_data` are handled.
"""
function Regridder(
    src_data::AbstractArray,
    src_grid,
    extrapolation,
    staggering;
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
            "Dimension $dim_idx has less than two coordinates, which is not supported",
        )

        # NaN and missing values in coordinates are not supported
        any(x -> ismissing(x) || isnan(x), coords) && error(
            "Missing or NaN values for dimension $dim_idx is not supported",
        )

        issorted(coords, lt = <=) ||
            error("Dimension $dim_idx is not strictly increasing")
    end

    for (dim_idx, (coords, extp)) in enumerate(zip(src_grid, extrapolation))
        # Check the periodic boundary condition if it exists for each dimension
        extp isa Periodic || continue
        extp.period > 0 || error(
            "Period of dimension $dim_idx ($(extp.period)) is not positive",
        )
        span = coords[end] - coords[begin]
        (span <= extp.period || span ≈ extp.period) || error(
            "Span of dimension $dim_idx ($span) is larger than the period ($(extp.period))",
        )
    end

    return Regridder(src_data, src_grid, extrapolation, staggering, method)
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
    dest_data = Array{_dest_eltype(regridder)}(undef, map(length, dest_grid))
    regrid!(dest_data, regridder, dest_grid)
    return dest_data
end

"""
    regrid!(
        dest_data,
        regridder::Regridder{N},
        dest_grid::Tuple,
    ) where {N}

Regrid `dest_data` in-place using `regridder` to match `dest_grid`.
"""
function regrid!(
    dest_data::AbstractArray,
    regridder::Regridder{N},
    dest_grid::Tuple,
) where {N}
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

    # Create empty cache for storing results of stencil computation
    # The cache are preallocated per-dimension buffers of 1D stencils
    # `(index1, index2, weight)`, one per destination coordinate
    # The weight type is promoted from the element types of the destination
    # grid, so the precision of the destination coordinates is not lost when
    # computing the weights
    FT = float(
        promote_type(
            _src_eltype(regridder),
            map(nonmissingtype ∘ eltype, dest_grid)...,
        ),
    )
    stencil_cache = ntuple(_ -> Tuple{Int, Int, FT}[], Val(N))

    # Precompute the cache for each dimension
    for (dim_idx, (src_coords, dest_coords)) in
        enumerate(zip(regridder.src_grid, dest_grid))
        _fill_dim_cache!(
            stencil_cache[dim_idx],
            dest_coords,
            src_coords,
            regridder.extrapolation[dim_idx],
            regridder.staggering[dim_idx],
        )
    end

    (; src_data, method) = regridder
    @inbounds for J in CartesianIndices(dest_data)
        stencils = ntuple(d -> stencil_cache[d][J[d]], Val(N))
        dest_data[J] = _apply_stencils(FT, src_data, stencils, method)
    end
    return nothing
end

"""
    _fill_dim_cache!(
        cache::Vector{Tuple{Int, Int, FT}},
        dest_coords::AbstractVector,
        src_coords::AbstractVector,
        extrapolation::AbstractExtrapolation,
        staggering::Union{Node, Center},
    )

Fill out the stencil `cache` for a single dimension.

Each destination coordinate is converted to the weight type `FT`, so the
stencil weights are computed at the precision of `FT`.
"""
function _fill_dim_cache!(
    cache::Vector{Tuple{Int, Int, FT}},
    dest_coords::AbstractVector,
    src_coords::AbstractVector,
    extrapolation::AbstractExtrapolation,
    staggering::Union{Node, Center},
) where {FT}
    resize!(cache, length(dest_coords))
    for (j, x) in enumerate(dest_coords)
        cache[j] = _stencil_1d(FT(x), src_coords, extrapolation, staggering)
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
        v = src_data[idx...]
        # A missing value propagates through the weighted sum
        ismissing(v) && return missing
        acc += w * v
    end
    # This is needed for the interpolate and interpolate! to return the same
    # type as regrid and regrid!
    return T(acc)
end

"""
    _apply_stencils(
        ::Type{T},
        src_data::AbstractArray,
        stencils,
        method::NaNLinearInterpolation,
    ) where {T}

Compute a single value for `NaN`-aware linear interpolation for a single point using the
`stencils` and `src_data`.

For more details, see `NaNLinearInterpolation`.
"""
@inline function _apply_stencils(
    ::Type{T},
    src_data::AbstractArray,
    stencils,
    method::NaNLinearInterpolation,
) where {T}
    N = length(stencils)
    # Iterate over all 2^N corners and compute a weighted average
    corners = CartesianIndices(ntuple(_ -> 2, Val(N)))
    acc = zero(T)
    nan_weight = zero(T)
    @inbounds for corner in corners
        idx = ntuple(
            d -> corner[d] == 1 ? stencils[d][1] : stencils[d][2],
            Val(N),
        )
        w = one(stencils[1][3])
        for d in 1:N
            w *= corner[d] == 1 ? stencils[d][3] : one(w) - stencils[d][3]
        end
        v = src_data[idx...]
        # A missing value propagates through the weighted sum
        ismissing(v) && return missing
        v_isnan = isnan(v)
        nan_weight += ifelse(v_isnan, w, zero(w))
        acc += w * ifelse(v_isnan, zero(v), v)
    end
    nan_weight > method.threshold && return T(NaN)
    # The corner weights sum to one, so the sum of the weights of the non-NaN
    # values is 1 - nan_weight (NaN if all values are NaN, from dividing by 0)
    return T(acc / (one(T) - nan_weight))
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
        ::Center,
)

Return the stencil `(index1, index2, w)` for `x` and `coords` when `x` is within
half a cell of the first or last coordinate and throw an error otherwise.

The stencil linearly extends the two nearest coordinates, so the weight is
outside of `[0, 1]` there.
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
    in_left_cell = coords[left] - half_width_left <= x < coords[left]
    if in_left_cell
        w = (coords[left + 1] - x) / (coords[left + 1] - coords[left])
        return (left, left + 1, w)
    end
    in_right_cell = coords[right] < x <= coords[right] + half_width_right
    if in_right_cell
        w = (coords[right] - x) / (coords[right] - coords[right - 1])
        return (right - 1, right, w)
    end
    throw(
        DomainError(
            x,
            "x is outside the domain [$(coords[left] - half_width_left), $(coords[right] + half_width_right)]",
        ),
    )
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

For `Center`, the value is constant everywhere outside the first and last
coordinates, including within half a cell of them.
"""
@inline function _boundary_stencil(
    x::Number,
    coords::AbstractVector,
    ::Flat,
    ::Union{Center, Node},
)
    left = firstindex(coords)
    right = lastindex(coords)
    FT = float(nonmissingtype(eltype(coords)))
    # NaN compares false with <
    isnan(x) && return (left, left, FT(NaN))
    w = one(FT)
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
@inline function _boundary_stencil(
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
        regridder::Regridder{N},
        point::Tuple{Vararg{Number, N}},
    ) where {N}

Interpolate a single `point` using the `regridder`.
"""
function interpolate(
    regridder::Regridder{N},
    point::Tuple{Vararg{Number, N}},
) where {N}
    FT = _src_eltype(regridder)
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
    return _apply_stencils(FT, regridder.src_data, stencils, regridder.method)
end

"""
    interpolate(regridder::Regridder, points::AbstractArray)

Interpolate multiple `points` using the `regridder`.
"""
function interpolate(regridder::Regridder, points::AbstractArray)
    return map(point -> interpolate(regridder, point), points)
end

"""
    _src_eltype(regridder::Regridder)

Return the float type promoted from the element types of the source data and
source grid of `regridder`, excluding `Missing`.
"""
function _src_eltype(regridder::Regridder)
    return float(
        promote_type(
            nonmissingtype(eltype(regridder.src_data)),
            map(nonmissingtype ∘ eltype, regridder.src_grid)...,
        ),
    )
end

"""
    _dest_eltype(regridder::Regridder)

Return the eltype of the data produced by `regridder`.
"""
function _dest_eltype(regridder::Regridder)
    return promote_type(
        _src_eltype(regridder),
        float(eltype(regridder.src_data)),
    )
end

end
