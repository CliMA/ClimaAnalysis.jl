module Numerics

import ..Utils: _isequispaced
import NaNStatistics: nansum

"""
    _integrate_lon(data::AbstractArray, lon::AbstractVector; dims)

Integrate out longitude from `data`. `data` has to be discretized on `lon`.

`dims` indicates which axis of `data` is longitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
longitude should be degrees.
"""
function _integrate_lon(data::AbstractArray, lon::AbstractVector; dims)
    length(lon) == 1 &&
        error("Cannot integrate when longitude is a single point")
    _isequispaced(lon) ?
    int_weights = _integration_weights_lon_equispaced(lon) :
    int_weights = _integration_weights_lon_left(lon)
    return _integrate_over_generic_dim(data, lon, dims, int_weights)
end

"""
    _integrate_lat(data::AbstractArray, lat::AbstractVector; dims)

Integrate out latitude from `data`. `data` has to be discretized on `lat`.

`dims` indicates which axis of `data` is latitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
latitude should be degrees.
"""
function _integrate_lat(data::AbstractArray, lat::AbstractVector; dims)
    length(lat) == 1 &&
        error("Cannot integrate when latitude is a single point")
    _isequispaced(lat) ?
    int_weights = _integration_weights_lat_equispaced(lat) :
    int_weights = _integration_weights_lat_left(lat)
    return _integrate_over_generic_dim(data, lat, dims, int_weights)
end

"""
    _integrate_dim(data::AbstractArray, dim::AbstractVector; dims)

Integrate out a dimension `data`. `data` has to be discretized on `dim`.

`dims` indicates which axis of `data` is `dim`.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints.
"""
function _integrate_dim(data::AbstractArray, dim::AbstractVector; dims)
    length(dim) == 1 && error("Cannot integrate when z is a single point")
    _isequispaced(dim) ?
    int_weights = _integration_weights_generic_equispaced(dim) :
    int_weights = _integration_weights_generic_left(dim)
    return _integrate_over_generic_dim(data, dim, dims, int_weights)
end

"""
    _integrate_over_generic_dim(
    data::AbstractArray,
    angle_arr::AbstractVector,
    angle_idx,
    int_weights::AbstractVector,
)

Integrate out a generic dimension from `data` using the weights `int_weights`.
`data` has to be discretized on `angle_arr`.

`angle_idx` indicates which axis of `data` is angle.
"""
function _integrate_over_generic_dim(
    data::AbstractArray,
    angle_arr::AbstractVector,
    angle_idx,
    int_weights,
)
    # Reshape to add extra dimensions for int_weights for broadcasting if needed
    size_to_reshape =
        (i == angle_idx ? length(int_weights) : 1 for i in 1:ndims(data))
    int_weights = reshape(int_weights, size_to_reshape...)
    int_on_angle = nansum(data .* int_weights, dims = angle_idx)
    return int_on_angle
end

"""
    _integration_weights_lon_left(lon)

Return integration weights for rectangular integration using left endpoints for integrating
along longitude.
"""
function _integration_weights_lon_left(lon)
    # This is where we use the assumption that units are degrees
    d_lon = deg2rad.(diff(lon))
    # We are doing integration using the left endpoints, so we weight the rightmost endpoint
    # zero so that it make no contribution to the integral
    push!(d_lon, zero(eltype(d_lon)))
    return d_lon
end

"""
    _integration_weights_lat_left(lat)

Return integration weights for rectangular integration using left endpoints for integrating
along latitude.
"""
function _integration_weights_lat_left(lat)
    d_lat = deg2rad.(diff(lat))
    # We are doing integration using the left endpoints, so we weight the rightmost endpoint
    # zero so that it make no contribution to the integral
    push!(d_lat, zero(eltype(d_lat)))
    cos_lat = cosd.(lat)
    return d_lat .* cos_lat
end

"""
    _integration_weights_lon_equispaced(lon)

Return integration weights for rectangular integration when the points are equispaced for
integrating along longitude.
"""
function _integration_weights_lon_equispaced(lon)
    # This is where we use the assumption that units are degrees
    # Use fill to make a zero dimensional array so reshaping is possible
    return fill(deg2rad(lon[begin + 1] - lon[begin]))
end

"""
    _integration_weights_lat_equispaced(lat)

Return integration weights for rectangular integration when the points are equispaced for
integrating along latitude.
"""
function _integration_weights_lat_equispaced(lat)
    d_lat = deg2rad.(lat[begin + 1] - lat[begin])
    cos_lat = cosd.(lat)
    return d_lat .* cos_lat
end

function _integration_weights_generic_left(arr)
    # This is where we use the assumption that units are degrees
    d_arr = diff(arr)
    # We are doing integration using the left endpoints, so we weight the rightmost endpoint
    # zero so that it make no contribution to the integral
    push!(d_arr, zero(eltype(d_arr)))
    return d_arr
end

"""
    _integration_weights_generic_equispaced(arr)

Return integration weights for rectangular integration when the points are equispaced for
integrating along a generic dimension.
"""
function _integration_weights_generic_equispaced(arr)
    # This is where we use the assumption that units are degrees
    # Use fill to make a zero dimensional array so reshaping is possible
    return fill(arr[begin + 1] - arr[begin])
end

#=
Why not use Interpolations.jl?

One of the most expensive operation in terms of time and memory is
`Var.resampled_as(src_var, dest_var)` which uses an interpolation to resample `data` in
`src_var` to match the `dims` in `dest_var`. However, interpolating using Interpolations.jl
is expensive. One must first intialize an interpolant, whose size in memory is at least as
big as the memory of `src_var.data`. Then, each evaluation using the interpolant
costs some amount of allocation on the heap. As a result, interpolating is extremely slow
for a large number of points.

To solve this, we write our interpolation routine in the Numerics module that support the following:
1. No extra dependencies
2. No allocation when evaluating a point
3. Support boundary conditions for periodic, flat, and throw on an irregular grid
5. Comparable or better performance to Interpolations.jl
=#

"""
    linear_interpolate(point::NTuple{N, FT1},
                       axes::NTuple{N, Vector},
                       data::AbstractArray{FT2, N},
                       extp_conds::NTuple{N, NTuple{4, Function}}) where {N, FT1, FT2}

Linear interpolate `data` on `axes` and return the value at `point`. Extrapolation is
handled by `extp_conds`.
"""
function linear_interpolate(
    point::NTuple{N},
    axes,
    data::AbstractArray{FT2, N},
    extp_conds,
) where {N, FT2}
    # Get a new point as determined by the extrapolation condition
    point = extp_to_point(point, axes, extp_conds)

    # Find which cell contain the point
    cell_indices_for_axes = find_cell_indices_for_axes(point, axes)
    val = zero(FT2)

    # Compute the denominator of the formula for linear interpolation
    # (in 1D, this is x_1 - x_0)
    bottom_term = compute_bottom_term(axes, cell_indices_for_axes)

    # Iterate through all 2^N points
    @inbounds for bits in 0:(2^N - 1)
        term = one(FT2)
        bound_indices = get_indices(cell_indices_for_axes, bits)
        sign = get_sign(cell_indices_for_axes, bits)
        # Weight is the value at each of the points of the cell
        weight = data[get_complement_indices(cell_indices_for_axes, bits)...]
        @inbounds for (dim_idx, bound_idx) in pairs(bound_indices)
            val_minus_x2_or_x1 = point[dim_idx] - axes[dim_idx][bound_idx]
            term *= val_minus_x2_or_x1
        end
        term *= sign * weight
        val += term
    end
    return val / bottom_term
end

function linear_interpolate(
    point::Number,
    axes,
    data::AbstractArray{FT2, N},
    extp_conds,
) where {N, FT2}
    point = Tuple(point...)
    return linear_interpolate(point, axes, data, extp_conds)
end

function linear_interpolate(
    point::AbstractVector,
    axes,
    data::AbstractArray{FT2, N},
    extp_conds,
) where {N, FT2}
    point = Tuple(coord for coord in point)
    return linear_interpolate(point, axes, data, extp_conds)
end

function linear_interpolate(
    point::Tuple,
    axes,
    data::AbstractArray{FT2, N},
    extp_conds,
) where {N, FT2}
    point = promote(point...)
    return linear_interpolate(point, axes, data, extp_conds)
end


"""
    compute_bottom_term(axes, cell_indices_for_axes::NTuple{N})

Compute the bottom term when computing in linear interpolation.

Consider the formula for 1D linear interpolation which is
    y_0 * ((x_1 - x) / (x_1 - x_0)) + y_1 * ((x - x_0) / (x_1 - x_0))
for interpolating the point (x, y) on line between (x_0, y_0) and (x_1, y_1). This function
computes x_1 - x_0.
"""
function compute_bottom_term(axes, cell_indices_for_axes::NTuple{N}) where {N}
    return reduce(
        *,
        ntuple(
            dim_idx ->
                axes[dim_idx][cell_indices_for_axes[dim_idx][end]] -
                axes[dim_idx][cell_indices_for_axes[dim_idx][begin]],
            N,
        ),
    )
end

"""
    get_complement_indices(indices::NTuple{N, Tuple{I, I}}, bits) where {N, I}

Given a tuple consisting of 2-tuple, return a tuple of one element from each tuple according
to bits. The elements in the tuple are the complement of those found by `get_indices`.
"""
function get_complement_indices(indices::NTuple{N}, bits) where {N}
    return ntuple(dim -> if (bits & (1 << (dim - 1)) != 0)
        indices[dim][2]
    else
        indices[dim][1]
    end, N)
end

"""
    get_indices(indices::NTuple{N, Tuple{I, I}}, bits) where {N, I}

Given a tuple consisting of 2-tuple, return a tuple of one element from each tuple according to bits.
"""
function get_indices(indices::NTuple{N}, bits) where {N}
    return ntuple(dim -> if (bits & (1 << (dim - 1)) != 0)
        indices[dim][1]
    else
        indices[dim][2]
    end, N)
end

"""
    get_sign(indices::NTuple{N, Tuple{I, I}}, bits) where {N, I}

Given a tuple consisting of 2-tuple, compute the appropriate sign when interpolating.
"""
function get_sign(_indices::NTuple{N}, bits) where {N}
    return reduce(*, ntuple(dim -> if (bits & (1 << (dim - 1)) != 0)
        1
    else
        -1
    end, N))
end

"""
    find_cell_indices_for_axes(point::NTuple{N, FT},
                               axes::NTuple{N, A}) where {N, FT, A <:AbstractVector}

Given a point and axes, find the indices of the N-dimensional hyperrectangle, where the
point lives in.
"""
function find_cell_indices_for_axes(point::NTuple{N}, axes) where {N}
    return ntuple(
        dim_idx -> find_cell_indices_for_ax(point[dim_idx], axes[dim_idx]),
        N,
    )
end

"""
    find_cell_indices_for_ax(val::FT1, ax::AbstractVector{FT2}) where {FT1, FT2}

Given `val` and an `ax`, find the indices of the cell, where `val` lives in.
"""
function find_cell_indices_for_ax(
    val::FT1,
    ax::AbstractVector{FT2},
) where {FT1, FT2}
    len_of_ax = length(ax)
    (val == ax[begin]) && return (1, 2)
    (val == ax[end]) && return (len_of_ax - 1, len_of_ax)
    idx = searchsortedfirst(ax, val)
    return (idx - 1, idx)
end

"""
    extp_to_point(point::NTuple{N, FT1}, axes::NTuple{N, Vector}, extp_conds) where {N, FT1}

Return a new point to evaluate at according to the extrapolation conditions.
"""
function extp_to_point(point::NTuple{N}, axes, extp_conds) where {N}
    return ntuple(
        idx -> extp_conds[idx].get_val_for_point(point[idx], axes[idx]),
        N,
    )
end

"""
    extp_cond_throw()

Return extrapolation condition for throwing an error when the point is out of bounds.

The first and last nodes are not co-located. For example, if the axis is [1.0, 2.0, 3.0]
and the data is [4.0, 5.0, 6.0], then the value at 3.0 is 6.0 and not 4.0.
"""
function extp_cond_throw()
    get_val_for_point(val, ax) = begin
        val < ax[begin] && return error("Out of bounds error with $val in $ax")
        val > ax[end] && return error("Out of bounds error with $val in $ax")
        return val
    end
    return (; get_val_for_point = get_val_for_point)
end

"""
    extp_cond_flat()

Return flat extrapolation condition.
"""
function extp_cond_flat()
    get_val_for_point(val, ax) = begin
        val < ax[begin] && return typeof(val)(ax[begin])
        val > ax[end] && return typeof(val)(ax[end])
        return val
    end
    return (; get_val_for_point)
end

"""
    extp_cond_periodic()

Return periodic extrapolation condtion.
"""
function extp_cond_periodic()
    get_val_for_point(val, ax) = begin
        if (val < ax[begin]) || (val > ax[end])
            width = ax[end] - ax[begin]
            new_val = mod(val - ax[begin], width) + ax[begin]
            return typeof(val)(new_val)
        end
        return val
    end
    return (; get_val_for_point)
end

end
