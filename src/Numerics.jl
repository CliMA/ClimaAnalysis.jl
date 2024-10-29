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

end
