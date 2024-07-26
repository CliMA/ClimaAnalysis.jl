"""
The `Atmos` module contains functions that are primarily useful when working
with atmospheric simulations.
"""
module Atmos

import Interpolations as Intp
import ..OutputVar, ..arecompatible
import ..short_name
import ..altitude_name

"""
    _resample_column!(dest, var1d, origin_pressure, target_pressure)

Linearly interpolate `var1d` from `origin_pressure` to `target_pressure`.

Note: Values outside of the range are linearly extrapolated
"""
function _resample_column!(dest, var1d, origin_pressure, target_pressure)
    # Interpolations.jl require increasing knots, but pressure is decreasing, so
    # we have to reverse it
    var1d_of_P = Intp.extrapolate(
        Intp.interpolate(
            (reverse(origin_pressure),),
            reverse(var1d),
            Intp.Gridded(Intp.Linear()),
        ),
        Intp.Line(),
    )

    dest .= var1d_of_P.(target_pressure)
    return nothing
end

"""
    to_pressure_coordinates(var::OutputVar, pressure::OutputVar; target_pressure = nothing)

Change the vertical dimension of `var` to be in pressure coordinates.

If `target_pressure` is nothing, the target pressure levels are computed by
linearly sampling the interval `minimum(pressure), maximum(pressure)`. Then, for
each column in `var`, the values are linearly interpolate onto this new grid.

`target_pressure` can be set to a `Vector` to specify custom pressure levels.

The return value is a new `OutputVar` where the vertical dimension is pressure.

> :important: Values outside of the range are linearly extrapolated, so do not
  trust them too much!

"""
function to_pressure_coordinates(
    var::OutputVar,
    pressure::OutputVar;
    target_pressure = nothing,
)
    arecompatible(var, pressure) ||
        error("Pressure and variable are not compatible")

    z_name = altitude_name(var)
    z_index = var.dim2index[z_name]
    pressure_name = short_name(pressure)

    # First, we construct the target pressure grid. For this, we take the
    # extrema of pressure and divide the interval linearly with the same number
    # of points we originally had in z

    if isnothing(target_pressure)
        # TODO: Pick this more sensibly
        # TODO: Make it go from max to min? (This is not supported by Interpolations.jl...)
        target_pressure = range(
            minimum(pressure.data),
            maximum(pressure.data),
            length = length(var.dims[z_name]),
        )
    end

    # Then, we prepare the output variable
    ret_attributes = copy(var.attributes)
    TypeOfDims = typeof(var.dims)
    ret_dims = TypeOfDims(
        k != z_name ? k => v : pressure_name => target_pressure for
        (k, v) in var.dims
    )
    TypeOfDimAttributes = typeof(var.dim_attributes)
    ret_dim_attributes = TypeOfDimAttributes(
        k != z_name ? k => v : pressure_name => pressure.attributes for
        (k, v) in var.dim_attributes
    )

    num_dims = ndims(var.data)

    # Account for possible custom target_pressure
    ret_size = ntuple(
        i -> (i != z_index ? size(var.data)[i] : length(target_pressure)),
        num_dims,
    )
    ret_data = zeros(ret_size...)

    # We have to loop over all the possible columns
    ranges = [1:size(var.data)[i] for i in 1:num_dims if i != z_index]

    # Iterate over the Cartesian product of these ranges
    for idx in CartesianIndices(Tuple(ranges))
        indices = ntuple(i -> (i == z_index ? Colon() : idx[i]), num_dims)

        _resample_column!(
            view(ret_data, indices...),
            view(var.data, indices...),
            view(pressure.data, indices...),
            target_pressure,
        )
    end

    return OutputVar(ret_attributes, ret_dims, ret_dim_attributes, ret_data)
end

end
