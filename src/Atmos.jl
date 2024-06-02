"""
The `Atmos` module contains functions that are primarily useful when working
with atmospheric simulations.
"""
module Atmos

import Interpolations as Intp
import ..OutputVar
import ..short_name

"""

    TODO
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
        Intp.Linear(),
    )

    dest .= var1d_of_P.(target_pressure)
    return nothing
end

"""

TODO
"""
function to_pressure_coordinates(var::OutputVar, pressure::OutputVar)
    # TODO: Add error checking:
    # - `var` and `pressure` have to be compatible

    ALTITUDE_NAMES = Set(["z", "z_physical", "z_reference"])

    # Pick the correct longitude name and check that we have an altitude variable
    z_name = ""
    for possible_z_name in ALTITUDE_NAMES
        haskey(var.dims, possible_z_name) && (z_name = possible_z_name; break)
    end
    z_name != "" || error("var does not have altitude among its dimensions")
    z_index = var.dim2index[z_name]
    pressure_name = short_name(pressure)

    # First, we construct the target pressure grid. For this, we take the
    # extrema of pressure and divide the interval linearly with the same number
    # of points we originally had in z

    # TODO: Pick this more sensibly
    # TODO: Make it go from max to min? (This is not supported by Interpolations.jl...)
    target_pressure = range(
        minimum(pressure.data),
        maximum(pressure.data),
        length = length(var.dims[z_name]),
    )

    # Then, we prepare the output variable
    ret_attributes = copy(var.attributes)
    TypeOfDims = typeof(var.dims)
    ret_dims = TypeOfDims(
        k != z_name ? k => v : pressure_name => target_pressure for
        (k, v) in var.dims
    )
    TypeOfDimAttributes = typeof(var.dim_attributes)
    ret_dim_attributes = TypeOfDims(
        k != z_name ? k => v : pressure_name => pressure.attributes for
        (k, v) in var.dim_attributes
    )

    ret_data = zeros(size(var.data)...)

    # We have to loop over all the possible columns
    num_dims = ndims(var.data)
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
