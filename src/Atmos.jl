"""
The `Atmos` module contains functions that are primarily useful when working
with atmospheric simulations.
"""
module Atmos

import OrderedCollections: OrderedDict
import Interpolations as Intp
import ..OutputVar, ..arecompatible
import ..short_name
import ..altitude_name
import ..Var
import ..Numerics

export global_rmse_pfull

"""
    _resample_column!(dest, var1d, origin_pressure, target_pressure)

Linearly interpolate `var1d` from `origin_pressure` to `target_pressure`.

Note: Values outside of the range are linearly extrapolated
"""
function _resample_column!(dest, var1d, origin_pressure, target_pressure)
    reverse_origin_pressure = reverse(origin_pressure)
    isunique =
        length(unique(reverse_origin_pressure)) ==
        length(reverse_origin_pressure)
    !(isunique && issorted(reverse_origin_pressure)) &&
        error("P(z) is not bijective, cannot resample column")
    # Interpolations.jl require increasing knots, but pressure is decreasing, so
    # we have to reverse it
    var1d_of_P = Intp.extrapolate(
        Intp.interpolate(
            (reverse_origin_pressure,),
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
Plvl(z) = P0 * exp(-z / H_EARTH), where H_EARTH = 7000.0 and P0 = 1e5, following
a simple hydrostatic model for the atmosphere.

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

    # First, we construct the target pressure grid. For this, we use the mapping
    # Plvl(z) = P0 * exp(-z / H_EARTH), where H_EARTH = 7000.0 and P0 = 1e5.

    if isnothing(target_pressure)
        H_EARTH = 7000.0
        P0 = 1e5
        Plvl(z) = P0 * exp(-z / H_EARTH)

        # Reverse vector because Interpolations.jl require increasing knots
        target_pressure = reverse(Plvl.(var.dims[z_name]))
    end

    # Then, we prepare the output variable
    ret_attributes = copy(var.attributes)
    TypeOfDims = typeof(var.dims)
    ret_dims = TypeOfDims(
        k != z_name ? k => v : pressure_name => target_pressure for
        (k, v) in var.dims
    )

    TypeOfDimAttributes = OrderedDict{
        keytype(var.dim_attributes),
        Union{valtype(var.dim_attributes), typeof(pressure.attributes)},
    }
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

"""
    global_rmse_pfull(sim::OutputVar,
                      obs::OutputVar;
                      sim_pressure = nothing,
                      obs_pressure = nothing,
                      target_pressure = nothing)

Return the global RMSE in pressure space as a number.

The arguments `sim` and `obs` are defined as 3-dimensional `OutputVar`s with
the dimensions longitude, latitude, and altitude or pressure. If altitude is a
dimension, then `sim` and `obs` are converted to pressure coordinates using
[`to_pressure_coordinates`](@ref) through the keywords `sim_pressure` and
`obs_pressure`. Specific pressure levels can chosen by supplying a vector of
pressure levels for the keyword `target_pressure`.

Resampling is done automatically by resampling `obs` on `sim`.

The global RMSE is calculated by computing the mean squared error over longitude
and latitude, integrating and normalizing over pressure, and taking the square
root of the resulting number.
"""
function global_rmse_pfull(
    sim::OutputVar,
    obs::OutputVar;
    sim_pressure = nothing,
    obs_pressure = nothing,
    target_pressure = nothing,
)
    # Convert to pressure coordinates if necessary
    if Var.has_altitude(sim)
        isnothing(sim_pressure) && error(
            "sim_pressure is not supplied; cannot convert to pressure coordinates",
        )
        sim = to_pressure_coordinates(
            sim,
            sim_pressure,
            target_pressure = target_pressure,
        )
    end

    if Var.has_altitude(obs)
        isnothing(obs_pressure) && error(
            "obs_pressure is not supplied; cannot convert to pressure coordinates",
        )
        obs = to_pressure_coordinates(
            obs,
            obs_pressure,
            target_pressure = target_pressure,
        )
    end
    # Check if presure is a dimension
    Var.has_pressure(sim) || error("sim does not have pressure as a dimension")
    Var.has_pressure(obs) || error("obs does not have pressure as a dimension")

    # Follow the same template for squared_error
    Var._check_sim_obs_units_consistent(sim, obs, 3)
    obs_resampled = Var.resampled_as(obs, sim)

    squared_error = (sim - obs_resampled) * (sim - obs_resampled)

    # Compute global mse and global rmse and store it as an attribute
    integrated_squared_error_var = Var.integrate_lonlat(squared_error)
    integrated_squared_error = integrated_squared_error_var.data
    ones_var = OutputVar(
        squared_error.attributes,
        squared_error.dims,
        squared_error.dim_attributes,
        ones(size(squared_error.data)),
    )
    normalization = Var.integrate_lonlat(ones_var).data
    # Doing ./ because both quantities are vectors
    mse = integrated_squared_error ./ normalization

    pfull = Var.pressures(integrated_squared_error_var)
    # Use first because the return type is an array
    length_pfull =
        first(Numerics._integrate_dim(ones(size(pfull)), pfull; dims = 1))

    # Normalize and compute RMSE
    # Use first because the return type is an array
    return sqrt(
        (1 / length_pfull) *
        first(Numerics._integrate_dim(mse, pfull; dims = 1)),
    )
end

end
