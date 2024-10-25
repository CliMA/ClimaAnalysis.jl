module Var

import Dates
import NCDatasets
import OrderedCollections: OrderedDict

import Interpolations as Intp
import Statistics: mean
import NaNStatistics: nanmean

import ..Numerics
import ..Utils:
    nearest_index,
    seconds_to_prettystr,
    squeeze,
    split_by_season,
    time_to_date,
    date_to_time,
    _data_at_dim_vals,
    _isequispaced

export OutputVar,
    read_var,
    average_lat,
    weighted_average_lat,
    average_lon,
    average_x,
    average_y,
    average_xy,
    average_time,
    is_z_1D,
    slice,
    window,
    arecompatible,
    center_longitude!,
    short_name,
    long_name,
    units,
    dim_units,
    range_dim,
    reordered_as,
    resampled_as,
    has_units,
    convert_units,
    integrate_lonlat,
    integrate_lon,
    integrate_lat,
    isempty,
    split_by_season,
    bias,
    global_bias,
    squared_error,
    global_mse,
    global_rmse,
    set_units,
    set_dim_units!,
    shift_to_start_of_previous_month,
    apply_landmask,
    apply_oceanmask,
    make_lonlat_mask,
    replace

"""
    Representing an output variable
"""
struct OutputVar{T <: AbstractArray, A <: AbstractArray, B, C, ITP}

    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Array that contains all the data"
    data::A

    "Dictionary that maps dimension name to its array index"
    dim2index::Dict{String, Int}

    "Array that maps name array index to the dimension name"
    index2dim::Vector{String}

    "Interpolant from Interpolations.jl, used to evaluate the OutputVar onto any given point."
    interpolant::ITP
end

"""
    _make_interpolant(dims, data)

Make a linear interpolant from `dims`, a dictionary mapping dimension name to array and
`data`, an array containing data. Used in constructing a `OutputVar`.

If any element of the arrays in `dims` is a Dates.DateTime, then no interpolant is returned.
Interpolations.jl does not support interpolating on dates. If the longitudes span the entire
range and are equispaced, then a periodic boundary condition is added for the longitude
dimension. If the latitudes span the entire range and are equispaced, then a flat boundary
condition is added for the latitude dimension. In all other cases, an error is thrown when
extrapolating outside of `dim_array`.
"""
function _make_interpolant(dims, data)
    # If any element is DateTime, then return nothing for the interpolant because
    # Interpolations.jl do not support DateTimes
    for dim_array in values(dims)
        eltype(dim_array) <: Dates.DateTime && return nothing
    end

    # We can only create interpolants when we have 1D dimensions
    if isempty(dims) || any(d -> ndims(d) != 1 || length(d) == 1, values(dims))
        return nothing
    end

    # Dimensions are all 1D, check that they are compatible with data
    size_data = size(data)
    for (dim_index, (dim_name, dim_array)) in enumerate(dims)
        dim_length = length(dim_array)
        data_length = size_data[dim_index]
        if dim_length != data_length
            error(
                "Dimension $dim_name has inconsistent size with provided data ($dim_length != $data_length)",
            )
        end
    end

    # Find boundary conditions for extrapolation
    extp_bound_conds = (
        _find_extp_bound_cond(dim_name, dim_array) for
        (dim_name, dim_array) in dims
    )

    dims_tuple = tuple(values(dims)...)
    extp_bound_conds_tuple = tuple(extp_bound_conds...)
    return Intp.extrapolate(
        Intp.interpolate(dims_tuple, data, Intp.Gridded(Intp.Linear())),
        extp_bound_conds_tuple,
    )
end

"""
    _find_extp_bound_cond(dim_name, dim_array)

Find the appropriate boundary condition for the `dim_name` dimension.
"""
function _find_extp_bound_cond(dim_name, dim_array)
    min_of_dim, max_of_dim = extrema(dim_array)
    dim_size = max_of_dim - min_of_dim
    dsize = dim_array[begin + 1] - dim_array[begin]

    # If the dimension array span the entire range and is equispaced, then add the
    # appropriate boundary condition
    # We do not handle the cases when the array is not equispaced
    (
        conventional_dim_name(dim_name) == "longitude" &&
        _isequispaced(dim_array) &&
        isapprox(dim_size + dsize, 360.0)
    ) && return Intp.Periodic()
    (
        conventional_dim_name(dim_name) == "latitude" &&
        _isequispaced(dim_array) &&
        isapprox(dim_size + dsize, 180.0)
    ) && return Intp.Flat()
    return Intp.Throw()
end

function OutputVar(attribs, dims, dim_attribs, data)
    index2dim = keys(dims) |> collect
    dim2index =
        Dict([dim_name => index for (index, dim_name) in enumerate(keys(dims))])

    # TODO: Make this lazy: we should compute the spline the first time we use
    # it, not when we create the object
    itp = _make_interpolant(dims, data)

    function _maybe_process_key_value(k, v)
        k != "units" && return k => v
        return k => _maybe_convert_to_unitful(v)
    end

    # Recreating an object to ensure that the type is correct
    if !isempty(attribs)
        attribs = Dict(_maybe_process_key_value(k, v) for (k, v) in attribs)
    end

    # TODO: Support units for dimensions too

    return OutputVar(
        attribs,
        OrderedDict(dims),
        OrderedDict(dim_attribs),
        data,
        dim2index,
        index2dim,
        itp,
    )
end

function OutputVar(dims, data)
    return OutputVar(Dict{String, Any}(), dims, Dict{String, Any}(), data)
end

"""
    OutputVar(path,
              short_name = nothing;
              new_start_date = nothing,
              shift_by = identity)

Read the NetCDF file in `path` as an `OutputVar`.

If `short_name` is `nothing`, automatically find the name.

Dates in the time dimension are automatically converted to seconds with respect to the first
date in the time dimension array or the `new_start_date`. The parameter `new_start_date` can
be any string parseable by the [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/)
module or a `Dates.DateTime` object. The parameter `shift_by` is a function that takes in
Dates.DateTime elements and return Dates.DateTime elements. The start date is added to the
attributes of the `OutputVar`. The parameter `shift_by` is a function that takes in
`Dates.DateTime` elements and returns `Dates.DateTime` elements. This function is applied to
each element of the time array. Shifting the dates and converting to seconds is done in that
order.
"""
function OutputVar(
    path::String,
    short_name = nothing;
    new_start_date = nothing,
    shift_by = identity,
)
    var = read_var(path; short_name)
    # Check if it is possible to convert dates to seconds in the time dimension
    if (
        has_time(var) && (
            eltype(times(var)) <:
            Union{Dates.DateTime, NCDatasets.CFTime.AbstractCFDateTime}
        )
    )
        var = _dates_to_seconds(
            read_var(path; short_name),
            new_start_date = new_start_date,
            shift_by = shift_by,
        )
    end
    return var
end

"""
    read_var(path::String; short_name = nothing)

Read the `short_name` variable in the given NetCDF file.

When `short_name` is `nothing`, automatically identify the name of the variable. If multiple
variables are present, the last one in alphabetical order is chosen.

When `units` is among the attributes, try to parse it and convert it into an
[`Unitful`](https://painterqubits.github.io/Unitful.jl) object. `OutputVar`s with `Unitful`
support automatic unit conversions.

If you want to access `units` as a string, look at [`units`](@ref) function.

Example
=========

```julia
simdir = SimDir("my_output")
read_var(simdir.variable_paths["hu"]["inst"])

read_var("my_netcdf_file.nc", short_name = "ts")
```
"""
function read_var(path::String; short_name = nothing)
    NCDatasets.NCDataset(path) do nc
        # First, if short_name is nothing, we have to identify the name of the variable by
        # finding what is not a dimension
        unordered_dims = NCDatasets.dimnames(nc)
        isnothing(short_name) &&
            (short_name = pop!(setdiff(keys(nc), unordered_dims)))

        dims =
            map(NCDatasets.dimnames(nc[short_name])) do dim_name
                return dim_name => Array(nc[dim_name])
            end |> OrderedDict
        attribs = Dict(k => v for (k, v) in nc[short_name].attrib)

        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[short_name])

        return OutputVar(attribs, dims, dim_attribs, data)
    end
end

"""
    read_var(paths::Vector{String}; short_name = nothing)

Read the `short_name` variable in the vector of NetCDF files.

This function automatically aggregates the NetCDF files along the time dimension. If this
is not possible, then an error is thrown. The attributes of the OutputVar is the same as the
attributes of the first NetCDF file.

When `short_name` is `nothing`, automatically identify the name of the variable. If multiple
variables are present, the last one in alphabetical order is chosen.

When `units` is among the attributes, try to parse it and convert it into an
[`Unitful`](https://painterqubits.github.io/Unitful.jl) object. `OutputVar`s with `Unitful`
support automatic unit conversions.

If you want to access `units` as a string, look at [`units`](@ref) function.

Example
=========

```julia
simdir = SimDir("my_output")
read_var(simdir.variable_paths["hu"]["inst"])

read_var(["my_netcdf_file1.nc", "my_netcdf_file2.nc"], short_name = "ts")
```
"""
function read_var(paths::Vector{String}; short_name = nothing)
    # Assume that `paths` contain all the datasets we want to try to stitch together
    # Assume that `paths` is sorted (will rely on simdir to return it sorted)

    # If there is only one path, then use `read_var(path::String)`
    # We do this because SimDir always store file paths in a vector even if there is only one
    # file path
    length(paths) == 1 && return read_var(first(paths))

    # Helper function to find the name of the time dimension
    function find_time_name(dim_names)
        for dim_name in dim_names
            dim_name in TIME_NAMES && return dim_name
        end
        return nothing
    end
    ncs = (NCDatasets.NCDataset(path) for path in paths)
    time_names = map(nc -> find_time_name(NCDatasets.dimnames(nc)), ncs)

    # Check if each NCDataset got a time dimension
    if nothing in time_names
        error(
            "Time dimension does not exist in one of the NetCDF files. Provide the specific folder in SimDir",
        )
    end

    # Check if the name for the time dimension is the same across alll NetCDF files
    if length(unique(time_names)) != 1
        error(
            " Names of the time dimension are not the same across all NetCDF files. Provide the specific folder in SimDir",
        )
    end

    time_dim_name = first(time_names)

    # Check if the start dates are the same
    unordered_dims = (NCDatasets.dimnames(nc) for nc in ncs)
    short_names = (
        !isnothing(short_name) ? short_name :
        pop!(setdiff(keys(nc), unordered_dim)) for
        (unordered_dim, nc) in zip(unordered_dims, ncs)
    )
    start_dates = (
        get(nc[short_name].attrib, "start_date", "") for
        (nc, short_name) in zip(ncs, short_names)
    )
    if length(unique(start_dates)) != 1
        error(
            "Start dates are not the same across all the NetCDF files. Provide the specific folder in SimDir",
        )
    end

    # Check if it make sense to concat time dimension arrays
    time_dim_arrays = (Array(nc[time_dim_name]) for nc in ncs)
    time_dim_concat = vcat(time_dim_arrays...)
    if !issorted(time_dim_concat)
        error(
            "Time dimension is not in increasing order after aggregating datasets. Provide the specific folder in SimDir",
        )
    end

    NCDatasets.NCDataset(paths, aggdim = time_dim_name) do nc
        # First, if short_name is nothing, we have to identify the name of the variable by
        # finding what is not a dimension
        unordered_dims = NCDatasets.dimnames(nc)
        isnothing(short_name) &&
            (short_name = pop!(setdiff(keys(nc), unordered_dims)))

        dims =
            map(NCDatasets.dimnames(nc[short_name])) do dim_name
                return dim_name => Array(nc[dim_name])
            end |> OrderedDict
        attribs = Dict(k => v for (k, v) in nc[short_name].attrib)

        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[short_name])

        return OutputVar(attribs, dims, dim_attribs, data)
    end
end

"""
    short_name(var::OutputVar)

Return the `short_name` of the given `var`, if available.

If not available, return an empty string.
"""
function short_name(var::OutputVar)
    get(var.attributes, "short_name", "")
end

"""
    long_name(var::OutputVar)

Return the `long_name` of the given `var`, if available.

If not available, return an empty string.
"""
function long_name(var::OutputVar)
    get(var.attributes, "long_name", "")
end

"""
    units(var::OutputVar)

Return the `units` of the given `var`, if available.

If not available, return an empty string.
"""
function units(var::OutputVar)
    string(get(var.attributes, "units", ""))
end

"""
    has_units(var::OutputVar)

Return whether the given `var` has `units` or not.
"""
function has_units(var::OutputVar)
    return haskey(var.attributes, "units")
end

# Implemented in ClimaAnalysisUnitfulExt
function _maybe_convert_to_unitful end

"""
    Var.convert_units(var, new_units; conversion_function = nothing)

Return a new `OutputVar` with converted physical units of `var` to `new_units`, if possible.

Automatic conversion happens when the units for `var` and `new_units` are both parseable.
When `var` does not have units (see also [`Var.has_units`](@ref)) or has no parseable units,
a conversion function `conversion_function` is required.

`conversion_function` has to be a function that takes one data point and returns the
transformed value.

Being parseable means that `Unitful.uparse` can parse the expression. Please, refer to the
documentation for [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) for more
information.

Examples
=======

Let us set up a trivial 1D `OutputVar` with units of meters per second and automatically
convert it to centimeters per second.

```jldoctest example1
julia> values = 0:100.0 |> collect;

julia> data = copy(values);

julia> attribs = Dict("long_name" => "speed", "units" => "m/s");

julia> dim_attribs = Dict{String, Any}();

julia> var = ClimaAnalysis.OutputVar(attribs, Dict("distance" => values), dim_attribs, data);

julia> ClimaAnalysis.has_units(var)
true

julia> var_cms = ClimaAnalysis.convert_units(var, "cm/s");

julia> extrema(var_cms.data)
(0.0, 10000.0)
```

Not all the units can be properly parsed, for example, assuming `bob=5lol`

```jldoctest example1
julia> attribs = Dict("long_name" => "speed", "units" => "bob/s");

julia> var_bob = ClimaAnalysis.OutputVar(attribs, Dict("distance" => values), dim_attribs, data);

julia> var_lols = ClimaAnalysis.convert_units(var, "lol/s", conversion_function = (x) -> 5x);

julia> extrema(var_lols.data)
(0.0, 500.0)
```

Failure to specify the `conversion_function` will produce an error.
"""
function convert_units end

"""
    set_units(var::OutputVar, units::AbstractString)

Set `units` for data in `var`.

!!! warning "Override existing units"
    If units already exist, this will override the units for data in `var`. To convert
    units, see [`Var.convert_units`](@ref)
"""
function set_units(var::OutputVar, units::AbstractString)
    converted_var = convert_units(var, units, conversion_function = identity)
    return converted_var
end

"""
    set_dim_units!(var::OutputVar, dim_name::AbstractString, units::AbstractString)

Set `units` for the `dim_name` dimension in `var`.

!!! warning "Override existing units"
    If units already exist for the dimension, this will override the units for the dimension
    in `var`.
"""
function set_dim_units!(
    var::OutputVar,
    dim_name::AbstractString,
    units::AbstractString,
)
    !haskey(var.dims, dim_name) &&
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    if haskey(var.dim_attributes, dim_name)
        push!(var.dim_attributes[dim_name], "units" => units)
    else
        var.dim_attributes[dim_name] = Dict("units" => units)
    end
    return nothing
end

"""
    is_z_1D(var::OutputVar)

Return whether the given `var`iable has an altitude dimension that is 1D.

When topography is present, the altitude dimension in the output variable is typically
multidimensional. The dimensions are (X, Y, Z), where (X, Y) are the horizontal dimensions.
In this case, `dims["z"]` is essentially a map that identifies the physical altitude of the
given point.
"""
function is_z_1D(var::OutputVar)
    has_altitude(var) || error("Variable does not have an altitude dimension")

    return length(size(altitudes(var))) == 1
end

"""
    isempty(var::OutputVar)

Determine whether an OutputVar is empty.
"""
function Base.isempty(var::OutputVar)
    # Do not include :interpolant because var.interpolant is Nothing if data is
    # zero dimensional and empty and isempty(Nothing) throws an error
    return map(
        fieldname -> isempty(getproperty(var, fieldname)),
        filter(x -> x != :interpolant, fieldnames(OutputVar)),
    ) |> all
end

function Base.copy(var::OutputVar)
    fields = fieldnames(OutputVar)
    # We have nested containers and we have to make sure we hand ownership off,
    # so we deepcopy
    return OutputVar([deepcopy(getfield(var, field)) for field in fields]...)
end

"""
    _reduce_over(reduction::F, dim::String, var::OutputVar, args...; kwargs...)

Apply the given reduction over the given dimension.

`reduction` has to support the `dims` key. Additional arguments are passed to `reduction`.

The return type is an `OutputVar` with the same attributes, the new data, and the dimension
dropped.

Example
=========

Average over latitudes

```julia
import Statistics: mean
long = 0.:180. |> collect
lat = 0.:90. |> collect
data = reshape(1.:91*181., (181, 91))
dims = Dict(["lat" => lat, "long" => long])
var = OutputVar(dims, data)
_reduce_over(mean, "lat", var)
```
"""
function _reduce_over(
    reduction::F,
    dim::String,
    var::OutputVar,
    args...;
    kwargs...,
) where {F <: Function}
    dim_index = var.dim2index[dim]

    # squeeze removes the unnecessary singleton dimension
    data = squeeze(
        reduction(var.data, args..., dims = dim_index, kwargs...),
        dims = (dim_index,),
    )

    # If we reduce over a dimension, we have to remove it
    dims = copy(var.dims)
    dim_attributes = copy(var.dim_attributes)
    pop!(dims, dim)
    haskey(var.dim_attributes, dim) && pop!(dim_attributes, dim)
    return OutputVar(copy(var.attributes), dims, dim_attributes, copy(data))
end

"""
    average_lat(var::OutputVar; ignore_nan = true, weighted = false)

Return a new OutputVar where the values on the latitudes are averaged arithmetically.

When `weighted` is `true`, weight the average over `cos(lat)`.
"""
function average_lat(var; ignore_nan = true, weighted = false)
    if weighted
        var = copy(var)
        lats = latitudes(var)
        abs(maximum(lats)) >= 0.5π ||
            @warn "Detected latitudes are small. If units are radians, results will be wrong"

        weights_1d = cosd.(lats)
        lat_index = var.dim2index[latitude_name(var)]
        weights = ones(size(var.data))
        # Create a bitmask for the NaN's, we use this to remove weights in the normalization (with nanmean)
        nan_mask = ifelse.(isnan.(var.data), NaN, 1)
        for index in CartesianIndices(weights)
            index_tuple =
                ntuple(d -> d == lat_index ? Colon() : index[d], ndims(weights))

            weights[index_tuple...] .= weights_1d
            weights[index_tuple...] ./=
                nanmean(nan_mask[index_tuple...] .* weights_1d)
        end
        var.data .*= weights
    end

    reduced_var =
        _reduce_over(ignore_nan ? nanmean : mean, latitude_name(var), var)
    weighted &&
        haskey(var.attributes, "long_name") &&
        (reduced_var.attributes["long_name"] *= " weighted")
    _update_long_name_generic!(reduced_var, var, latitude_name(var), "averaged")


    return reduced_var
end

"""
    weighted_average_lat(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values on the latitudes are averaged arithmetically
with weights of `cos(lat)`.
"""
weighted_average_lat(var; ignore_nan = true) =
    average_lat(var; ignore_nan, weighted = true)

"""

    average_lon(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values on the longitudes are averaged arithmetically.
"""
function average_lon(var; ignore_nan = true)
    reduced_var =
        _reduce_over(ignore_nan ? nanmean : mean, longitude_name(var), var)
    _update_long_name_generic!(
        reduced_var,
        var,
        longitude_name(var),
        "averaged",
    )
    return reduced_var
end

"""
    average_x(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `x` dimension are averaged arithmetically.
"""
function average_x(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "x", var)
    _update_long_name_generic!(reduced_var, var, "x", "averaged")
    return reduced_var
end

"""
    average_y(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `y` dimension are averaged arithmetically.
"""
function average_y(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "y", var)
    _update_long_name_generic!(reduced_var, var, "y", "averaged")
    return reduced_var
end

"""
    average_xy(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along both horizontal dimensions `x` and `y`
are averaged arithmetically.
"""
function average_xy(var; ignore_nan = true)
    fn = ignore_nan ? nanmean : mean
    reduced_var = _reduce_over(fn, "x", _reduce_over(fn, "y", var))

    first_x, last_x = range_dim(var, "x")
    first_y, last_y = range_dim(var, "y")
    units_x = dim_units(var, "x")
    units_y = dim_units(var, "y")
    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged horizontally over x ($first_x to $last_x$units_x) and y ($first_y to $last_y$units_y)"
    end

    return reduced_var
end

"""
    average_time(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values are averaged arithmetically in time.
"""
function average_time(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, time_name(var), var)
    _update_long_name_generic!(reduced_var, var, time_name(var), "averaged")
    return reduced_var
end

"""
    dim_units(var::OutputVar, dim_name)

Return the `units` of the given `dim_name` in `var`, if available.

If not available, return an empty string.
"""
function dim_units(var::OutputVar, dim_name)
    !haskey(var.dims, dim_name) &&
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    # Double get because var.dim_attributes is a dictionry whose values are dictionaries
    get(get(var.dim_attributes, dim_name, Dict()), "units", "")
end

"""
    range_dim(var::OutputVar, dim_name)

Return the range of the dimension `dim_name` in `var`.

Range here is a tuple with the minimum and maximum of `dim_name`.
"""
function range_dim(var::OutputVar, dim_name)
    !haskey(var.dims, dim_name) &&
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    first_elt = first(var.dims[dim_name])
    last_elt = last(var.dims[dim_name])
    return first_elt, last_elt
end

"""
    _update_long_name_generic!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_name,
    operation_name,
)

Used by reductions (e.g., average) to update the long name of `reduced_var` by describing
the operation being used to reduce the data and the associated units.
"""
function _update_long_name_generic!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_name,
    operation_name,
)
    dim_of_units = dim_units(var, dim_name)
    first_elt, last_elt = range_dim(var, dim_name)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " $operation_name over $dim_name ($first_elt to $last_elt$dim_of_units)"
    end
    return nothing
end

"""
    center_longitude!(var::OutputVar, lon::Real)

Shift the longitudes in `var` so that `lon` is the center one.

This is useful to center the global projection to the 180 meridian instead of the 0.
"""
function center_longitude!(var, lon)
    lon_name = longitude_name(var)

    old_center_lon_index = nearest_index(var.dims[lon_name], lon)
    half_index = div(length(var.dims[lon_name]), 2)

    # half_index = old_index + shift => shift = half_index - old_index
    shift = half_index - old_center_lon_index

    # We do not use circshift! because it can lead to unpredictable problems when mutating
    shifted_lon = circshift(var.dims[lon_name], shift)
    var.dims[lon_name] = shifted_lon

    lon_dim_index = var.dim2index[lon_name]

    # Prepare the shift tuple for the data array: do not shift, except for the dimension
    # corresponding to the longitude
    shift_tuple = zeros(length(var.dims))
    shift_tuple[lon_dim_index] = shift

    shifted_data = circshift(var.data, shift)
    var.data .= shifted_data
end

"""
    _slice_general(var::OutputVar, val, dim_name)

Return a new OutputVar by selecting the available index closest to the given `val` for the
given dimension
"""
function _slice_general(var, val, dim_name)
    haskey(var.dims, dim_name) ||
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")

    nearest_index_val = nearest_index(var.dims[dim_name], val)
    _slice_over(data; dims) = selectdim(data, dims, nearest_index_val)
    reduced_var = _reduce_over(_slice_over, dim_name, var)

    # Let's try adding this operation to the long_name, if possible (ie, if the correct
    # attributes are available)
    if haskey(var.attributes, "long_name") &&
       haskey(var.dim_attributes, dim_name) &&
       haskey(var.dim_attributes[dim_name], "units")
        dim_array = var.dims[dim_name]
        dim_units = var.dim_attributes[dim_name]["units"]
        cut_point = dim_array[nearest_index_val]
        if (dim_name == "time" || dim_name == "t") && dim_units == "s"
            # Dimension is time and units are seconds. Let's convert them to something nicer
            pretty_timestr = seconds_to_prettystr(cut_point)
            reduced_var.attributes["long_name"] *= " $dim_name = $pretty_timestr"
        else
            reduced_var.attributes["long_name"] *= " $dim_name = $cut_point $dim_units"
        end
        reduced_var.attributes["slice_$dim_name"] = "$cut_point"
        reduced_var.attributes["slice_$(dim_name)_units"] = dim_units
    end
    return reduced_var
end

"""
    slice(var::OutputVar, kwargs...)

Return a new OutputVar by slicing across dimensions as defined by the keyword arguments.

Example
===========
```julia
slice(var, lat = 30, lon = 20, time = 100)
```
"""
function slice(var; kwargs...)
    sliced_var = var
    for (dim_name, val) in kwargs
        sliced_var = _slice_general(sliced_var, val, String(dim_name))
    end
    return sliced_var
end

"""
    window(var::OutputVar, dim_name; left = nothing, right = nothing)

Return a new OutputVar by selecting the values of the given `dim`ension that are between
`left` and `right`.

If `left` and/or `right` are `nothing`, assume beginning (or end) of the array.

Example
===========
```julia
window(var, 'lat', left = -50, right = 50)
```
"""
function window(var, dim_name; left = nothing, right = nothing)
    haskey(var.dims, dim_name) ||
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")

    nearest_index_left =
        isnothing(left) ? 1 : nearest_index(var.dims[dim_name], left)
    nearest_index_right =
        isnothing(right) ? length(var.dims[dim_name]) :
        nearest_index(var.dims[dim_name], right)

    (nearest_index_right >= nearest_index_left) ||
        error("Right window value has to be larger than left one")

    # Take only what's between nearest_index_left and nearest_index_right
    reduced_data = selectdim(
        var.data,
        var.dim2index[dim_name],
        nearest_index_left:nearest_index_right,
    )

    dims = copy(var.dims)
    reduced_dim = var.dims[dim_name][nearest_index_left:nearest_index_right]
    dims[dim_name] = reduced_dim

    dim_attributes = copy(var.dim_attributes)
    return OutputVar(copy(var.attributes), dims, dim_attributes, reduced_data)
end

"""
    (x::OutputVar)(target_coord)

Interpolate variable `x` onto the given `target_coord` coordinate using
multilinear interpolation.

Extrapolation is now allowed and will throw a `BoundsError` in most cases.

If any element of the arrays of the dimensions is a Dates.DateTime, then interpolation is
not possible. Interpolations.jl do not support making interpolations for dates. If the
longitudes span the entire range and are equispaced, then a periodic boundary condition is
added for the longitude dimension. If the latitudes span the entire range and are
equispaced, then a flat boundary condition is added for the latitude dimension. In all other
cases, an error is thrown when extrapolating outside of the array of the dimension.

Example
=======
```jldoctest
julia> import ClimaAnalysis;

julia> time = 100.0:110.0 |> collect;

julia> z = 0.0:20.0 |> collect;

julia> data = reshape(1.0:(11 * 21), (11, 21));

julia> var2d = ClimaAnalysis.OutputVar(Dict("time" => time, "z" => z), data); var2d.([[105., 10.], [105.5, 10.5]])
2-element Vector{Float64}:
 116.0
 122.0
```
"""
function (x::OutputVar)(target_coord)
    isnothing(x.interpolant) && error(
        "Splines cannot be constructed because one (or more) of the dimensions of variable is not 1D",
    )
    return x.interpolant(target_coord...)
end

"""
    arecompatible(x::OutputVar, y::OutputVar)

Return whether two `OutputVar` are defined on the same physical space

This is accomplished by comparing `dims` and `dim_attributes` (the latter because they might contain information about the units).

We assume that:
- `dim2index` and `index2dim` where correctly created and they reflect `dims`
- `data` is also consistent with `dims`,

We also *do not* check units for `data`.
"""
function arecompatible(x::OutputVar, y::OutputVar)
    x_dims = collect(keys(x.dims))
    y_dims = collect(keys(y.dims))
    x_units = (dim_units(x, dim_name) for dim_name in x_dims)
    y_units = (dim_units(y, dim_name) for dim_name in y_dims)

    for (x_dim, x_unit, y_dim, y_unit) in zip(x_dims, x_units, y_dims, y_units)
        x_unit == "" && @warn "Missing units for dimension $x_dim in x"
        y_unit == "" && @warn "Missing units for dimension $y_dim in y"
        x_unit != y_unit && return false
    end
    return x.dims == y.dims
end

"""
    _check_dims_consistent(x::OutputVar, y::OutputVar)

Check if the number, name, and unit of dimensions in `x` and `y` are consistent.

If the unit for a dimension is missing, then the unit is not consistent for that dimension.
"""
function _check_dims_consistent(x::OutputVar, y::OutputVar)
    # Check if the number of dimensions is the same
    x_num_dims = length(x.dims)
    y_num_dims = length(y.dims)
    x_num_dims != y_num_dims && error(
        "Number of dimensions do not match between x ($x_num_dims) and y ($y_num_dims)",
    )

    # Check if the dimensions agree with each other
    conventional_dim_name_x = conventional_dim_name.(keys(x.dims))
    conventional_dim_name_y = conventional_dim_name.(keys(y.dims))
    mismatch_conventional_dim_name =
        conventional_dim_name_x .!= conventional_dim_name_y
    any(mismatch_conventional_dim_name) && error(
        "Dimensions do not agree between x ($conventional_dim_name_x) and y ($conventional_dim_name_y)",
    )

    x_dims = collect(keys(x.dims))
    y_dims = collect(keys(y.dims))
    x_units = [dim_units(x, dim_name) for dim_name in x_dims]
    y_units = [dim_units(y, dim_name) for dim_name in y_dims]

    # Check for any missing units (missing units are represented with an empty string)
    missing_x = (x_units .== "")
    missing_y = (y_units .== "")
    (any(missing_x) && any(missing_y)) && error(
        "Units for dimensions $(x_dims[missing_x]) are missing in x and units for dimensions $(y_dims[missing_y]) are missing in y",
    )
    any(missing_x) &&
        error("Units for dimensions $(x_dims[missing_x]) are missing in x")
    any(missing_y) &&
        error("Units for dimensions $(x_dims[missing_y]) are missing in y")

    # Check if units match between dimensions
    not_consistent_units = (x_units .!= y_units)
    any(not_consistent_units) && error(
        "Units for dimensions $(x_dims[not_consistent_units]) in x is not consistent with units for dimensions $(y_dims[not_consistent_units]) in y",
    )
    return nothing
end

"""
    reordered_as(src_var::OutputVar, dest_var::OutputVar)

Reorder the dimensions in `src_var` to match the ordering of dimensions in `dest_var`.
"""
function reordered_as(src_var::OutputVar, dest_var::OutputVar)
    # Get the conventional dim names for both src_var and dest_var
    conventional_dim_name_src = conventional_dim_name.(keys(src_var.dims))
    conventional_dim_name_dest = conventional_dim_name.(keys(dest_var.dims))

    # Check if the dimensions are the same (order does not matter)
    Set(conventional_dim_name_src) == Set(conventional_dim_name_dest) || error(
        "Dimensions are not the same between src ($conventional_dim_name_src) and dest ($conventional_dim_name_dest)",
    )

    # Find permutation indices to reorder dims
    reorder_indices =
        indexin(conventional_dim_name_dest, conventional_dim_name_src)

    # Reorder dims, dim_attribs, and data, but not attribs
    ret_dims = deepcopy(src_var.dims)
    ret_dims = OrderedDict(collect(ret_dims)[reorder_indices])
    ret_attribs = deepcopy(src_var.attributes)

    # Cannot assume that every dimension is present in dim_attribs so we loop to reorder the
    # best we can and merge with src_var.dim_attributes to add any remaining pairs to
    # ret_dim_attribs
    ret_dim_attribs = empty(src_var.dim_attributes)
    src_var_dim_attribs = src_var.dim_attributes |> deepcopy
    src_var_dim_names = collect(keys(src_var.dims))
    for idx in reorder_indices
        dim_name = src_var_dim_names[idx]
        haskey(src_var_dim_attribs, dim_name) &&
            (ret_dim_attribs[dim_name] = src_var_dim_attribs[dim_name])
    end
    merge!(ret_dim_attribs, src_var_dim_attribs)

    ret_data = copy(src_var.data)
    ret_data = permutedims(ret_data, reorder_indices)
    return OutputVar(ret_attribs, ret_dims, ret_dim_attribs, ret_data)
end

"""
    resampled_as(src_var::OutputVar, dest_var::OutputVar)

Resample `data` in `src_var` to `dims` in `dest_var`.

The resampling performed here is a 1st-order linear resampling.
"""
function resampled_as(src_var::OutputVar, dest_var::OutputVar)
    src_var = reordered_as(src_var, dest_var)
    _check_dims_consistent(src_var, dest_var)

    src_resampled_data =
        [src_var(pt) for pt in Base.product(values(dest_var.dims)...)]

    # Construct new OutputVar to return
    src_var_ret_dims = empty(src_var.dims)

    # Loop because names could be different in src_var compared to dest_var
    # (e.g., `long` in one and `lon` in the other)
    for (dim_name, dim_data) in zip(keys(src_var.dims), values(dest_var.dims))
        src_var_ret_dims[dim_name] = copy(dim_data)
    end
    scr_var_ret_attribs = deepcopy(src_var.attributes)
    scr_var_ret_dim_attribs = deepcopy(src_var.dim_attributes)
    return OutputVar(
        scr_var_ret_attribs,
        src_var_ret_dims,
        scr_var_ret_dim_attribs,
        src_resampled_data,
    )
end

"""
    integrate_lonlat(var::OutputVar)

Integrate `data` in `var` on longitude and latitude with a first-order scheme. `data` has to
be discretized on longitude and latitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints for integrating
longitude and latitude. The units for longitude and latitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lonlat(var::OutputVar)
    var_integrate_lon = var |> integrate_lon
    # Update long name so that we get "...integrated over lat... and integrated over lon..."
    # instead of "...integrated over lat... integrated over lon..."
    if haskey(var_integrate_lon.attributes, "long_name")
        var_integrate_lon.attributes["long_name"] *= " and"
    end
    return var_integrate_lon |> integrate_lat
end

"""
    integrate_lon(var::OutputVar)

Integrate `data` in `var` on longitude with a first-order scheme. `data` has to be
discretized on longitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
longitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lon(var::OutputVar)
    has_longitude(var) || error("var does not has longitude as a dimension")
    lon_name = longitude_name(var)
    return _integrate_over_angle(var, Numerics._integrate_lon, lon_name)
end

"""
    integrate_lat(var::OutputVar)

Integrate `data` in `var` on latitude with a first-order scheme. `data` has to be
discretized on latitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
latitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lat(var::OutputVar)
    has_latitude(var) || error("var does not has latitude as a dimension")
    lat_name = latitude_name(var)
    return _integrate_over_angle(var, Numerics._integrate_lat, lat_name)
end

"""
    _integrate_over_angle(var::OutputVar, integrate_on, angle_dim_name)

Integrate `data` in `var` on `angle_dim_name` with a first-order scheme. `data` has to be
discretized on `angle_dim_name`.

`angle_dim_name` is the name of the angle that is being integrated over. `integrate_on` is a
function that compute the integration for data over `angle_dim_name`.
"""
function _integrate_over_angle(var::OutputVar, integrate_on, angle_dim_name)
    # Enforce constraint that unit is degree because we compute integration weights assuming
    # degrees (see functions _integration_weights_lon_left and
    # _integration_weights_lon_equispaced for examples in Numerics.jl)
    deg_unit_names = [
        "degrees",
        "degree",
        "deg",
        "degs",
        "°",
        "degrees_north",
        "degrees_east",
    ]
    angle_dim_unit = dim_units(var, angle_dim_name)
    lowercase(angle_dim_unit) in deg_unit_names ||
        error("The unit for $angle_dim_name is missing or is not degree")

    integrated_var = _reduce_over(
        integrate_on,
        angle_dim_name,
        var,
        var.dims[angle_dim_name],
    )

    _update_long_name_generic!(
        integrated_var,
        var,
        angle_dim_name,
        "integrated",
    )
    return integrated_var
end


"""
    split_by_season(var::OutputVar)

Return a vector of four `OutputVar`s split by season.

The months of the seasons are March to May, June to August, September to November, and
December to February. The order of the vector is MAM, JJA, SON, and DJF. If there are no
dates found for a season, then the `OutputVar` for that season will be an empty `OutputVar`.

The function will use the start date in `var.attributes["start_date"]`. The unit of time is
expected to be second. Also, the interpolations will be inaccurate in time intervals
outside of their respective season for the returned `OutputVar`s.
"""
function split_by_season(var::OutputVar)
    # Check time exists and unit is second
    has_time(var) || error("Time is not a dimension in var")
    dim_units(var, time_name(var)) == "s" ||
        error("Unit for time is not second")

    # Check start date exists
    haskey(var.attributes, "start_date") ?
    start_date = Dates.DateTime(var.attributes["start_date"]) :
    error("Start date is not found in var")

    season_dates = split_by_season(time_to_date.(start_date, times(var)))
    season_times =
        (date_to_time.(start_date, season) for season in season_dates)

    # Split data according to seasons
    season_data = (
        collect(
            _data_at_dim_vals(
                var.data,
                times(var),
                var.dim2index[time_name(var)],
                season_time,
            ),
        ) for season_time in season_times
    )

    # Construct an OutputVar for each season
    return map(season_times, season_data) do time, data
        if isempty(time)
            dims = empty(var.dims)
            data = similar(var.data, 0)
            return OutputVar(dims, data)
        end
        ret_dims = deepcopy(var.dims)
        ret_attribs = deepcopy(var.attributes)
        ret_dim_attribs = deepcopy(var.dim_attributes)
        ret_dims[time_name(var)] = time
        OutputVar(ret_attribs, ret_dims, ret_dim_attribs, data)
    end
end

"""
    _check_sim_obs_units_consistent(sim::OutputVar, obs::OutputVar)

Check if the number of dimensions are two, the `data` in `sim` and `obs` is missing units or
not, and if the units of data are the same in `sim` and `obs`.

This function does not check if the dimensions are longitude and latitude in `sim` and `obs`
because `integrate_lonlat` (in `bias` and `squared_error`) handles that. The function also
does not check if the units of dimensions in `sim` and `obs` are the same because
`resampled_as` (in `bias` and `squared_error`) handles that.
"""
function _check_sim_obs_units_consistent(sim::OutputVar, obs::OutputVar)
    # Check number of dimensions
    sim_num_dims = length(sim.dims)
    obs_num_dims = length(obs.dims)
    ((sim_num_dims != 2) || (obs_num_dims != 2)) && error(
        "There are not only two dimensions in sim ($sim_num_dims) or obs ($obs_num_dims).",
    )

    # Check units for data is not missing
    sim_data_units = units(sim)
    obs_data_units = units(obs)
    sim_data_units == "" && error("Unit is missing in data for sim")
    obs_data_units == "" && error("Unit is missing in data for obs")

    # Check if units of data match between sim and obs
    sim_data_units == obs_data_units || error(
        "Units do not match between the data in sim ($sim_data_units) and obs ($obs_data_units)",
    )
    return nothing
end

"""
    bias(sim::OutputVar, obs::OutputVar; mask = nothing)

Return a `OutputVar` whose data is the bias (`sim.data - obs.data`) and compute the global
bias of `data` in `sim` and `obs` over longitude and latitude. The result is stored in
`var.attributes["global_bias"]`.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`. Attributes in `sim` and `obs` will be thrown away. The long name and short
name of the returned `OutputVar` will be updated to reflect that a bias is computed.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`global_bias`](@ref), [`squared_error`](@ref), [`global_mse`](@ref),
[`global_rmse`](@ref).
"""
function bias(sim::OutputVar, obs::OutputVar; mask = nothing)
    _check_sim_obs_units_consistent(sim, obs)

    # Resample obs on sim to ensure the size of data in sim and obs are the same and the
    # dims are the same
    obs_resampled = resampled_as(obs, sim)

    # Compute bias
    bias = sim - obs_resampled

    # Do this because we do not want to store global bias as a string and unit could be Unitful
    ret_attributes = Dict{keytype(bias.attributes), Any}(bias.attributes)

    # Add units back for bias
    ret_attributes["units"] = units(sim)

    # Add short and long name
    ret_attributes["short_name"] = "sim-obs"
    ret_attributes["long_name"] = "SIM - OBS"
    if !isempty(short_name(sim))
        ret_attributes["short_name"] *= "_" * short_name(sim)
        ret_attributes["long_name"] *= " " * short_name(sim)
    end

    # Compute global bias and store it as an attribute
    !isnothing(mask) && (bias = mask(bias))
    integrated_bias = integrate_lonlat(bias).data
    ones_var = OutputVar(
        bias.attributes,
        bias.dims,
        bias.dim_attributes,
        ones(size(bias.data)),
    )
    !isnothing(mask) && (ones_var = mask(ones_var))
    normalization = integrate_lonlat(ones_var).data
    # Do ./ instead of / because we are dividing between zero dimensional arrays
    global_bias = integrated_bias ./ normalization
    ret_attributes["global_bias"] = global_bias
    return OutputVar(ret_attributes, bias.dims, bias.dim_attributes, bias.data)
end

"""
    global_bias(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global bias of `data` in `sim` and `obs` over longitude and latitude.

This function is currently only implemented for `OutputVar`s with only the dimensions
longitude and latitude. Units must be supplied for data and dimensions in `sim` and `obs`.
The units for longitude and latitude should be degrees. Resampling is done automatically by
resampling `obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`bias`](@ref), [`squared_error`](@ref), [`global_mse`](@ref),
[`global_rmse`](@ref).
"""
function global_bias(sim::OutputVar, obs::OutputVar; mask = nothing)
    bias_var = bias(sim, obs, mask = mask)
    return bias_var.attributes["global_bias"]
end

"""
    squared_error(sim::OutputVar, obs::OutputVar; mask = nothing)

Return a `OutputVar` whose data is the squared error (`(sim.data - obs.data)^2`) and compute
the global mean squared error (MSE) and global root mean squared error (RMSE) of `data` in
`sim` and `obs` over longitude and latitude. The result is stored in `var.attributes["mse"]`
and `var.attributes["rmse"]`.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`. Attributes in `sim` and `obs` will be thrown away. The long name and short
name of the returned `OutputVar` will be updated to reflect that a squared error is computed.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`global_mse`](@ref), [`global_rmse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function squared_error(sim::OutputVar, obs::OutputVar; mask = nothing)
    _check_sim_obs_units_consistent(sim, obs)

    # Resample obs on sim to ensure the size of data in sim and obs are the same and the
    # dims are the same
    obs_resampled = resampled_as(obs, sim)

    # Compute squared error
    # Do not use ^ since ^ is not defined between a OutputVar and Real
    squared_error = (sim - obs_resampled) * (sim - obs_resampled)

    # Do this because we do not want to store global mse and rmse as strings
    ret_attributes = Dict{String, Any}(squared_error.attributes)

    # Add units back for bias
    # Always treat as a string since the string representation of something type Unitful is
    # not always parseable as a Unitful.Unit (see:
    # https://github.com/PainterQubits/Unitful.jl/issues/412)
    ret_attributes["units"] = "($(units(sim)))^2"

    # Add short and long name
    ret_attributes["short_name"] = "(sim-obs)^2"
    ret_attributes["long_name"] = "(SIM - OBS)^2"
    if !isempty(short_name(sim))
        ret_attributes["short_name"] *= "_" * short_name(sim)
        ret_attributes["long_name"] *= " " * short_name(sim)
    end

    # Compute global mse and global rmse and store it as an attribute
    !isnothing(mask) && (squared_error = mask(squared_error))
    integrated_squared_error = integrate_lonlat(squared_error).data
    ones_var = OutputVar(
        squared_error.attributes,
        squared_error.dims,
        squared_error.dim_attributes,
        ones(size(squared_error.data)),
    )
    !isnothing(mask) && (ones_var = mask(ones_var))
    normalization = integrate_lonlat(ones_var).data
    # Do ./ instead of / because we are dividing between zero dimensional arrays
    mse = integrated_squared_error ./ normalization
    ret_attributes["global_mse"] = mse
    ret_attributes["global_rmse"] = sqrt(mse)
    return OutputVar(
        ret_attributes,
        squared_error.dims,
        squared_error.dim_attributes,
        squared_error.data,
    )
end

"""
    global_mse(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global mean squared error (MSE) of `data` in `sim` and `obs` over longitude and
latitude.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`squared_error`](@ref), [`global_rmse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function global_mse(sim::OutputVar, obs::OutputVar; mask = nothing)
    squared_error_var = squared_error(sim, obs, mask = mask)
    return squared_error_var.attributes["global_mse"]
end

"""
    global_rmse(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global root mean squared error (RMSE) of `data` in `sim` and `obs` over longitude
and latitude.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`squared_error`](@ref), [`global_mse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function global_rmse(sim::OutputVar, obs::OutputVar; mask = nothing)
    squared_error_var = squared_error(sim, obs, mask = mask)
    return squared_error_var.attributes["global_rmse"]
end

"""
    _dates_to_seconds(var::OutputVar;
                      new_start_date = nothing,
                      shift_by = identity)

Convert dates in time dimension to seconds with respect to the first date in the time
dimension array or the `new_start_date`.

Dates in the time dimension are automatically converted to seconds with respect to the first
date in the time dimension array or the `new_start_date`. The parameter `new_start_date` can
be any string parseable by the [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/)
module or a `Dates.DateTime` object. The parameter `shift_by` is a function that takes in
Dates.DateTime elements and return Dates.DateTime elements. The start date is added to the
attributes of the `OutputVar`. The parameter `shift_by` is a function that takes in
`Dates.DateTime` elements and returns `Dates.DateTime` elements. This function is applied to
each element of the time array. Shifting the dates and converting to seconds is done in that
order.

Note that this function only works for the time dimension and will not work for the date
dimension.
"""
function _dates_to_seconds(
    var::OutputVar;
    new_start_date = nothing,
    shift_by = identity,
)
    has_time(var) || error(
        "Converting from dates to seconds is only supported for the time dimension",
    )
    eltype(times(var)) <:
    Union{Dates.DateTime, NCDatasets.CFTime.AbstractCFDateTime} ||
        error("Type of time dimension is not dates")

    # Reinterpret everything as DateTime.Dates
    time_arr = copy(times(var))
    if eltype(time_arr) <: NCDatasets.CFTime.AbstractCFDateTime
        time_arr = map(t -> reinterpret(Dates.DateTime, t), time_arr)
    end

    # Preprocess time_arr by shifting dates
    if !isnothing(shift_by)
        time_arr .= shift_by.(time_arr)
    end

    # Convert from dates to seconds using the first date in the time dimension array as the
    # start date or the new_start_date
    start_date = isnothing(new_start_date) ? time_arr[begin] : new_start_date

    # Handle the case if start_date is a DateTime or string; if it is the latter, then try
    # to parse it as a DateTime
    start_date isa AbstractString && (start_date = Dates.DateTime(start_date))
    time_arr = map(date -> date_to_time(start_date, date), time_arr)

    # Remake OutputVar
    ret_attribs = deepcopy(var.attributes)
    ret_attribs["start_date"] = string(start_date) # add start_date as an attribute
    ret_dim_attribs = deepcopy(var.dim_attributes)
    ret_dim_attribs[time_name(var)]["units"] = "s" # add unit
    var_dims = deepcopy(var.dims)
    ret_dims_generator = (
        conventional_dim_name(dim_name) == "time" ? dim_name => time_arr :
        dim_name => dim_data for (dim_name, dim_data) in var_dims
    )
    ret_dims = OrderedDict(ret_dims_generator...)
    ret_data = copy(var.data)
    return OutputVar(ret_attribs, ret_dims, ret_dim_attribs, ret_data)
end

"""
    shift_to_start_of_previous_month(var::OutputVar)

Shift the times in the time dimension to the start of the previous month.

After applying this function, the start date in the attributes correspond to the first
element in the time array.

This function is helpful in ensuring consistency in dates between simulation and
observational data. One example of this is when adjusting monthly averaged data. For
instance, suppose that data on 2010-02-01 in the `OutputVar` corresponds to the monthly
average for January. This function shifts the times so that 2010-01-01 will correspond to
the monthly average for January.

Note that this function only works for the time dimension and will not work for the date
dimension.
"""
function shift_to_start_of_previous_month(var::OutputVar)
    # Check if time dimension exists, floats are in the array, and unit of data is
    # second
    has_time(var) || error("Time is not a dimension of var")
    eltype(times(var)) <: Dates.DateTime && ("Dates found in time array")
    dim_units(var, time_name(var)) != "s" &&
        error("Unit of data is not in second")

    # Convert to seconds to dates
    date_arr =
        Dates.Second.(times(var)) .+
        Dates.DateTime.(var.attributes["start_date"])

    # Apply transformations (find first day of month and subtract one month)
    date_arr .=
        date_arr .|> Dates.firstdayofmonth .|> date -> date - Dates.Month(1)

    # Check for duplicate dates
    unique_dates = unique(date_arr)
    unique_dates != date_arr && error(
        "Dates are not unique after applying shift_to_start_of_previous_month",
    )

    # Convert from dates to seconds
    start_date = date_arr[begin]
    time_arr = map(date -> date_to_time(start_date, date), date_arr)

    ret_attribs = deepcopy(var.attributes)
    ret_attribs["start_date"] = string(start_date)
    ret_dims = deepcopy(var.dims)
    ret_dims["time"] = time_arr
    ret_dim_attributes = deepcopy(var.dim_attributes)
    ret_data = copy(var.data)
    return OutputVar(ret_attribs, ret_dims, ret_dim_attributes, ret_data)
end

"""
    apply_landmask(var::OutputVar)

Apply a land mask to `var` by zeroing any data whose coordinates are located on land.
"""
function apply_landmask(var::OutputVar)
    _apply_lonlat_mask(var, LAND_MASK)
end

"""
    apply_oceanmask(var::OutputVar)

Apply an ocean mask to `var` by zeroing any data whose coordinates are in the ocean.
"""
function apply_oceanmask(var::OutputVar)
    _apply_lonlat_mask(var, OCEAN_MASK)
end

"""
    _apply_lonlat_mask(var, mask::AbstractString)

Apply a mask using the NCDataset found at `mask`.

The dimensions of the mask should only contain longitude and latitude and are in that order.
"""
function _apply_lonlat_mask(var, mask::AbstractString)
    mask_var = OutputVar(mask)
    apply_mask = make_lonlat_mask(mask_var)
    return apply_mask(var)
end

"""
    make_lonlat_mask(var;
                     set_to_val = nothing,
                     true_val = NaN,
                     false_val = 1.0)

Return a masking function that takes in a `OutputVar` and mask the data using `var.data`.

The parameter `set_to_val` is a function that takes in an element of `var.data` and return
`True` or `False`. If `set_to_nan` returns `True`, then the value will be `true_val` in the
mask and if `set_to_nan` returns `False`, then the value will be `false_val` in the mask.

If `set_to_nan` is `nothing`, then no transformation is done and `var.data` is used. This is
helpful if `var.data` is already an array of NaNs and ones or zeros and ones.
"""
function make_lonlat_mask(
    var;
    set_to_val = nothing,
    true_val = NaN,
    false_val = 1.0,
)
    # Check if lon and lat are the only dimensions
    has_longitude(var) ||
        error("Longitude does not exist as a dimension in var")
    has_latitude(var) || error("Latitude does not exist as a dimension in var")
    length(var.dims) == 2 ||
        error("Number of dimensions ($(length(var.dims))) is not two")

    # Preprocess data for the mask if needed
    if !isnothing(set_to_val)
        true_indices = findall(set_to_val, var.data)
        false_indices = findall(!set_to_val, var.data)
        cleaned_up_data = var.data |> copy
        cleaned_up_data[true_indices] .= true_val
        cleaned_up_data[false_indices] .= false_val
        var = OutputVar(
            var.attributes |> deepcopy,
            var.dims |> deepcopy,
            var.dim_attributes |> deepcopy,
            cleaned_up_data,
        )
    end

    # Make a dummy OutputVar so we can reorder on it
    lon = [42.0]
    lat = [42.0]
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict{String, String}()
    dim_attribs = OrderedDict{String, Dict{String, String}}()
    reordered_var = OutputVar(attribs, dims, dim_attribs, data)

    mask_var = reordered_as(var, reordered_var)

    return apply_lonlat_mask(input_var) = begin
        # Check if lon and lat exist
        has_longitude(input_var) ||
        error("Longitude does not exist as a dimension in var")
        has_latitude(input_var) ||
        error("Latitude does not exist as a dimension in var")

        # Resample so that the mask match up with the grid of var
        # Round because linear resampling is done and we want the mask to be only ones and zeros
        mask_arr =
            [
                mask_var(pt) for pt in Base.product(
                    input_var.dims[longitude_name(input_var)],
                    input_var.dims[latitude_name(input_var)],
                )
            ] .|> round

        # Reshape data for broadcasting
        lon_idx = input_var.dim2index[longitude_name(input_var)]
        lat_idx = input_var.dim2index[latitude_name(input_var)]
        lon_length = input_var.dims[longitude_name(mask_var)] |> length
        lat_length = input_var.dims[latitude_name(mask_var)] |> length
        if lon_idx > lat_idx
            mask_arr = transpose(mask_arr)
        end
        size_to_reshape = (
            if i == lon_idx
                lon_length
            elseif i == lat_idx
                lat_length
            else
                1
            end for i in 1:ndims(input_var.data)
        )

        # Mask data
        mask_arr = reshape(mask_arr, size_to_reshape...)
        masked_data = input_var.data .* mask_arr

        # Remake OutputVar with new data
        ret_attribs = deepcopy(input_var.attributes)
        ret_dims = deepcopy(input_var.dims)
        ret_dim_attributes = deepcopy(input_var.dim_attributes)
        return OutputVar(ret_attribs, ret_dims, ret_dim_attributes, masked_data)
    end
end

"""
    replace(var::OutputVar, old_new::Pair...)

Return a `OutputVar` where, for each pair `old  => new`, all occurences of `old` are
replaced by `new` in `Var.data`

This function is useful if there are `NaN`s or `missing` values in the data. For instance,
you want to use the ocean mask, but there are `NaN`s in the ocean. You can replace all the
`NaN` and `missing` values with 0.0 and apply the ocean mask afterward.
"""
function Base.replace(var::OutputVar, old_new::Pair...)
    replaced_data = replace(var.data, old_new...)

    # Remake OutputVar with replaced_data
    ret_attribs = deepcopy(var.attributes)
    ret_dims = deepcopy(var.dims)
    ret_dim_attributes = deepcopy(var.dim_attributes)
    return OutputVar(ret_attribs, ret_dims, ret_dim_attributes, replaced_data)
end

"""
    overload_binary_op(op)

Add methods to overload the given binary `op`erator for `OutputVars` and `Real`s.

Attributes that are not `short_name`, `long_name`, are discarded in the process.
"""
macro overload_binary_op(op)
    quote
        function Base.$op(x::OutputVar, y::OutputVar)
            arecompatible(x, y) || error("Input OutputVars are not compatible")

            ret_attributes = Dict{String, Any}()

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(x.attributes, attr) && haskey(y.attributes, attr)
                    ret_attributes[attr] = string(
                        x.attributes[attr],
                        " ",
                        string($op),
                        " ",
                        y.attributes[attr],
                    )
                end
            end

            ret_dims = x.dims
            ret_dim_attributes = x.dim_attributes

            ret_data = @. $op(x.data, y.data)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
        function Base.$op(x::OutputVar, y::Real)
            ret_attributes = empty(x.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(x.attributes, attr)
                    ret_attributes[attr] =
                        string(x.attributes[attr], " ", string($op), " ", y)
                end
            end

            ret_dims = deepcopy(x.dims)
            ret_dim_attributes = deepcopy(x.dim_attributes)

            ret_data = @. $op(x.data, y)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
        function Base.$op(x::Real, y::OutputVar)
            ret_attributes = empty(y.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(y.attributes, attr)
                    ret_attributes[attr] =
                        string(x, " ", string($op), " ", y.attributes[attr])
                end
            end

            ret_dims = deepcopy(y.dims)
            ret_dim_attributes = deepcopy(y.dim_attributes)

            ret_data = @. $op(x, y.data)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
    end
end

@overload_binary_op (+)
@overload_binary_op (-)
@overload_binary_op (*)
@overload_binary_op (/)

include("outvar_dimensions.jl")
include("constants.jl")

end
