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
    _data_at_dim_vals

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
    resampled_as,
    has_units,
    convert_units,
    integrate_lonlat,
    integrate_lon,
    integrate_lat,
    isempty,
    split_by_season

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

function OutputVar(attribs, dims, dim_attribs, data)
    index2dim = keys(dims) |> collect
    dim2index =
        Dict([dim_name => index for (index, dim_name) in enumerate(keys(dims))])

    # We can only create interpolants when we have 1D dimensions
    if isempty(index2dim) ||
       any(d -> ndims(d) != 1 || length(d) == 1, values(dims))
        itp = nothing
    else
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

        dims_tuple = tuple(values(dims)...)

        # TODO: Make this lazy: we should compute the spline the first time we use
        # it, not when we create the object
        itp = Intp.extrapolate(
            Intp.interpolate(dims_tuple, data, Intp.Gridded(Intp.Linear())),
            Intp.Throw(),
        )
    end

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
    OutputVar(path, short_name = nothing)

Read the NetCDF file in `path` as an `OutputVar`.

If `short_name` is `nothing`, automatically find the name.
"""
function OutputVar(path::String, short_name = nothing)
    return read_var(path; short_name)
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

Extrapolation is now allowed and will throw a `BoundsError`.

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
    return (x.dims == y.dims) && (x.dim_attributes == y.dim_attributes)
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
    resampled_as(src_var::OutputVar, dest_var::OutputVar)

Resample `data` in `src_var` to `dims` in `dest_var`.

The resampling performed here is a 1st-order linear resampling.
"""
function resampled_as(src_var::OutputVar, dest_var::OutputVar)
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

            ret_data = $op(x.data, y.data)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
        function Base.$op(x::OutputVar, y::Real)
            ret_var = copy(x)
            @. ret_var.data = $op(x.data, y)
            empty!(ret_var.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(x.attributes, attr)
                    ret_var.attributes[attr] =
                        string(x.attributes[attr], " ", string($op), " ", y)
                end
            end

            return ret_var
        end
        function Base.$op(x::Real, y::OutputVar)
            ret_var = copy(y)
            @. ret_var.data = $op(x, y.data)
            empty!(ret_var.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(y.attributes, attr)
                    ret_var.attributes[attr] =
                        string(x, " ", string($op), " ", y.attributes[attr])
                end
            end

            return ret_var
        end
    end
end

@overload_binary_op (+)
@overload_binary_op (-)
@overload_binary_op (*)
@overload_binary_op (/)

include("outvar_dimensions.jl")

end
