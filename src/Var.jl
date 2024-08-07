module Var

import NCDatasets
import OrderedCollections: OrderedDict

import Interpolations as Intp
import Statistics: mean
import NaNStatistics: nanmean

import ..Utils: nearest_index, seconds_to_prettystr, squeeze

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
    range_dim

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
    read_var(path::String)

Read a variable in the given NetCDF file.

Example
=========

```julia
simdir = SimDir("my_output")
read_var(simdir.variable_paths["hu"]["inst"])
```
"""
function read_var(path::String)
    NCDatasets.NCDataset(path) do nc
        # First, we have to identify the name of the variable by finding what is
        # not a dimension
        unordered_dims = NCDatasets.dimnames(nc)
        var_name = pop!(setdiff(keys(nc), unordered_dims))
        dims = map(NCDatasets.dimnames(nc[var_name])) do dim_name
            return dim_name => Array(nc[dim_name])
        end |> OrderedDict
        attribs = Dict(k => v for (k, v) in nc[var_name].attrib)
        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[var_name])
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
    get(var.attributes, "units", "")
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
    haskey(var.dims, "z") ||
        error("Variable does not have an altitude dimension")

    return length(size(var.dims["z"])) == 1
end

function Base.copy(var::OutputVar)
    fields = fieldnames(OutputVar)
    # We have nested containers and we have to make sure we hand ownership off,
    # so we deepcopy
    return OutputVar([deepcopy(getfield(var, field)) for field in fields]...)
end

"""
    _reduce_over(reduction::F, dim::String, var::OutputVar)

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
    data = squeeze(reduction(var.data, dims = dim_index), dims = (dim_index,))

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
    _update_long_name_average!(reduced_var, var, latitude_name(var))


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
    _update_long_name_average!(reduced_var, var, longitude_name(var))
    return reduced_var
end

"""
    average_x(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `x` dimension are averaged arithmetically.
"""
function average_x(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "x", var)
    _update_long_name_average!(reduced_var, var, "x")
    return reduced_var
end

"""
    average_y(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `y` dimension are averaged arithmetically.
"""
function average_y(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "y", var)
    _update_long_name_average!(reduced_var, var, "y")
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
    _update_long_name_average!(reduced_var, var, time_name(var))
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
    _update_long_name_average!(reduced_var::OutputVar, var::OutputVar, dim_name)

Used by averaging functions to update the long name of `reduced_var` by describing what the
average is being taken over and the associated units. 
"""
function _update_long_name_average!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_name,
)
    dim_of_units = dim_units(var, dim_name)
    first_elt, last_elt = range_dim(var, dim_name)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged over $dim_name ($first_elt to $last_elt$dim_of_units)"
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
