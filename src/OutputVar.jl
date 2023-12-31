import NCDatasets
import OrderedCollections: OrderedDict

import Statistics: mean

import .Utils: nearest_index

export OutputVar,
    read_var,
    average_lat,
    average_lon,
    average_time,
    slice,
    slice_time,
    slice_x,
    slice_x,
    slice_y,
    slice_lon,
    slice_lat,
    window

"""
    Representing an output variable
"""
struct OutputVar{
    T <: AbstractArray,
    A <: AbstractArray,
    B,
    C,
    FP <: Union{<:AbstractString, Nothing},
}

    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Array that contains all the data"
    data::A

    "File associated to this variable"
    file_path::FP

    "Dictionary that maps dimension name to its array index"
    dim2index::Dict{String, Int}

    "Array that maps name array index to the dimension name"
    index2dim::Vector{String}

end

function OutputVar(attribs, dims, dim_attribs, data, path)
    index2dim = keys(dims) |> collect
    dim2index =
        Dict([dim_name => index for (index, dim_name) in enumerate(keys(dims))])

    return OutputVar(
        attribs,
        OrderedDict(dims),
        OrderedDict(dim_attribs),
        data,
        path,
        dim2index,
        index2dim,
    )
end


function OutputVar(dims, data)
    return OutputVar(
        Dict{String, Any}(),
        dims,
        Dict{String, Any}(),
        data,
        nothing,
    )
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
        dims = map(NCDatasets.dimnames(nc)) do dim_name
            return dim_name => Array(nc[dim_name])
        end |> OrderedDict
        var_name = pop!(setdiff(keys(nc), keys(dims)))
        attribs = Dict(k => v for (k, v) in nc[var_name].attrib)
        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[var_name])
        return OutputVar(attribs, dims, dim_attribs, data, path)
    end
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
    data = Utils.squeeze(
        reduction(var.data, dims = dim_index),
        dims = (dim_index,),
    )

    # If we reduce over a dimension, we have to remove it
    dims = copy(var.dims)
    dim_attributes = copy(var.dim_attributes)
    pop!(dims, dim)
    haskey(var.dim_attributes, dim) && pop!(dim_attributes, dim)
    return OutputVar(
        copy(var.attributes),
        dims,
        dim_attributes,
        data,
        var.file_path,
    )
end

"""
    average_lat(var::OutputVar)

Return a new OutputVar where the values on the latitudes are averaged arithmetically.
"""
function average_lat(var)
    reduced_var = _reduce_over(mean, "lat", var)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged over latitudes"
    end
    return reduced_var
end

"""
    average_lon(var::OutputVar)

Return a new OutputVar where the values on the longitudes are averaged arithmetically.
"""
function average_lon(var)
    reduced_var = _reduce_over(mean, "lon", var)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged over longitudes"
    end

    return reduced_var
end

"""
    average_time(var::OutputVar)

Return a new OutputVar where the values are averaged arithmetically in time.
"""
function average_time(var)
    reduced_var = _reduce_over(mean, "time", var)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged over time"
    end

    return reduced_var
end

"""
    slice_general(var::OutputVar, val, dim_name)

Return a new OutputVar by selecting the available index closest to the given `val` for the
given dimension
"""
function slice_general(var, val, dim_name)
    haskey(var.dims, dim_name) ||
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")

    nearest_index_val = nearest_index(var.dims[dim_name], val)
    _slice_over(data; dims) = selectdim(data, dims, nearest_index_val)
    reduced_var = _reduce_over(_slice_over, dim_name, var)

    # Let's try adding this operation to the long_name, if possible (ie, if the correct
    # attributes are available)
    try
        dim_array = var.dims[dim_name]
        dim_units = var.dim_attributes[dim_name]["units"]
        cut_point = dim_array[nearest_index_val]
        reduced_var.attributes["long_name"] *= " $dim_name = $cut_point $dim_units"
    catch
    end
    return reduced_var
end

"""
    slice_time(var::OutputVar, time)

Return a new OutputVar by selecting the available snapshot closest to the given `time`.
"""
slice_time(var, time) = slice_general(var, time, "time")

"""
    slice_z(var::OutputVar, z)

Return a new OutputVar by selecting the available date closest to the given `z`.
"""
slice_z(var, z) = slice_general(var, z, "z")

"""
    slice_x(var::OutputVar, x)

Return a new OutputVar by selecting the available date closest to the given `x`.
"""
slice_x(var, x) = slice_general(var, x, "x")

"""
    slice_y(var::OutputVar, y)

Return a new OutputVar by selecting the available date closest to the given `y`.
"""
slice_y(var, y) = slice_general(var, y, "y")

"""
    slice_lon(var::OutputVar, lon)

Return a new OutputVar by selecting the available date closest to the given `lon`.
"""
slice_lon(var, lon) = slice_general(var, lon, "lon")

"""
    slice_lat(var::OutputVar, lat)

Return a new OutputVar by selecting the available date closest to the given `lat`.
"""
slice_lat(var, lat) = slice_general(var, lat, "lon")

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
        sliced_var = slice_general(sliced_var, val, String(dim_name))
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
    dim_attributes = copy(var.dim_attributes)
    return OutputVar(
        copy(var.attributes),
        dims,
        dim_attributes,
        reduced_data,
        var.file_path,
    )
end
