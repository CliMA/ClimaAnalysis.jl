import NCDatasets
import OrderedCollections: OrderedDict

import Statistics: mean

struct OutputVar{
    T <: AbstractArray,
    A <: AbstractArray,
    B,
    C,
    FP <: Union{String, Nothing},
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
        dims = map(nc.dim) do dim
            dim_name = dim[1]
            return dim_name => Array(nc[dim_name])
        end |> OrderedDict
        var_name = pop!(setdiff(keys(nc), keys(dims)))
        attribs = Dict(nc[var_name].attrib)
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

```jldoctest
julia> import Statistics: mean
julia> long = 0.:180. |> collect
julia> lat = 0.:90. |> collect
julia> data = reshape(1.:91*181., (181, 91))
julia> dims = Dict(["lat" => lat, "long" => long])
julia> var = OutputVar(dims, data)

julia> reduce_over(mean, "lat", var)
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
    data = reduction(var.data, dims = dim_index) |> Utils.squeeze

    # If we reduce over a dimension, we have to remove it
    dims = copy(var.dims)
    dim_attributes = copy(var.dim_attributes)
    pop!(dims, dim)
    haskey(var.dim_attributes, dim) && pop!(dim_attributes, dim)
    return OutputVar(var.attributes, dims, dim_attributes, data, var.file_path)
end

"""
    average_lat(var::OutputVar)

Return a new OutputVar where the values on the latitudes are averaged arithmetically.
"""
average_lat(var) = _reduce_over(mean, "lat", var)

"""
    average_lon(var::OutputVar)

Return a new OutputVar where the values on the longitudes are averaged arithmetically.
"""
average_lon(var) = _reduce_over(mean, "lon", var)

"""
    average_time(var::OutputVar)

Return a new OutputVar where the values are averaged arithmetically in time.
"""
average_time(var) = _reduce_over(mean, "time", var)

"""
    slice_general(var::OutputVar, val, dim_name)

Return a new OutputVar by selecting the available index closest to the given `val` for the
given dimension
"""
function slice_general(var, val, dim_name)
    nearest_index = Utils.nearest_index(var.dims[dim_name], val)
    _slice_over(data; dims) = selectdim(data, dims, nearest_index)
    _reduce_over(_slice_over, dim_name, var)
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
