import NCDatasets

struct OutputVar{T <: AbstractArray, A <: AbstractArray, B}

    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::Dict{String, T}

    "Array that contains all the data"
    var::A

    "File associated to this variable"
    file_path::String
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
        end |> Dict
        var_name = pop!(setdiff(keys(nc), keys(dims)))
        attribs = Dict(nc[var_name].attrib)
        var = Array(nc[var_name])
        return OutputVar(attribs, dims, var, path)
    end
end
