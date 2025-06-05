module Template

import ClimaAnalysis
import OrderedCollections: OrderedDict

# Alternative approaches to lazy evaluation exist, such as using broadcasted objects (see
# LazyBroadcast.jl or DontMaterialize.jl). However, this approach was not pursued because
# handling dictionaries is cumbersome (require using Ref), and requiring @. for every lazy
# evaluation makes reading the code more difficult.
"A struct to hold a lazy evaluation"
struct LazyEval{F <: Function, T <: Tuple, NT <: NamedTuple}
    "Function to lazily evaluate"
    f::F

    "Positional arguments to pass to f"
    args::T

    "Keyword arguments to pass to f"
    kwargs::NT
end

"""
    make_lazy(f, args...; kwargs...)

Delay function call to `f`.

Not all positional arguments and keyword arguments need to be supplied to `f`
with this function. Any additional arguments and keyword arguments can be
supplied when evaluating `LazyEval`.
"""
function make_lazy(f, args...; kwargs...)
    return LazyEval(f, args, (; kwargs...))
end

"""
    (l::LazyEval)(args...; kwargs...)

Evaluate a `LazyEval` object by calling its stored function `f` with the stored arguments
and any additional arguments provided.
"""
function (l::LazyEval)(args...; kwargs...)
    l.f(l.args..., args...; merge(l.kwargs, kwargs)...)
end

"A struct to hold an uninitialized OutputVar"
struct TemplateVar
    "Stores attribute functions to be added to the OutputVar"
    attributes_fn::Vector{LazyEval}

    "Maps each conventional dimension name to a tuple of lazy functions: one for generating
    the dimension values and one for generating the dimension attributes"
    dims_fn::OrderedDict{String, Tuple{LazyEval, LazyEval}}

    "Stores lazy functions that generate data"
    data_fn::Vector{LazyEval}
end

"""
    TemplateVar()

Intialize a `TemplateVar`.

A `TemplateVar` is an uninitialized `OutputVar`.

# Examples

This example shows how to add attributes, dimensions, and data to a `TemplateVar` and
initialize it.
```julia
using ClimaAnalysis.Template
var =
    Template.TemplateVar() |>
    Template.add_attribs(long_name = "Test", short_name = "test") |>
    Template.add_time_dim(name = "t", dim = collect(0.0:2.0)) |>
    Template.add_lon_dim(units = "degrees") |>
    Template.add_lat_dim(axis = "Y") |>
    Template.one_to_n_data(collected = true) |>
    Template.initialize
```
"""
function TemplateVar()
    return TemplateVar(
        LazyEval[],
        OrderedDict{String, Tuple{LazyEval, LazyEval}}(),
        LazyEval[],
    )
end

"""
    initialize(var::TemplateVar)

Initialize a `TemplateVar` into an `OutputVar`.
"""
function initialize(var::TemplateVar)
    # Add attributes
    attribs =
        isempty(var.attributes_fn) ? Dict{String, Any}() :
        merge((make_attrib() for make_attrib in var.attributes_fn)...)

    # Add dimensions and dimension attributes
    if isempty(var.dims_fn)
        dims = OrderedDict{String, Vector}()
        dim_attribs = OrderedDict{String, Dict}()
    else
        dims = OrderedDict(
            (make_dim() for (make_dim, _) in values(var.dims_fn))...,
        )
        dim_attribs = OrderedDict(
            (
                make_dim_attrib() for
                (_, make_dim_attrib) in values(var.dims_fn)
            )...,
        )
    end

    # Add data
    dim_sizes = Tuple(length(dim) for dim in values(dims))
    if length(var.data_fn) == 0 && length(dims) == 0
        data = empty([0.0])
    elseif length(var.data_fn) == 0
        # default is the one to n array
        data = reshape(1:prod(dim_sizes), dim_sizes...)
    else
        # If there are more than one LazyEval for initializing data, use the last one. It is
        # not easy to modify the data in-place, since you may want to use an UnitRange for
        # the data which doesn't support indexed assignment
        data = var.data_fn[end](; dim_sizes = dim_sizes)
    end
    return ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
end

"""
    add_attribs(; attribs...)

Return the function `add_attribs` with the keyword arguments `attribs`.

Designed to be used with the pipe operator (`|>`).
"""
function add_attribs(; attribs...)
    return var -> add_attribs!(var; attribs...)
end

"""
    add_attribs!(var::TemplateVar; attribs...)

Add attributes and return `var`.

Designed to be used with function composition.
"""
function add_attribs!(var::TemplateVar; attribs...)
    if !isempty(attribs)
        function add_attribs(attribs)
            return Dict(string(key) => value for (key, value) in attribs)
        end
        push!(var.attributes_fn, make_lazy(add_attribs, attribs))
    end
    return var
end

"""
    macro generate_add_dim(dim_name,
                           default_dim_name,
                           dim,
                           units,
                           possible_dim_names)

Macro to generate functions `add_DIMNAME_dim` and `add_DIMNAME_dim!`.
"""
macro generate_add_dim(
    dim_name,
    default_dim_name,
    dim,
    units,
    possible_dim_names,
)
    func_name! = Symbol("add_", dim_name, "_dim!")
    func_name = Symbol("add_", dim_name, "_dim")
    conventional_name = ClimaAnalysis.conventional_dim_name(default_dim_name)
    return quote
        # Define add_dim
        function $(esc(func_name))(;
            name = $default_dim_name,
            dim = $dim,
            units = $units,
            dim_attribs...,
        )
            name in $possible_dim_names || error(
                "$name is not a name for the $($conventional_name) dimension",
            )
            return var -> $(esc(func_name!))(
                var;
                name = name,
                dim = dim,
                units = units,
                dim_attribs...,
            )
        end

        # Define add_dim!
        function $(esc(func_name!))(
            var::TemplateVar;
            name = $dim_name,
            dim = $dim,
            units = $units,
            dim_attribs...,
        )
            name in $possible_dim_names || error(
                "$name is not a name for the $($conventional_name) dimension",
            )
            return add_dim!(var, name, dim; units = units, dim_attribs...)
        end
    end
end

@generate_add_dim(
    "z",
    "z",
    collect(0.0:10.0),
    "m",
    ClimaAnalysis.Var.ALTITUDE_NAMES
)
@generate_add_dim(
    "lat",
    "latitude",
    collect(range(-90.0, 90.0, 181)),
    "degrees_east",
    ClimaAnalysis.Var.LATITUDE_NAMES
)
@generate_add_dim(
    "lon",
    "longitude",
    collect(range(-180.0, 180.0, 361)),
    "degrees_north",
    ClimaAnalysis.Var.LONGITUDE_NAMES
)
@generate_add_dim(
    "time",
    "time",
    collect(0.0:10.0),
    "seconds",
    ClimaAnalysis.Var.TIME_NAMES
)
@generate_add_dim(
    "pfull",
    "pfull",
    collect(0.0:10.0),
    "Pa",
    ClimaAnalysis.Var.PRESSURE_NAMES
)

"""
    add_dim(dim_name, dim; dim_attribs...)

Add `dim_name` dimension with array `dim` and attributes `dim_attribs`.

There are the convenience functions `add_time_dim`, `add_lon_dim`, `add_lat_dim`,
`add_pfull_dim`, and `add_z_dim` which uses default values for `dim_name`, `dim`, and
`units`.
"""
function add_dim(dim_name, dim; dim_attribs...)
    return var -> add_dim!(var, dim_name, dim; dim_attribs...)
end

"""
    add_dim!(var::TemplateVar, dim_name, dim; dim_attribs...)

Add `dim_name` dimension with array `dim` and attributes `dim_attribs`.

There are the convenience functions `add_time_dim!`, `add_lon_dim!`, `add_lat_dim!`,
`add_pfull_dim!`, and `add_z_dim!` which uses default values for `dim_name`, `dim`, and
`units`.
"""
function add_dim!(var::TemplateVar, dim_name, dim; dim_attribs...)
    function add_dim_attribs(dim_name, dim_attribs)
        dim_attribs = Dict(string(key) => value for (key, value) in dim_attribs)
        return dim_name => dim_attribs
    end

    conventional_name = ClimaAnalysis.conventional_dim_name(dim_name)
    var.dims_fn[conventional_name] = (
        make_lazy((dim_name, dim) -> dim_name => dim, dim_name, dim),
        make_lazy(add_dim_attribs, dim_name, dim_attribs),
    )
    return var
end

"""
    ones_data(; data_type = Float64)

Add ones data to `TemplateVar`.

Designed to be used with the pipe operator (`|>`).
"""
function ones_data(; data_type = Float64)
    return var -> ones_data!(var, data_type = data_type)
end

"""
    ones_data!(var::TemplateVar; data_type = Float64)

Add ones data to `var`.

Designed to be used with function composition.
"""
function ones_data!(var::TemplateVar; data_type = Float64)
    push!(
        var.data_fn,
        make_lazy(
            (data_type; dim_sizes) -> ones(data_type, dim_sizes...),
            data_type,
        ),
    )
    return var
end

"""
    zeros_data(; data_type = Float64)

Add zeros data to `TemplateVar`.

Designed to be used with the pipe operator (`|>`).
"""
function zeros_data(; data_type = Float64)
    return var -> zeros_data!(var, data_type = data_type)
end

"""
    zeros_data!(var::TemplateVar; data_type = Float64)

Add zeros data to `var`.

Designed to be used with function composition.
"""
function zeros_data!(var::TemplateVar; data_type = Float64)
    push!(
        var.data_fn,
        make_lazy(
            (data_type; dim_sizes) -> zeros(data_type, dim_sizes...),
            data_type,
        ),
    )
    return var
end

"""
    one_to_n_data(; data_type = Float64, collected = false)

Add data of `1:n` to `TemplateVar` where `n` is the product of the sizes of the dimensions.

If `collected = false`, then `collect` is not called on the data and if `collected = true`,
then `collect` is called on the data.

Designed to be used with the pipe operator (`|>`).
"""
function one_to_n_data(; data_type = Float64, collected = false)
    return var ->
        one_to_n_data!(var, data_type = data_type, collected = collected)
end

"""
    one_to_n_data!(var::TemplateVar; data_type = Float64, collected = false)

Add data of `1:n` to `var` where `n` is the product of the sizes of the dimensions.

If `collected = true`, then `collect` is not called on the data and if `collected = false`,
then `collect` is called on the data.

Designed to be used with function composition.
"""
function one_to_n_data!(
    var::TemplateVar;
    data_type = Float64,
    collected = false,
)
    if collected
        push!(
            var.data_fn,
            make_lazy(
                (data_type; dim_sizes) -> collect(
                    reshape(data_type(1):data_type(prod(dim_sizes)), dim_sizes),
                ),
                data_type,
            ),
        )
    else
        push!(
            var.data_fn,
            make_lazy(
                (data_type; dim_sizes) ->
                    reshape(data_type(1):data_type(prod(dim_sizes)), dim_sizes),
                data_type,
            ),
        )
    end
    return var
end

"""
    add_data(; data)

Add `data` to `TemplateVar`.

Designed to be used with the pipe operator (`|>`).
"""
function add_data(; data)
    return var -> add_data!(var; data = data)
end

"""
    add_data!(var::TemplateVar; data)

Add `data` to `var`.

Designed to be used with function composition.
"""
function add_data!(var::TemplateVar; data)
    push!(var.data_fn, make_lazy((data; dim_sizes) -> data, data))
    return var
end


"""
    make_template_var(dims::String...; attribs...)

Make a `TemplateVar` with the specified dimensions in `dims` and additional `attribs`.

The order of `dims` dictates the order of the dimensions in the `OutputVar` after
initialization. The names of `dims` is not passed to the `TemplateVar`.
"""
function make_template_var(dims::String...; attribs...)
    dims = (ClimaAnalysis.conventional_dim_name(dim) for dim in dims)
    dims_to_fn = Dict(
        "time" => add_time_dim!,
        "longitude" => add_lon_dim!,
        "latitude" => add_lat_dim!,
        "pressure" => add_pfull_dim!,
        "altitude" => add_z_dim!,
    )

    var = TemplateVar()
    var = add_attribs!(var; attribs...)
    for dim in dims
        var = dims_to_fn[dim](var)
    end
    return var
end

end
