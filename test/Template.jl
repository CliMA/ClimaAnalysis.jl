module Template

import ClimaAnalysis
import OrderedCollections: OrderedDict

"A struct to hold a lazy evaluation"
struct LazyEval

    "Function to do lazily evaluate"
    f::Function

    "Positional arguments to pass to f"
    args::Tuple

    "Keyword arguments to pass to f"
    kwargs::Tuple
end

"""
    make_lazy(f, args...; kwargs...)

Delay function call to `f`.

Not all positional arguments and keyword arguments need to be supplied to `f`
when making the function call lazy and any remaining arguments can be passed
with with `call`.
"""
function make_lazy(f, args...; kwargs...)
    return LazyEval(f, tuple(args...), tuple(kwargs...))
end

"""
    call(comp::LazyEval, args...; kwargs...)

Evaluate `comp` with positional arguments and keyword arguments stored in
`comp`.

Any other arguments needed to evaluate `comp` can be passed with `args` and `kwargs`.
"""
function call(comp::LazyEval, args...; kwargs...)
    return comp.f(comp.args..., args...; comp.kwargs..., kwargs...)
end

"A struct to hold an uninitialized OutputVar"
struct TemplateVar
    attributes_fn::Vector{LazyEval}
    # TODO: Make dims and dim_attributes_fn into one thing
    dims_fn::Vector{LazyEval}
    dim_attributes_fn::Vector{LazyEval}
    data_fn::Vector{LazyEval}
end

"""
    TemplateVar()

Intialize a `TemplateVar`.

A `TemplateVar` is an uninitialized `OutputVar`.
"""
function TemplateVar()
    return TemplateVar(LazyEval[], Function[], Function[], Function[])
end

"""
    initialize(var::TemplateVar)

Initialize a `TemplateVar`.
"""
function initialize(var::TemplateVar)
    # Add attributes
    attribs = Dict{String, Any}()
    for fn! in var.attributes_fn
        call(fn!, attribs) # Can make this be a functor
    end

    # Add dimensions
    # TODO: Maybe enforce that there are only be one of each dimension and replace the lazyeval
    # when adding dimensions
    # For example, adding two dimensions can be weird
    # easiest way would be to add a tag, so maybe (dim_name, dims_fn)
    # but then, the attributes would be out of sync so we store everything together
    # so a vector of (dim_name, dims_fn, dims_attribs_fn)
    dims = OrderedDict{String, AbstractArray}()
    for fn! in var.dims_fn
        call(fn!, dims)
    end

    # Add attributes for dimensions
    dim_attribs = OrderedDict(dim_name => Dict() for dim_name in keys(dims))
    for fn! in var.dim_attributes_fn
        call(fn!, dim_attribs)
    end

    # Add data
    dim_sizes = Tuple(length(dim) for dim in values(dims))
    if length(var.data_fn) == 0
        # default is the zero matrix
        data = zeros(dim_sizes...)
    else
        # If there are more than one LazyComp for initializing data, use the last one
        # It is not easy to modify this since you may want to use an UnitRange for the data
        # which doesn't support indexed assignment
        data = call(var.data_fn[end]; dim_sizes = dim_sizes)
    end

    return ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
end

"""
    add_attribs(; attribs...)

Return the function `add_attribs` with the keyword arguments `attribs`.

Intended to be used with the pipe operator (`|>`).
"""
function add_attribs(; attribs...)
    return var -> add_attribs!(var; attribs...)
end

"""
    add_attribs!(var::TemplateVar; attribs...)

Add attributes and return `var`.

Intended to be used with function composition.
"""
function add_attribs!(var::TemplateVar; attribs...)
    if !isempty(attribs)
        function add_attribs(attribs_dict, attribs)
            for (key, value) in attribs
                attribs_dict[string(key)] = value
            end
            return nothing
        end
        push!(var.attributes_fn, make_lazy(add_attribs, attribs))
    end
    return var
end

"""
    add_lat_dim(; dim_name = "latitude",
                  dim_array = collect(range(-90.0, 90.0, 181)),
                  units = "degrees",
                  lat_attribs...)

Return the function `add_lat_dim!` with all keyword arguments filled out.

Intended to be used with the pipe operator (`|>`).
"""
function add_lat_dim(;
    dim_name = "latitude",
    dim_array = collect(range(-90.0, 90.0, 181)),
    units = "degrees",
    lat_attribs...,
)
    dim_name in ClimaAnalysis.Var.LATITUDE_NAMES ||
        error("$dim_name is not a name for the latitude dimension")
    return var -> add_lat_dim!(
        var;
        dim_name = dim_name,
        dim_array = dim_array,
        units = units,
        lat_attribs...,
    )
end

"""
    add_lat_dim!(var::TemplateVar;
                 dim_name = "latitude",
                 dim_array = collect(range(-90.0, 90.0, 181)),
                 units = "degrees", # todo: Change this
                 lat_attribs...,
)

Add latitude dimension and attributes and return `var`.

Intended to be used with function composition.
"""
function add_lat_dim!(
    var::TemplateVar;
    dim_name = "latitude",
    dim_array = collect(range(-90.0, 90.0, 181)),
    units = "degrees", # todo: Change this
    lat_attribs...,
)
    dim_name in ClimaAnalysis.Var.LATITUDE_NAMES ||
        error("$dim_name is not a name for the latitude dimension")
    return add_dim!(var, dim_name, dim_array; units = units, lat_attribs...)
end

"""
    add_lon_dim(; dim_name = "longitude",
                dim_array = collect(range(-180.0, 180.0, 361)),
                units = "degrees", # TODO: Change this
                lon_attribs...)

Return the function `add_lon_dim!` with all keyword arguments filled out.

Intended to be used with the pipe operator (`|>`).
"""
function add_lon_dim(;
    dim_name = "longitude",
    dim_array = collect(range(-180.0, 180.0, 361)),
    units = "degrees", # TODO: Change this
    lon_attribs...,
)
    dim_name in ClimaAnalysis.Var.LONGITUDE_NAMES ||
        error("$dim_name is not a name for the longitude dimension")
    return var -> add_lon_dim!(
        var;
        dim_name = dim_name,
        dim_array = dim_array,
        units = units,
        lon_attribs...,
    )
end

"""
    add_lon_dim!(var::TemplateVar;
                 dim_name = "longitude",
                 dim_array = collect(range(-180.0, 180.0, 361)),
                 units = "degrees", # todo: Change this
                 lon_attribs...)

Add longitude dimension and attributes and return `var`.

Intended to be used with function composition.
"""
function add_lon_dim!(
    var::TemplateVar;
    dim_name = "longitude",
    dim_array = collect(range(-180.0, 180.0, 361)),
    units = "degrees", # todo: Change this
    lon_attribs...,
)
    dim_name in ClimaAnalysis.Var.LONGITUDE_NAMES ||
        error("$dim_name is not a name for the longitude dimension")
    return add_dim!(var, dim_name, dim_array; units = units, lon_attribs...)
end

function add_time_dim(;
    dim_name = "time",
    dim_array = collect(0.0:10.0),
    units = "seconds",
    time_attribs...,
)
    dim_name in ClimaAnalysis.Var.TIME_NAMES ||
        error("$dim_name is not a name for the time dimension")
    return var -> add_time_dim!(
        var;
        dim_name = dim_name,
        dim_array = dim_array,
        units = units,
        time_attribs...,
    )
end

function add_time_dim!(
    var::TemplateVar;
    dim_name = "time",
    dim_array = collect(0.0:10.0),
    units = "seconds",
    time_attribs...,
)
    dim_name in ClimaAnalysis.Var.TIME_NAMES ||
        error("$dim_name is not a name for the time dimension")
    return add_dim!(var, dim_name, dim_array; units = units, time_attribs...)
end

# TODO: Add pressure dimension

function add_dim!(var::TemplateVar, dim_name, dim_array; dim_attribs...)
    push!(
        var.dims_fn,
        make_lazy(
            (dim_name, dim_array, dims) -> dims[dim_name] = dim_array,
            dim_name,
            dim_array,
        ),
    )
    if !isempty(dim_attribs)
        function add_dim_attribs(dim_name, dim_attribs, dim_attribs_dict)
            for (key, value) in dim_attribs
                dim_attribs_dict[dim_name][string(key)] = value
            end
            return nothing
        end
        push!(
            var.dim_attributes_fn,
            make_lazy(add_dim_attribs, dim_name, dim_attribs),
        )
    end
    return var
end

"""
    ones_data(; data_type = Float64)
"""
function ones_data(; data_type = Float64)
    return var -> ones_data!(var, data_type = data_type)
end

"""
    ones_data!(var::TemplateVar; data_type = Float64)
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

TODO
"""
function zeros_data(; data_type = Float64)
    return var -> zeros_data!(var, data_type = data_type)
end

"""
    zeros_data!(var::TemplateVar; data_type = Float64)

TODO
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
    ones_to_n_data(; data_type = Float64, lazy = true)
"""
function ones_to_n_data(; data_type = Float64, lazy = true)
    return var ->
        ones_to_n_data!(var::TemplateVar, data_type = data_type, lazy = lazy)
end

"""
    ones_to_n_data!(var; data_type = Float64, lazy = true)

Add data of `1:n` to `var` where `n` is the product of the sizes of the dimensions.

If `lazy = true`, then `collect` is not called on the data and if `lazy = false`, then
`collect` is called on the data.
"""
function ones_to_n_data!(var; data_type = Float64, lazy = true)
    if lazy
        push!(
            var.data_fn,
            make_lazy(
                (data_type; dim_sizes) ->
                    reshape(data_type(1):data_type(prod(dim_sizes)), dim_sizes),
                data_type,
            ),
        )
    else
        push!(
            var.data_fn,
            make_lazy(
                (data_type; dim_sizes) -> reshape(
                    collect(data_type(1):data_type(prod(dim_sizes))),
                    dim_sizes,
                ),
                data_type,
            ),
        )
    end
    return var
end

"""
    make_template_var(dims::String...)

Make a `TemplateVar` with the specified dimensions in `dims`.

The order of `dims` dicates the order of the dimensions in the `OutputVar` after
initialization.
"""
function make_template_var(dims::String...)
    # TODO: Convert to conventional dim names
    dims_to_fn = Dict(
        "time" => add_time_dim!,
        "lon" => add_lon_dim!,
        "lat" => add_lat_dim!,
    )

    # Maybe use reduce here or some functional programming thing
    var = TemplateVar()
    for dim in dims
        var = dims_to_fn[dim](var)
    end
    return var
end

end

var = Template.make_template_var("lat", "lon", "time") |> Template.initialize

# var =
#     Template.TemplateVar() |>
#     Template.add_attribs(; long_name = "hi") |>
#     Template.add_time_dim(; dim_array = collect(1.0:30.0)) |>
#     Template.add_lon_dim(; units = "deg") |>
#     Template.add_lat_dim(; units = "deg") |>
#     Template.ones_to_n_data(lazy = true) |>
#     Template.initialize


# lat = collect(range(-90.0, 90.0, 181))
# lon = collect(range(-180.0, 180.0, 360))
# data = ones(length(lat), length(lon))
# dims = OrderedDict(["lat" => lat, "lon" => lon])
# attribs = Dict("long_name" => "hi")
# dim_attribs = OrderedDict([
#     "lat" => Dict("units" => "deg"),
#     "lon" => Dict("units" => "deg"),
# ])
# var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
