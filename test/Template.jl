module Template

import ClimaAnalysis
import OrderedCollections: OrderedDict

struct LazyComp
    f::Function
    args::Tuple
    kwargs::Tuple
end

function make_lazy(f, args...; kwargs...)
    return LazyComp(f, tuple(args...), tuple(kwargs...))
end

function call(comp::LazyComp, args...; kwargs...)
    return comp.f(args..., comp.args...; kwargs..., comp.kwargs...)
end

struct TemplateVar
    attributes_fn::Vector{LazyComp}
    dims_fn::Vector{LazyComp}
    dim_attributes_fn::Vector{LazyComp}
    data_fn::Vector{LazyComp}
    # TODO: Maybe add a way to keep track of which dimensions are added so that one can reorder them?
end

"""
    TemplateVar()

Intialize a `TemplateVar`.
"""
function TemplateVar()
    return TemplateVar(LazyComp[], Function[], Function[], Function[])
end

function initialize(var::TemplateVar)
    # Add attributes
    attribs = Dict{String, Any}()
    for fn! in var.attributes_fn
        call(fn!, attribs) # Can make this be a functor
    end

    # Add dimensions
    dims = OrderedDict{String, AbstractArray}()
    for fn! in var.dims_fn
        call(fn!, dims)
    end

    # Add attributes for dimensions
    dim_attribs = OrderedDict(dim_name => Dict() for dim_name in keys(dims))
    for fn! in var.dim_attributes_fn
        call(fn!, dim_attribs)
    end

    # Add data (TODO: Still a little unsure about this)
    dim_sizes = Tuple(length(dim) for dim in values(dims))
    if length(var.data_fn) == 0
        data = zeros(dim_sizes...)
    else
        data = call(var.data_fn[end]; dim_sizes = dim_sizes) # I think the most conservative thing to do is to always use the last function
    end

    return ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
end

"""
    add_attribs(; attribs...)

Return the function `add_attribs` with the keyword arguments `attribs`.

Intended to be used by the piping operator.
"""
function add_attribs(; attribs...)
    return var -> add_attribs!(var; attribs...)
end

"""
    add_attribs!(var; attribs...)

TODO: Suppose to be used with function composition
"""
function add_attribs!(var; attribs...)
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
    add_lat_dim

Return the function `add_lat_dim!` with all keyword arguments filled out.
Intended to be used with piping (|>)
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
    add_lat_dim!

Supposed to be used with functional composition
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

function add_dim!(var::TemplateVar, dim_name, dim_array; dim_attribs...)
    push!(
        var.dims_fn,
        make_lazy(
            (dims, dim_name, dim_array) -> dims[dim_name] = dim_array,
            dim_name,
            dim_array,
        ),
    )
    if !isempty(dim_attribs)
        function add_dim_attribs(dim_attribs_dict, dim_name, dim_attribs)
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

function ones_data(; data_type = Float64)
    return var -> ones_data!(var, data_type = data_type)
end

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

function zeros_data(; data_type = Float64)
    return var -> zeros_data!(var, data_type = data_type)
end

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

function ones_to_n_data(; data_type = Float64, lazy = true)
    return var ->
        ones_to_n_data!(var::TemplateVar, data_type = data_type, lazy = lazy)
end

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
end

var =
    Template.TemplateVar() |>
    Template.add_attribs(; long_name = "hi") |>
    Template.add_time_dim(; dim_array = collect(1.0:30.0)) |>
    Template.add_lon_dim(; units = "deg") |>
    Template.add_lat_dim(; units = "deg") |>
    Template.ones_to_n_data(lazy = true) |>
    Template.initialize


# lat = collect(range(-90.0, 90.0, 181))
# lon = collect(range(180.0, 180.0, 360))
# data = ones(length(lat), length(lon))
# dims = OrderedDict(["lat" => lat, "lon" => lon])
# attribs = Dict("long_name" => "hi")
# dim_attribs = OrderedDict([
#     "lat" => Dict("units" => "deg"),
#     "lon" => Dict("units" => "deg"),
# ])
# var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
