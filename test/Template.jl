module Template

import ClimaAnalysis
import OrderedCollections: OrderedDict

struct TemplateVar
    attributes_fn::Vector{Function}
    dims_fn::Vector{Function}
    dim_attributes_fn::Vector{Function}
    data_fn::Vector{Function}
    # data_type::Base.RefValue{DataType} # I think I am going to remove this
    # TODO: Maybe add a way to keep track of which dimensions are added so that one can reorder them?
end

"""
    TemplateVar()

Intialize a `TemplateVar`.
"""
function TemplateVar()
    return TemplateVar(
        Function[],
        Function[],
        Function[],
        Function[],
        # Ref(Float64),
    )
end

function initialize(var::TemplateVar)
    # Add attributes
    attribs = Dict{String, Any}()
    for fn! in var.attributes_fn
        fn!(attribs)
    end

    # Add dimensions
    dims = OrderedDict{String, AbstractArray}()
    for fn! in var.dims_fn
        fn!(dims)
    end

    # Add attributes for dimensions
    dim_attribs = OrderedDict(dim_name => Dict() for dim_name in keys(dims))
    for fn! in var.dim_attributes_fn
        fn!(dim_attribs)
    end

    # Add data (TODO: Still a little unsure about this)
    dim_sizes = Tuple(length(dim) for dim in values(dims))
    if length(var.data_fn) == 0
        data = zeros(dim_sizes...)
    else
        data = var.data_fn[end](; dim_sizes = dim_sizes) # I think the most conservative thing to do is to always use the last function
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
        add_attribs = let attribs = attribs
            function add_attribs(attribs_dict)
                for (key, value) in attribs
                    attribs_dict[string(key)] = value
                end
                return nothing
            end
        end
        push!(var.attributes_fn, add_attribs)
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

function add_dim!(var::TemplateVar, dim_name, dim_array; dim_attribs...)
    push!(var.dims_fn, dims -> dims[dim_name] = dim_array)
    if !isempty(dim_attribs)
        add_dim_attribs = let dim_attribs = dim_attribs, dim_name = dim_name
            function add_dim_attribs(dim_attribs_dict)
                for (key, value) in dim_attribs
                    dim_attribs_dict[dim_name][string(key)] = value
                end
                return nothing
            end
        end
        push!(var.dim_attributes_fn, add_dim_attribs)
    end
    return var
end

function ones_data(; data_type = Float64)
    return var -> ones_data!(var, data_type = data_type)
end

function ones_data!(var::TemplateVar; data_type = Float64)
    push!(var.data_fn, (; dim_sizes) -> ones(data_type, dim_sizes...))
    return var
end

function zeros_data(; data_type = Float64)
    return var -> zeros_data!(var, data_type = data_type)
end

function zeros_data!(var::TemplateVar; data_type = Float64)
    push!(var.data_fn, (; dim_sizes) -> zeros(data_type, dim_sizes...))
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
            (; dim_sizes) ->
                reshape(data_type(1):data_type(prod(dim_sizes)), dim_sizes),
        )
    else
        push!(
            var.data_fn,
            (; dim_sizes) -> reshape(
                collect(data_type(1):data_type(prod(dim_sizes))),
                dim_sizes,
            ),
        )
    end
    return var
end
end

var =
    Template.TemplateVar() |>
    Template.add_attribs(; long_name = "hi") |>
    Template.add_lat_dim(; units = "deg") |>
    Template.add_lon_dim(; units = "deg") |>
    Template.ones_data() |>
    Template.initialize


lat = collect(range(-90.0, 90.0, 181))
lon = collect(range(180.0, 180.0, 360))
data = ones(length(lat), length(lon))
dims = OrderedDict(["lat" => lat, "lon" => lon])
attribs = Dict("long_name" => "hi")
dim_attribs = OrderedDict([
    "lat" => Dict("units" => "deg"),
    "lon" => Dict("units" => "deg"),
])
var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
