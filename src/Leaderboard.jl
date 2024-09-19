module Leaderboard

import OrderedCollections: OrderedDict
import NaNStatistics: nanmedian

export RMSEVariable,
    model_names,
    category_names,
    rmse_units

"""
    Holding root mean squared errors over multiple categories and models for a single
    variable.
"""
struct RMSEVariable{
    FT <: AbstractFloat,
    S1 <: AbstractString,
    S2 <: AbstractString,
}
    "Short name of variable of interest"
    short_name::AbstractString

    "Map model name (like `CliMA`) to index"
    model2index::OrderedDict{S1, Int}

    "Map category name (like `ANN` and seasons) to index"
    category2index::OrderedDict{S2, Int}

    "Arrray of RMSEs (rows correspond to model names and columns correspond to categories)"
    RMSEs::Array{FT, 2}

    "Map model_name to units"
    units::Dict{S1, String}
end

"""
    RMSEVariable(short_name::String,
                 model_names::Vector{String},
                 category_names::Vector{String},
                 RMSEs,
                 units::Dict)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, the categories in `category_names`, the root mean squared errors in `RMSEs`,
and `units`.
"""
function RMSEVariable(
    short_name::String,
    model_names::Vector{String},
    category_names::Vector{String},
    RMSEs,
    units::Dict,
)
    # Check if the dimensions of model_names and category_names match the dimensions of RMSE array
    (length(model_names), length(category_names)) != size(RMSEs) && error(
        "The size of RMSEs ($(size(RMSEs))) does not fit the number of model names ($(length(model_names))) and categories ($(length(category_names)))",
    )

    # Check if RMSE is negative
    any(RMSEs .< 0.0) && error("RMSEs cannot be negative")

    # Check for uniqueness
    length(unique(model_names)) == length(model_names) ||
        error("Model names are not unique")
    length(unique(category_names)) == length(category_names) ||
        error("Category names are not unique")
    model2index = OrderedDict(model_names |> enumerate |> collect .|> reverse)
    category2index =
        OrderedDict(category_names |> enumerate |> collect .|> reverse)

    # Add missing model_name and key in model_names so that each model present in
    # model_names is missing units or has an unit
    for model_name in model_names
        !haskey(units, model_name) && (units[model_name] = "")
    end

    # Delete model_name in units if they do not appear in model_names. We do not want to
    # have unnecessary units in the Dict
    for key in keys(units)
        !(key in model_names) && delete!(units, key)
    end

    # Check number of model to unit pairs match the number of models
    length(units) != length(model_names) && error(
        "The number of unit for each model ($length(units)) is not equal to the number of models ($length(model_names))",
    )
    return RMSEVariable(short_name, model2index, category2index, RMSEs, units)
end

"""
    RMSEVariable(short_name, model_names::Vector{String})

Construct a RMSEVariable with the `short_name` of the variable and the names of the
models in `model_names`.

The categories default to "ANN", "DJF", "MAM", "JJA", "SON". The root mean square errors
default to `NaN`. The unit for each model is missing which is denoted by an empty string.
"""
function RMSEVariable(short_name, model_names::Vector{String})
    category_names = ["ANN", "DJF", "MAM", "JJA", "SON"]
    RMSEs = fill(NaN, length(model_names), length(category_names))
    units = Dict{valtype(model_names), String}([
        (model_name, "") for model_name in model_names
    ])
    return RMSEVariable(short_name, model_names, category_names, RMSEs, units)
end

"""
    RMSEVariable(short_name, model_names::Vector{String}, units::Dict)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, and provided units in the dictionary `units` that map model name to unit.

The categories default to "ANN", "DJF", "MAM", "JJA", "SON". The root mean square errors
default to `NaN`. Any missing model in the dictionary `units` will has missing unit which
is denoted by an empty string.
"""
function RMSEVariable(short_name, model_names::Vector{String}, units::Dict)
    category_names = ["ANN", "DJF", "MAM", "JJA", "SON"]
    RMSEs = fill(NaN, length(model_names), length(category_names))
    return RMSEVariable(short_name, model_names, category_names, RMSEs, units)
end

"""
    RMSEVariable(short_name,
                 model_names::Vector{String},
                 category_names::Vector{String},
                 units::Dict)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, the categories in `category_names`, and provided units in the dictionary
`units` that map model name to unit.

The root mean square errors default to `NaN`. Any missing model in the dictionary `units`
will has missing unit which is denoted by an empty string.
"""
function RMSEVariable(
    short_name,
    model_names::Vector{String},
    category_names::Vector{String},
    units::Dict,
)
    RMSEs = fill(NaN, length(model_names), length(category_names))
    return RMSEVariable(short_name, model_names, category_names, RMSEs, units)
end

"""
    RMSEVariable(short_name::String,
                 model_names::Vector{String},
                 category_names::Vector{String},
                 RMSEs,
                 units::String)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, the categories in `category_names`, the root mean squared errors in `RMSEs`,
and units which map each model name to `units`.

This is useful if all the models share the same unit.
"""
function RMSEVariable(
    short_name::String,
    model_names::Vector{String},
    category_names::Vector{String},
    RMSEs,
    units::String,
)
    units_dict = Dict(model_name => units for model_name in model_names)
    return RMSEVariable(
        short_name,
        model_names,
        category_names,
        RMSEs,
        units_dict,
    )
end

"""
    RMSEVariable(short_name, model_names::Vector{String}, units::String)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, and units which map each model name to `units`.

The categories default to "ANN", "DJF", "MAM", "JJA", "SON". The root mean square errors
default to `NaN`.

This is useful if all the models share the same unit.
"""
function RMSEVariable(short_name, model_names::Vector{String}, units::String)
    units_dict = Dict(model_name => units for model_name in model_names)
    return RMSEVariable(short_name, model_names, units_dict)
end

"""
    RMSEVariable(short_name,
                 model_names::Vector{String},
                 category_names::Vector{String},
                 units::String)

Construct a RMSEVariable with the `short_name` of the variable, the names of the models in
`model_names`, the categories in `category_names`, and units which map each model name to
`units`.

The root mean square errors default to `NaN`.

This is useful if all the models share the same unit.
"""
function RMSEVariable(
    short_name,
    model_names::Vector{String},
    category_names::Vector{String},
    units::String,
)
    RMSEs = fill(NaN, length(model_names), length(category_names))
    return RMSEVariable(short_name, model_names, category_names, RMSEs, units)
end

"""
    model_names(rmse_var::RMSEVariable)

Return all the model names in `rmse_var`.
"""
model_names(rmse_var::RMSEVariable) = rmse_var.model2index |> keys |> collect

"""
    category_names(rmse_var::RMSEVariable)

Return all the category names in `rmse_var`.
"""
category_names(rmse_var::RMSEVariable) =
    rmse_var.category2index |> keys |> collect

"""
    rmse_units(rmse_var::RMSEVariable)

Return all the unit of the models in `rmse_var`.
"""
rmse_units(rmse_var::RMSEVariable) = rmse_var.units

end
