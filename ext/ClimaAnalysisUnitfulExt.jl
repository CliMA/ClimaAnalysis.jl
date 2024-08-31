module ClimaAnalysisUnitfulExt

# ClimaAnalysisUnitfulExt is implemented as extension in case we decide to turn Unitful into
# an extension in the future, but, right now, everything is included directly in
# ClimaAnalysis.

import Unitful

import ClimaAnalysis: Var

"""
    _maybe_convert_to_unitful(value)

Try converting `value` to a `Uniftul` object. If unsuccessful, just return it.
"""
function Var._maybe_convert_to_unitful(value)
    value isa Unitful.Units && return value
    # This function in inherently type-unstable
    try
        return Unitful.uparse(value)
    catch exc
        # ParseError when it cannot be parsed
        # ArgumentError when symbols are not available
        if exc isa Base.Meta.ParseError || exc isa ArgumentError
            return value
        else
            rethrow(exc)
        end
    end
end

function _converted_data(data, conversion_function)
    return conversion_function.(data)
end

function _converted_data_unitful(data, old_units, new_units)
    # We add FT because sometimes convert changes the type and because
    # ustrip only reinterprets the given array
    FT = eltype(data)
    return FT.(Unitful.ustrip(Unitful.uconvert.(new_units, data * old_units)))
end

function Var.convert_units(
    var::Var.OutputVar,
    new_units::AbstractString;
    conversion_function = nothing,
)
    has_unitful_units =
        Var.has_units(var) && (var.attributes["units"] isa Unitful.Units)
    new_units_maybe_unitful = Var._maybe_convert_to_unitful(new_units)
    new_units_are_unitful = new_units_maybe_unitful isa Unitful.Units

    if has_unitful_units && new_units_are_unitful
        isnothing(conversion_function) ||
            @warn "Ignoring conversion_function, units are parseable."
        convert_function =
            data -> _converted_data_unitful(
                data,
                var.attributes["units"],
                new_units_maybe_unitful,
            )
    else
        isnothing(conversion_function) && error(
            "Conversion function required for var with non-parseable/absent units.",
        )
        convert_function = data -> _converted_data(data, conversion_function)
    end

    new_data = convert_function(var.data)
    new_attribs = copy(var.attributes)
    # The constructor will take care of converting new_units to Unitful
    new_attribs["units"] = new_units

    return Var.OutputVar(new_attribs, var.dims, var.dim_attributes, new_data)
end

end
