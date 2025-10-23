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

            keep_attributes = ("start_date",)

            for attr in keep_attributes
                if haskey(x.attributes, attr) &&
                   haskey(y.attributes, attr) &&
                   x.attributes[attr] == y.attributes[attr]
                    ret_attributes[attr] = x.attributes[attr]
                end
            end

            # This could add units when x and y do not have units in their attributes.
            # This is okay, because the units will be the empty string which is considered
            # missing by ClimaAnalysis.
            units = find_units_for_binary($op, x, y)
            ret_attributes["units"] = units

            ret_dims = x.dims
            ret_dim_attributes = x.dim_attributes

            ret_data = @. $op(x.data, y.data)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
        function Base.$op(x::OutputVar, y::Real)
            ret_attributes = empty(x.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(x.attributes, attr)
                    ret_attributes[attr] =
                        string(x.attributes[attr], " ", string($op), " ", y)
                end
            end

            keep_attributes = ("start_date", "units")

            for attr in keep_attributes
                if haskey(x.attributes, attr)
                    ret_attributes[attr] = x.attributes[attr]
                end
            end

            ret_dims = deepcopy(x.dims)
            ret_dim_attributes = deepcopy(x.dim_attributes)

            ret_data = @. $op(x.data, y)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
        function Base.$op(x::Real, y::OutputVar)
            ret_attributes = empty(y.attributes)

            specific_attributes = ("short_name", "long_name")

            for attr in specific_attributes
                if haskey(y.attributes, attr)
                    ret_attributes[attr] =
                        string(x, " ", string($op), " ", y.attributes[attr])
                end
            end

            keep_attributes = ("start_date", "units")

            for attr in keep_attributes
                if haskey(y.attributes, attr)
                    ret_attributes[attr] = y.attributes[attr]
                end
            end

            ret_dims = deepcopy(y.dims)
            ret_dim_attributes = deepcopy(y.dim_attributes)

            ret_data = @. $op(x, y.data)

            return OutputVar(
                ret_attributes,
                ret_dims,
                ret_dim_attributes,
                ret_data,
            )
        end
    end
end

"""
    find_units_for_binary(
        f::Function,
        x_var::OutputVar,
        y_var::OutputVar,
    )

Given a binary operator `f` and the `OutputVar`s `x` and `y`, find the appropriate units to
return.
"""
function find_units_for_binary(
    ::Union{typeof(+), typeof(-), typeof(max), typeof(min)},
    x_var::OutputVar,
    y_var::OutputVar,
)
    x_units = units(x_var)
    y_units = units(y_var)
    return x_units == y_units ? x_units : ""
end

function find_units_for_binary(
    f::Union{typeof(*), typeof(/)},
    x_var::OutputVar,
    y_var::OutputVar,
)
    x_units = units(x_var)
    y_units = units(y_var)
    # Empty string is considered missing units
    if x_units == "" || y_units == ""
        return ""
    end
    return "($x_units) $f ($y_units)"
end

@overload_binary_op (+)
@overload_binary_op (-)
@overload_binary_op (*)
@overload_binary_op (/)
@overload_binary_op max
@overload_binary_op min

"""
    @overload_unary_op op

Generate a method to overload the unary operator `op` for `OutputVar`.

Handles attributes `short_name`, `long_name`: Prepends the operator name, e.g., "log(Temperature)".
"""
macro overload_unary_op(op)
    esc(
        quote
            function Base.$op(x::OutputVar)
                # Start by copying all attributes from x
                ret_attributes = Dict{String, Any}()

                # Modify specific attributes like names
                specific_attributes = ("short_name", "long_name")
                for attr in specific_attributes
                    if haskey(x.attributes, attr)
                        # Prepend the operator name
                        ret_attributes[attr] =
                            string($(string(op)), "(", x.attributes[attr], ")")
                    end
                end

                keep_attributes = ("start_date",)
                for attr in keep_attributes
                    if haskey(x.attributes, attr)
                        ret_attributes[attr] = x.attributes[attr]
                    end
                end

                # Keep dimensions and their attributes (use deepcopy for safety)
                ret_dims = deepcopy(x.dims)
                ret_dim_attributes = deepcopy(x.dim_attributes)

                # Perform operation element-wise
                ret_data = @. $op(x.data)

                return OutputVar(
                    ret_attributes,
                    ret_dims,
                    ret_dim_attributes,
                    ret_data,
                )
            end
        end,
    )
end

@overload_unary_op log
@overload_unary_op log2
@overload_unary_op log10
@overload_unary_op exp
@overload_unary_op sin
@overload_unary_op cos
@overload_unary_op tan
@overload_unary_op sqrt
@overload_unary_op (-) # Unary minus
