export slice,
    window, NearestValue, MatchValue, Index, AbstractSelector, get_index

"""
    AbstractSelector

An object that determines which indices are selected.

`AbstractSelector`s have to provide one function, `get_index`

The function has to have the signature
`get_index(var, dim_name, idx_or_val, ::AbstractSelector)` and return a single index. You
can assume that `dim_name` is in `keys(var.dims)`.

The function `get_index` is used by [`slice`](@ref) and [`window`](@ref). For instance, if
you use `ClimaAnalysis.slice(var, time = 2)`, then `dim_name` is `time`, and `idx_or_val` is
`2`.

For example, to implement an `AbstractSelector` that always take the first index of the
dimension regardless of the value or index passed in, you can write the following:

```julia
struct FirstIndex <: ClimaAnalysis.AbstractSelector end

function ClimaAnalysis.get_index(var, dim_name, idx_or_val, ::FirstIndex)
    return 1
end

# Get the first time slice of var. The parameter 10 does not do anything.
ClimaAnalysis.slice(var, time = 10, by = FirstIndex())
```
"""
abstract type AbstractSelector end

"""
    NearestValue

Get the index of the nearest value.

If the dimension is not one dimesional, then an error is thrown.
"""
struct NearestValue <: AbstractSelector end

"""
    MatchValue

Get the index of the approximately matched value.

If the value does not exist, or the dimension is not one-dimensional, then an error is
thrown.
"""
struct MatchValue <: AbstractSelector end

"""
    Index

Select the index in the dimension.

If the index is out of bounds, then an error is thrown.
"""
struct Index <: AbstractSelector end

"""
    _slice_general(var, val, dim_name, by)

Return a new OutputVar by selecting the available index closest to the given `val` for the
given dimension.

When a time dimension is selected, `start_date` is among the `var`'s attributes, `val` can
also be a `Dates.DateTime`. In this case, the time associated to the given date is
automatically computed.

The `by` keyword argument is an `AbstractSelector`. See the `slice` function for more
information.
"""
function _slice_general(var, val, dim_name, by)
    dim_name_in_var = find_corresponding_dim_name_in_var(dim_name, var)

    val_is_date = val isa Dates.DateTime
    dim_is_time = conventional_dim_name(dim_name_in_var) == "time"
    (val_is_date && !dim_is_time) &&
        error("Dates are only supported with time dimension")
    # If val is a Date, let's compute the associated time
    val_is_date && (val = date_to_time(var, val))

    _slice_over = _create_slice_over(var, dim_name_in_var, val, by)
    idx = get_index(var, dim_name_in_var, val, by)
    reduced_var = _reduce_over(_slice_over, dim_name_in_var, var)

    # Let's try adding this operation to the long_name, if possible (ie, if the correct
    # attributes are available)
    if haskey(var.attributes, "long_name") &&
       haskey(var.dim_attributes, dim_name_in_var) &&
       haskey(var.dim_attributes[dim_name_in_var], "units")
        dim_array = var.dims[dim_name_in_var]
        dim_units = var.dim_attributes[dim_name_in_var]["units"]
        cut_point = dim_array[idx]

        dim_units_are_seconds = dim_units == "s"

        if dim_is_time && dim_units_are_seconds
            # Dimension is time and units are seconds. Let's convert them to something nicer
            pretty_timestr = seconds_to_prettystr(cut_point)
            reduced_var.attributes["long_name"] *= " $dim_name_in_var = $pretty_timestr"
        else
            reduced_var.attributes["long_name"] *= " $dim_name_in_var = $cut_point $dim_units"
        end
        reduced_var.attributes["slice_$dim_name_in_var"] = "$cut_point"
        reduced_var.attributes["slice_$(dim_name_in_var)_units"] = dim_units
    end
    return reduced_var
end

"""
    slice(var::OutputVar, by = NearestValue(), kwargs...)

Return a new `OutputVar` by slicing across dimensions as defined by the keyword arguments.

When a time dimension is selected and `start_date` is among `var`'s attributes, it is
possible to directly pass a `Dates.DateTime` and `slice` will automatically convert it into
the corresponding time.

The keyword argument `by` can be `ClimaAnalysis.NearestValue()`,
`ClimaAnalysis.MatchValue()`, or `ClimaAnalysis.Index()`.
- `NearestValue()`: The nearest value will be selected for slicing.
- `MatchValue()`: The approximately matched value will be selected for slicing.
- `Index()`: The index passed in will be used for slicing.

Example
===========
```julia
slice(var, lat = 30, lon = 20, time = 100)
```

If `var` has `start_date` among its attributes.
```julia
import Dates
slice(var, lat = 30, lon = 20, time = Dates.DateTime(2020, 12, 21))
```

!!! compat "Support for Dates and generalized dimension names"
    Calling `slice` with a `DateTime` and specifying a dimension by one of its name
    (instead of the actual name in the file) was introduced in ClimaAnalysis v0.5.17.

!!! compat "Support for AbstractSelector"
    Calling `slice` with the `by` keyword argument was introduced after ClimaAnalysis
    v0.5.18.

"""
function slice(var; by = NearestValue(), kwargs...)
    sliced_var = var
    for (dim_name, val_or_idx) in kwargs
        sliced_var =
            _slice_general(sliced_var, val_or_idx, String(dim_name), by)
    end
    return sliced_var
end

"""
    window(var::OutputVar, dim_name; left = nothing, right = nothing, by = NearestValue())

Return a new OutputVar by selecting the values of the given `dim`ension that are between
`left` and `right`.

If `left` and/or `right` are `nothing`, assume beginning (or end) of the array.

For the time dimension, `left` and `right` can also be `Dates.DateTime` (if `var` contains a
`start_date`).

The keyword argument `by` can be `ClimaAnalysis.NearestValue()`,
`ClimaAnalysis.MatchValue()`, or `ClimaAnalysis.Index()`.
- `NearestValue()`: The nearest value will be selected for windowing.
- `MatchValue()`: The approximately matched value will be selected for windowing.
- `Index()`: The index passed in will be used for windowing.

Example
===========
```julia
window(var, 'lat', left = -50, right = 50)
window(var, 'time', left = DateTime(2008), right = DateTime(2009))
```

!!! compat "Support for Dates and generalized dimension names"
    Calling `window` with a `DateTime` and specifying a dimension by one of its name
    (instead of the actual name in the file) was introduced in ClimaAnalysis v0.5.17.

!!! compat "Support for AbstractSelector"
    Calling `window` with the `by` keyword argument was introduced after ClimaAnalysis
    v0.5.18.

"""
function window(
    var,
    dim_name;
    left = nothing,
    right = nothing,
    by = NearestValue(),
)
    dim_name_in_var = find_corresponding_dim_name_in_var(dim_name, var)

    dim_is_time = conventional_dim_name(dim_name_in_var) == "time"

    # If left/right is a Date, let's compute the associated time
    function _maybe_convert_to_time(num)
        if num isa Dates.DateTime
            dim_is_time ||
                error("Dates are only supported with the time dimension")
            return date_to_time(var, num)
        else
            return num
        end
    end
    left, right = _maybe_convert_to_time(left), _maybe_convert_to_time(right)
    nearest_index_left, nearest_index_right =
        _get_window_indices(var, dim_name_in_var, left, right, by)

    (nearest_index_right >= nearest_index_left) ||
        error("Right window value has to be larger than left one")

    # Take only what's between nearest_index_left and nearest_index_right
    reduced_data = selectdim(
        var.data,
        var.dim2index[dim_name_in_var],
        nearest_index_left:nearest_index_right,
    )

    dims = copy(var.dims)
    reduced_dim =
        var.dims[dim_name_in_var][nearest_index_left:nearest_index_right]
    dims[dim_name_in_var] = reduced_dim

    dim_attributes = copy(var.dim_attributes)
    return OutputVar(copy(var.attributes), dims, dim_attributes, reduced_data)
end

"""
    get_index(var, dim_name, val, ::NearestValue)

Get the index corresponding to the nearest value to `val` in the `dim_name` dimension.

This function assumes that `dim_name` exists in `var`, and the `dim_name` dimension is
one-dimensional.
"""
function get_index(var, dim_name, val, ::NearestValue)
    dim = var.dims[dim_name]
    ndims(dim) > 1 && error(
        "$dim_name dimension is not one-dimensional. Use Index() instead of NearestValue()",
    )
    return nearest_index(dim, val)
end

"""
    get_index(var, dim_name, val, ::MatchValue)

Get the index that approximately match with `val` in the `dim_name` dimension.

If a `val` is not in the `dim_name` dimension, then an error is thrown.

This function assumes that `dim_name` exists in `var`, and the `dim_name` dimension is
one-dimensional.
"""
function get_index(var, dim_name, val, ::MatchValue)
    dim = var.dims[dim_name]
    ndims(dim) > 1 && error(
        "$dim_name dimension is not one-dimensional. Use Index() instead of MatchValue()",
    )
    idx = findfirst(x -> x ≈ val, dim)
    isnothing(idx) &&
        error("Cannot find $val in $dim_name dimension with values $dim")
    return idx
end

"""
    get_index(var, dim_name, idx, ::Index)

Get the index `idx`.

If `idx` is out of bounds, then an error is thrown.

This function assumes that `dim_name` exists in `var`.
"""
function get_index(var, dim_name, idx, ::Index)
    dim_idx = var.dim2index[dim_name]
    idx ∉ axes(var.data, dim_idx) &&
        error("Attempt to access $dim_name dimension at index $idx")
    return idx
end

"""
    _get_window_indices(var, dim_name, left, right, by::AbstractSelector)

Get the indices for the left and right values for `window`.

If `nothing` is provided for `left` or `right`, then the first or last point are used
respectively.
"""
function _get_window_indices(var, dim_name, left, right, by::AbstractSelector)
    left_idx =
        isnothing(left) ? firstindex(var.dims[dim_name]) :
        get_index(var, dim_name, left, by)
    right_idx =
        isnothing(right) ? lastindex(var.dims[dim_name]) :
        get_index(var, dim_name, right, by)
    return left_idx, right_idx
end

"""
    _create_slice_over(var, dim_name, idx_or_val, by::AbstractSelector)

Return a `slice_over` function which is used by slicing functions.
"""
function _create_slice_over(var, dim_name, idx_or_val, by::AbstractSelector)
    idx = get_index(var, dim_name, idx_or_val, by)
    function _slice_over(data; dims)
        if !(dims isa Integer) && !(length(dims) == 1)
            throw(
                ArgumentError(
                    "selectdim_as_reduction only supports a single dimension",
                ),
            )
        end

        dim = isa(dims, Integer) ? dims : only(dims)
        sliced = selectdim(data, dim, idx)

        # Insert singleton dimension back at the same position, so that it behaves like a reduction
        new_shape = size(sliced)
        new_shape_with_dim = ntuple(
            i -> (i < dim ? new_shape[i] : i == dim ? 1 : new_shape[i - 1]),
            ndims(data),
        )

        return reshape(sliced, new_shape_with_dim)
    end
    return _slice_over
end
