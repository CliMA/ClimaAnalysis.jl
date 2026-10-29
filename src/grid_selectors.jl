"""
    AbstractSelector

An object that determines which indices are selected.

TODO: Write docs for this for grid and dim
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
    TODO
"""
function get_index(dim::Dim, value, ::NearestValue)
    dim_values = coord_values(dim)
    eltype(dim_values) <: AbstractString && error(
        "Cannot find the nearest value in $(name(dim)) dimension because its values are labels. Use MatchValue() or Index() instead",
    )
    return nearest_index(dim_values, value)
end

function get_index(dim::Dim, value, ::MatchValue)
    dim_values = coord_values(dim)
    idx = findfirst(x -> _matches(x, value), dim_values)
    isnothing(idx) && error(
        "Cannot find $value in $(name(dim)) dimension with values $dim_values",
    )
    return idx
end

_matches(x::Number, value::Number) = x ≈ value
_matches(x, value) = isequal(x, value)

function get_index(dim::Dim, idx, ::Index)
    idx ∉ eachindex(coord_values(dim)) &&
        error("Attempt to access $(name(dim)) dimension at index $idx")
    return idx
end

function select(grid::Grid; by::AbstractSelector = NearestValue(), inds...)
    return grid[select_indices(grid; by, inds...)...]
end

function select_indices(
    grid::Grid;
    by::AbstractSelector = NearestValue(),
    inds...,
)
    selected_dim_names =
        [dim(grid, String(dim_name)).name for dim_name in keys(inds)]
    selected_indices = [
        _selector_indices(grid, by, dim_name, indices_or_vals) for
        (dim_name, indices_or_vals) in zip(selected_dim_names, values(inds))
    ]
    return map(grid.dims) do dim
        j = findfirst(==(name(dim)), selected_dim_names)
        isnothing(j) ? Colon() : selected_indices[j]
    end
end

function _selector_indices(
    grid::Grid,
    by::AbstractSelector,
    dim_name,
    indices_or_vals,
)
    # TODO: Fix this later
    # Low priority: dispatch per element, not `first` (fails on empty/CFTime input)
    # We cannot support `CartesianIndex` because of the line below, because it is difficult
    # to tell whether indices_or_vals is a Dates.DateTime or an iterable of Dates.DateTime
    val_is_date =
        indices_or_vals isa Dates.AbstractDateTime ||
        first(indices_or_vals) isa Dates.AbstractDateTime
    dim_is_time = conventional_dim_name(dim_name) == "time"
    (val_is_date && !dim_is_time) &&
        error("Dates are only supported with time dimension")
    val_is_date &&
        (indices_or_vals = _date_to_time.(Ref(grid), indices_or_vals))
    return get_index.(Ref(dim(grid, dim_name)), indices_or_vals, Ref(by))
end

function window_indices(
    grid::Grid,
    dim_name;
    left = nothing,
    right = nothing,
    by::AbstractSelector = NearestValue(),
)
    dim_is_time = conventional_dim_name(dim_name) == "time"

    # If left/right is a Date, let's compute the associated time
    function _maybe_convert_to_time(num)
        if num isa Dates.DateTime
            dim_is_time ||
                error("Dates are only supported with the time dimension")
            return _date_to_time(grid, num)
        else
            return num
        end
    end
    left, right = _maybe_convert_to_time(left), _maybe_convert_to_time(right)

    selected_dim = dim(grid, dim_name)
    dim_values = coord_values(selected_dim)
    left_idx =
        isnothing(left) ? firstindex(dim_values) :
        get_index(selected_dim, left, by)
    right_idx =
        isnothing(right) ? lastindex(dim_values) :
        get_index(selected_dim, right, by)
    right_idx >= left_idx ||
        error("Right window value has to be larger than left one")
    return left_idx:right_idx
end

# TODO: Rewrite this so that it deal with Dim instead of the Grid
# Future case might mean we deal with the Grid itself?
# get_index(dim::Dim, value, ::NearestValue)
# get_index(dim::Dim, value, ::MatchValue)
# get_index(dim::Dim, idx, ::Index)
