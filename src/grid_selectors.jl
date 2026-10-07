"""
    AbstractSelector

An object that determines which indices are selected.

`AbstractSelector`s have to provide one function, `get_index`

The function has to have the signature
`get_index(var, dim_name, idx_or_val, ::AbstractSelector)` and return a single index. You
can assume that `dim_name` is in `ClimaAnalysis.dim_names(var)`.

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

# Selection and structural indexing

function get_index(grid::Grid, dim_name, value, ::NearestValue)
    return nearest_index(coord_values(dim(grid, dim_name)), value)
end

function get_index(grid::Grid, dim_name, value, ::MatchValue)
    dim_values = coord_values(dim(grid, dim_name))
    idx = findfirst(x -> x ≈ value, dim_values)
    isnothing(idx) && error(
        "Cannot find $value in $dim_name dimension with values $dim_values",
    )
    return idx
end

function get_index(grid::Grid, dim_name, idx, ::Index)
    idx ∉ eachindex(coord_values(dim(grid, dim_name))) &&
        error("Attempt to access $dim_name dimension at index $idx")
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
    return get_index.(Ref(grid), Ref(dim_name), indices_or_vals, Ref(by))
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

    dim_values = coord_values(dim(grid, dim_name))
    left_idx =
        isnothing(left) ? firstindex(dim_values) :
        get_index(grid, dim_name, left, by)
    right_idx =
        isnothing(right) ? lastindex(dim_values) :
        get_index(grid, dim_name, right, by)
    right_idx >= left_idx ||
        error("Right window value has to be larger than left one")
    return left_idx:right_idx
end
