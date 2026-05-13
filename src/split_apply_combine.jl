export AbstractSplitOperation,
    AbstractApplyOperation,
    SplitApplyVar,
    GroupAll,
    SplitSeason,
    Reduce,
    combine

"""
    abstract type AbstractSplitOperation

Supertype for operations that partition an `OutputVar` into groups in the
split-apply-combine pattern.
"""
abstract type AbstractSplitOperation end

"""
    abstract type AbstractApplyOperation

Supertype for operations applied to each group in a split-apply-combine pattern.
"""
abstract type AbstractApplyOperation end

"""
    GroupAll <: AbstractSplitOperation

Specify the dimension over which to apply an operation in a split-apply-combine pattern.

`GroupAll` treats the entire dimension as a single group (i.e. no splitting occurs). It is
used in combination with a [`Reduce`](@ref) (or other apply operation) and [`combine`](@ref)
to transform an `OutputVar`.
"""
struct GroupAll <: AbstractSplitOperation
    "Dimension to treat as a single group"
    dim_name::String
end

"""
    SplitSeason <: AbstractSplitOperation end

Split the dates of the `OutputVar` into seasons, in chronological order, as part of the
split-apply-combine pattern.

If a season is empty, it will not be present as a group. The dates do not necessarily need
to be in sorted order.

# Example

This is a singleton type and can be constructed as follows.

```julia
ClimaAnalysis.SplitSeason()
```
"""
struct SplitSeason <: AbstractSplitOperation end

"""
    Reduce{F} <: AbstractApplyOperation

Specify a reduction function to apply to each group in a split-apply-combine pattern.

The reduction `f` must accept an `Array` and a `dims` keyword argument, and return an array
of the same number of dimensions where the size along `dims` is 1 (i.e., `f` does not drop
the reduced dimension).

It is used in combination with an [`GroupAll`](@ref) and [`combine`](@ref) to transform an
`OutputVar`.

# Example

You can pass the reduction function when constructing a `Reduce`.

```julia
ClimaAnalysis.Reduce(sum)
```
"""
struct Reduce{F} <: AbstractApplyOperation
    "The reduction function applied to each group along the split dimension"
    reduction::F
end

"""
    SplitApplyVar

Represents lazy split and apply operations before combining to produce the resulting
`OutputVar`.
"""
struct SplitApplyVar{
    OUTPUTVAR <: OutputVar,
    SPLIT <: AbstractSplitOperation,
    APPLY <: Union{AbstractApplyOperation, Nothing},
}
    "The `OutputVar` to split, apply, and combine."
    var::OUTPUTVAR

    "A split operation that partitions the dimension into groups of indices."
    split_op::SPLIT

    "An apply operation to apply to each group, or `nothing` if not yet set."
    apply_op::APPLY
end

"""
    (split_op::AbstractSplitOperation)(var::OutputVar)

Lazily applies `split_op` to `var` as part of the split-apply-combine pattern.

The result should be passed to an [`AbstractApplyOperation`](@ref) and then to
[`combine`](@ref) to produce the final `OutputVar`.
"""
function (split_op::AbstractSplitOperation)(var::OutputVar)
    return SplitApplyVar(var, split_op, nothing)
end

"""
    (split_op::AbstractSplitOperation)(::SplitApplyVar)

This throws a helpful error message that you cannot apply another split
operation.
"""
function (split_op::AbstractSplitOperation)(::SplitApplyVar)
    error("A split operation is already set")
end

"""
    (apply_op::AbstractApplyOperation)(split_apply_var::SplitApplyVar)

Lazily applies `apply_op` to `split_apply_var` as part of the split-apply-combine pattern.

The result should be passed to [`combine`](@ref) to produce the final `OutputVar`.

You cannot applies another `AbstractApplyOperation` to a `split_apply_var` that already have
an `AbstractApplyOperation` applied on it.
"""
function (apply_op::AbstractApplyOperation)(split_apply_var::SplitApplyVar)
    curr_apply_op = split_apply_var.apply_op
    isnothing(curr_apply_op) || error("An apply operation is already set")
    (; var, split_op) = split_apply_var
    return SplitApplyVar(var, split_op, apply_op)
end

"""
    combine(split_apply_var::SplitApplyVar)

Materialize the `OutputVar` with the split and apply operations.

This splits the `OutputVar` into groups via an `AbstractSplitOperation`, applies
an `AbstractApplyOperation` to each group, and combine the results together. The
split dimension is kept. The coordinate values of the split dimension in the
returned `OutputVar` are the first coordinate value from each group.
"""
function combine(split_apply_var::SplitApplyVar)
    (; var, split_op, apply_op) = split_apply_var
    isnothing(apply_op) && error("No apply operation is defined")
    return _combine(var, split_op, apply_op)
end

"""
    _combine(
        var::OutputVar,
        split_op::AbstractSplitOperation,
        apply_op::AbstractApplyOperation,
    )

Create groups using `split_op`, perform an operation on each group using `apply_op`, and
combine the results together.
"""
function _combine(
    var::OutputVar,
    split_op::AbstractSplitOperation,
    apply_op::AbstractApplyOperation,
)
    dim_idx, groups = _create_groups(var, split_op)
    return _apply_and_combine(var, apply_op, dim_idx, groups)
end

"""
    _create_groups(var::OutputVar, group_all::GroupAll)

Return the index of the grouped dimension and a single-element vector consisting of a vector
of all indices of that dimension.
"""
function _create_groups(var::OutputVar, group_all::GroupAll)
    (; dim_name) = group_all
    dim_name_in_var = find_corresponding_dim_name_in_var(dim_name, var)
    dim_idx = var.dim2index[dim_name_in_var]
    return dim_idx, [eachindex(var.dims[dim_name_in_var])]
end

"""
    _create_groups(var::OutputVar, ::SplitSeason)

Return the index of the time dimension and a vector of index vectors, one per non-empty
season.
"""
function _create_groups(var::OutputVar, ::SplitSeason)
    _check_time_dim(var)

    # Get dimension index
    dim_name_in_var = find_corresponding_dim_name_in_var("time", var)
    dim_idx = var.dim2index[dim_name_in_var]

    # Compute groups
    var_dates = dates(var)
    date_groups = ClimaAnalysis.split_by_season_across_time(var_dates)
    dim_indices_groups =
        [indexin(group, var_dates) for group in date_groups if !isempty(group)]
    return dim_idx, dim_indices_groups
end

"""
    _apply_and_combine(
        var::OutputVar,
        apply_op::Reduce,
        dim_idx,
        dim_indices_groups,
    )

Apply `apply_op` over the `dim_idx`th dimension for each group specified by
`dim_indices_groups` over the data of `var` and combine the results together.
"""
function _apply_and_combine(
    var::OutputVar,
    apply_op::Reduce,
    dim_idx,
    dim_indices_groups,
)
    data_groups = map(dim_indices_groups) do dim_indices
        index_tuple = ntuple(
            idx -> idx == dim_idx ? dim_indices : Colon(),
            ndims(var.data),
        )
        view(var.data, index_tuple...)
    end

    reduction = apply_op.reduction
    reduced = map(g -> reduction(g; dims = dim_idx), data_groups)
    ret_data = cat(reduced..., dims = dim_idx)

    # Note: For Reduce, the size of the resulting dimension is the same as the
    # number of groups
    # We will always get the first element of the dimension
    ret_dim_indices = (first(dim_indices) for dim_indices in dim_indices_groups)

    # New dimension to return
    dim_name_in_var = var.index2dim[dim_idx]
    dim = [var.dims[dim_name_in_var][i] for i in ret_dim_indices]

    # Make new dimensions for OutputVar
    ret_dims = deepcopy(var.dims)
    ret_dims[dim_name_in_var] = dim

    # Do not update the attributes or dimension attributes
    return ClimaAnalysis.remake(var, dims = ret_dims, data = ret_data)
end
