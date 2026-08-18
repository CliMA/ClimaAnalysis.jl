export AbstractSplitOperation,
    AbstractApplyOperation,
    SplitApplyVar,
    Group,
    GroupBy,
    GroupAll,
    SplitSeason,
    Bins,
    Reduce,
    WithCoords,
    ReduceWithCoords,
    Average,
    LatWeightedAverage,
    combine

"""
    abstract type AbstractSplitOperation

Supertype for operations that partition an `OutputVar` into groups in the
split-apply-combine pattern.

A subtype must be a struct with a `dim_name` field and implement
`_group_indices(var, op)`, returning
`(; idx, name, groups, labels, coords)`: the index and name of the grouped
dimension, the vector of index groups, the label of each group, and the
representative coordinate of each group.
"""
abstract type AbstractSplitOperation end

"""
    abstract type AbstractApplyOperation

Supertype for operations applied to each block in a split-apply-combine
pattern.
"""
abstract type AbstractApplyOperation end

"""
    GroupBy(dim_name, by; on = nothing)

Group the dimension `dim_name` by the labels `by.(values)` as part of the
split-apply-combine pattern.

Indices whose labels are `isequal` form a group and an index whose label is
`nothing` is dropped. By default, the values passed to `by` are the
coordinates of `dim_name`, except for a time dimension of an `OutputVar` with
a `start_date` attribute, whose coordinates are converted to
`Dates.DateTime`s first (the units of the time dimension must be seconds).
Pass `on`, a function `var -> vector`, to label the dimension by other values
instead.

Groups are ordered by their smallest coordinate (largest for a descending
dimension), so a monotone dimension stays monotone, and that coordinate
represents the group in the combined `OutputVar`.

Chaining another `GroupBy` groups several dimensions at once and the apply
operation runs once per block of the Cartesian product of the groups.
Grouping the same dimension twice is an error.

# Example

```julia
import Dates

var |> ClimaAnalysis.GroupBy("time", Dates.year)   # annual groups
var |> ClimaAnalysis.GroupBy("time", Dates.month)  # monthly climatology
var |> ClimaAnalysis.GroupBy("lat", l -> abs(l) < 30.0)  # tropics or not

# Group raw seconds instead of dates
var |> ClimaAnalysis.GroupBy(
    "time",
    t -> fld(t, 86_400.0),
    on = var -> var.dims["time"],
)
```
"""
struct GroupBy{F, C} <: AbstractSplitOperation
    "Dimension to group along"
    dim_name::String

    "Function mapping a value to a label, or `nothing` to drop the index"
    by::F

    "Function `var -> vector` of values to label; `nothing` uses coordinates"
    on::C
end

GroupBy(dim_name, by; on = nothing) = GroupBy(dim_name, by, on)

"""
    GroupAll(dim_name)

Treat the entire dimension `dim_name` as a single group (i.e. no splitting
occurs) as part of the split-apply-combine pattern.

It is used in combination with a [`Reduce`](@ref) (or other apply operation)
and [`combine`](@ref) to transform an `OutputVar`.

# Example

```julia
import Statistics: mean

var |>
ClimaAnalysis.GroupAll("time") |>
ClimaAnalysis.Reduce(mean) |>
ClimaAnalysis.combine
```
"""
GroupAll(dim_name) = GroupBy(
    dim_name,
    Returns(:all),
    var -> var.dims[find_corresponding_dim_name_in_var(dim_name, var)],
)

"""
    SplitSeason()

Split the dates of the `OutputVar` into seasons (DJF, MAM, JJA, and SON), in
chronological order for each year, as part of the split-apply-combine pattern.

The `OutputVar` must have a time dimension in units of seconds and a
`start_date` attribute. If a season is empty, it will not be present as a
group. The dates do not need to be in sorted order.

# Example

```julia
ClimaAnalysis.SplitSeason()
```
"""
SplitSeason() = GroupBy("time", find_season_and_year)

"""
    Bins(edges)

A labeler for [`GroupBy`](@ref) mapping a value to the index of its bin, or to
`nothing` (dropping the value) when it falls outside every bin.

Bins are `[edges[i], edges[i + 1])`, except the last bin, which also contains
`last(edges)`. `NaN` is never in a bin. The edges must be in increasing order.

# Example

```julia
var |> ClimaAnalysis.GroupBy("lat", ClimaAnalysis.Bins(-90.0:30.0:90.0))
```
"""
struct Bins{E <: AbstractVector}
    "Bin edges; bins are [edges[i], edges[i+1]), last bin closed on the right"
    edges::E

    function Bins(edges::E) where {E <: AbstractVector}
        length(edges) >= 2 || error("At least two bin edges are needed")
        issorted(edges) || error("Bin edges must be in increasing order")
        return new{E}(edges)
    end
end

function (b::Bins)(x)
    i = searchsortedlast(b.edges, x)
    i == length(b.edges) && return x == last(b.edges) ? i - 1 : nothing
    return i == 0 ? nothing : i
end

"""
    _labeling_values(var::OutputVar, dim_name)

Return the values labeled by a [`GroupBy`](@ref) without `on`: the coordinates
of `dim_name`, converted to `Dates.DateTime`s for a time dimension of a `var`
with a `start_date` attribute.
"""
function _labeling_values(var::OutputVar, dim_name)
    dated_time =
        conventional_dim_name(dim_name) == "time" &&
        haskey(var.attributes, "start_date")
    dated_time || return var.dims[dim_name]
    _check_time_dim(var)
    # Not dates(var): that returns the separate date dimension when one exists
    return time_to_date.(
        Dates.DateTime(var.attributes["start_date"]),
        var.dims[dim_name],
    )
end

"""
    _group_indices(var::OutputVar, op::GroupBy)

Return the grouping of the dimension `op.dim_name` of `var` as a `NamedTuple`
`(; idx, name, groups, labels, coords)`.

Groups are ordered by, and represented by, their extreme coordinate (minimum
for an ascending dimension, maximum for a descending one), so results do not
depend on the ordering of the input and a monotone dimension stays monotone.
"""
function _group_indices(var::OutputVar, op::GroupBy)
    # Resolved at pipe time
    name = op.dim_name
    coord = var.dims[name]
    ndims(coord) == 1 ||
        error("Cannot group along the multidimensional dimension $name")
    labeled = isnothing(op.on) ? _labeling_values(var, name) : op.on(var)
    length(labeled) == length(coord) || error(
        "Labeling vector for $name has length $(length(labeled)), " *
        "expected $(length(coord))",
    )

    labels_per_index = map(op.by, labeled)
    buckets = OrderedDict{eltype(labels_per_index), Vector{Int}}()
    for (i, label) in pairs(labels_per_index)
        isnothing(label) || push!(get!(() -> Int[], buckets, label), i)
    end
    isempty(buckets) && error("Grouping $name produced no groups")

    rev = issorted(coord, rev = true)
    pick = rev ? maximum : minimum
    labels = collect(keys(buckets))
    reps = [pick(view(coord, buckets[k])) for k in labels]
    # MergeSort keeps ties deterministic
    perm = sortperm(reps, rev = rev, alg = MergeSort)
    labels, reps = labels[perm], reps[perm]
    groups = _as_ranges([buckets[k] for k in labels])
    return (; idx = var.dim2index[name], name, groups, labels, coords = reps)
end

"""
    _as_ranges(groups)

Convert every group to a `UnitRange` when all groups are consecutive.

Converting all groups or none keeps the element type of `groups` uniform,
which keeps `_combine_blocks` type stable, and a `UnitRange` view is affine
where a `Vector{Int}` view is a gather.
"""
function _as_ranges(groups)
    all(g -> g == first(g):last(g), groups) || return groups
    return [first(g):last(g) for g in groups]
end

"""
    Reduce{F} <: AbstractApplyOperation

Specify a reduction function to apply to each block in a split-apply-combine
pattern.

The reduction `f` is called as `f(block; dims)`, with `dims` the `Tuple` of
the indices of the grouped dimensions, and must return an array of the same
number of dimensions where the size along each grouped dimension is 1 (i.e.,
`f` does not drop the reduced dimensions). Wrap `f` in [`WithCoords`](@ref)
for it to be called as `f(block, group)` instead, where `group` is a
[`Group`](@ref) carrying the coordinates of the block.

It is used in combination with an [`AbstractSplitOperation`](@ref) and
[`combine`](@ref) to transform an `OutputVar`.

# Example

```julia
ClimaAnalysis.Reduce(sum)
```
"""
struct Reduce{F} <: AbstractApplyOperation
    "The reduction function applied to each block"
    reduction::F
end

"""
    WithCoords(f)

Mark `f` as coordinate aware: an apply operation calls it as
`f(block, group)`, with `group` a [`Group`](@ref), instead of
`f(block; dims)`.

# Example

```julia
import NaNStatistics: nanmean

ClimaAnalysis.Reduce(
    ClimaAnalysis.WithCoords(
        (block, group) ->
            nanmean(block, cosd.(group.lat); dims = group.dim_indices),
    ),
)
```
"""
struct WithCoords{F}
    "Function called as `f(block, group)`"
    f::F
end

_call(f, block, group) = f(block; dims = group.dim_indices)
_call(f::WithCoords, block, group) = f.f(block, group)

_apply(op::Reduce, block, group) = _call(op.reduction, block, group)

"""
    ReduceWithCoords(f)

Specify a coordinate-aware reduction function to apply to each block in a
split-apply-combine pattern; shorthand for `Reduce(WithCoords(f))`.

The reduction is called as `f(block, group)`, with `group` a [`Group`](@ref),
and must return an array of the same number of dimensions as `block` where
the size along each grouped dimension is 1.

Note that the coordinates in `group` keep their own element type, so weighting
`Float32` data by `cosd.(group.lat)` with `Float64` latitudes promotes the
result to `Float64`. Use [`Average`](@ref) with `weights` to match the data's
precision instead.

# Example

```julia
import NaNStatistics: nanmean

ClimaAnalysis.ReduceWithCoords(
    (block, group) ->
        nanmean(block, cosd.(group.lat); dims = group.dim_indices),
)
```
"""
ReduceWithCoords(f) = Reduce(WithCoords(f))

"""
    Group

The grouped dimensions of a block, passed to every [`WithCoords`](@ref)
function.

`group.dim_indices` is the `Tuple` of the indices of the grouped dimensions.
The coordinates of each grouped dimension in the block are available by name
(e.g. `group.lat`), reshaped to broadcast against the block; `vec` recovers
the plain vector. Names are conventional, so `group.lat` and `group.latitude`
are the same coordinates whatever the dimension is called in the `OutputVar`.
Asking for a dimension that is not grouped is an error.

`group.coords` is the `NamedTuple` of all grouped coordinates, keyed by
conventional dimension name and ordered like `group.dim_indices`.
"""
struct Group{D <: Tuple, C <: NamedTuple}
    "Indices of the dimensions being reduced"
    dim_indices::D

    "Coordinates of the grouped dimensions, reshaped to broadcast against the block"
    coords::C
end

# Map every recognized dimension name to its conventional name (e.g. :lat and
# :latitude to :latitude); a const NamedTuple rather than a Dict so lookups
# with a literal symbol fold away
const _CONVENTIONAL = NamedTuple(
    Symbol(name) => Symbol(conventional_dim_name(name)) for
    name in Iterators.flatten((
        LONGITUDE_NAMES,
        LATITUDE_NAMES,
        TIME_NAMES,
        DATE_NAMES,
        ALTITUDE_NAMES,
        PRESSURE_NAMES,
    ))
)

function _conventional(name::Symbol)
    hasfield(typeof(_CONVENTIONAL), name) || return name
    return getfield(_CONVENTIONAL, name)
end

function Base.getproperty(group::Group, name::Symbol)
    name in fieldnames(Group) && return getfield(group, name)
    coords = getfield(group, :coords)
    key = _conventional(name)
    haskey(coords, key) || error(
        "$name is not a grouped dimension " *
        "(grouped: $(join(keys(coords), ", ")))",
    )
    return coords[key]
end

Base.propertynames(group::Group) =
    (fieldnames(Group)..., keys(getfield(group, :coords))...)

"""
    Average(; ignore_nan = true, weights = nothing)

Specify a mean, optionally weighted, to apply to each block in a
split-apply-combine pattern.

If `ignore_nan` is `true`, then `NaN`s in the data are excluded from the mean.
If `weights` is not `nothing`, then it is a function of a [`Group`](@ref)
returning an array broadcastable against the block, and the mean is weighted
by it. The weights are converted to the data's floating-point precision, so
averaging `Float32` data returns `Float32`. Weighting a dimension that is not
grouped is an error.

# Example

```julia
ClimaAnalysis.Average()
ClimaAnalysis.Average(weights = group -> cosd.(group.lat) .^ 2)
```
"""
function Average(; ignore_nan = true, weights = nothing)
    if isnothing(weights)
        return ignore_nan ? Reduce(nanmean) : Reduce(mean)
    end
    ignore_nan && return ReduceWithCoords(
        (A, g) ->
            nanmean(A, _match_eltype(weights(g), A); dims = g.dim_indices),
    )
    # Mirror average_lonlat's elementwise formula so the results match exactly
    return ReduceWithCoords() do A, g
        w = _match_eltype(weights(g), A)
        kept_dims = Tuple(filter(d -> d ∉ g.dim_indices, 1:ndims(A)))
        w_slice = dropdims(w, dims = kept_dims)
        norm = mapslices(A, dims = g.dim_indices) do slice
            T = float(eltype(slice))
            mask = ifelse.(isnan.(slice), T(NaN), one(T))
            mask .*= w_slice
            sum(mask)
        end
        return sum((A .* w) ./ norm, dims = g.dim_indices)
    end
end

# Weights follow the data's precision, so a Float32 variable stays Float32
_match_eltype(w, ::AbstractArray{T}) where {T <: AbstractFloat} =
    convert.(T, w)
_match_eltype(w, ::AbstractArray) = w

"""
    LatWeightedAverage(; ignore_nan = true)

Specify a `cos(latitude)`-weighted mean to apply to each block in a
split-apply-combine pattern; shorthand for
`Average(weights = group -> cosd.(group.lat))`.

Latitude must be one of the grouped dimensions and is assumed to be in
degrees. Grouping longitude and latitude together computes the joint weighted
mean matching [`weighted_average_lonlat`](@ref) (not
[`weighted_average_lat`](@ref)), except that the data's floating-point
precision is preserved.

# Example

```julia
var |>
ClimaAnalysis.GroupAll("lon") |>
ClimaAnalysis.GroupAll("lat") |>
ClimaAnalysis.LatWeightedAverage() |>
ClimaAnalysis.combine
```
"""
LatWeightedAverage(; ignore_nan = true) =
    Average(; ignore_nan, weights = g -> cosd.(g.lat))

"""
    SplitApplyVar

Represents lazy split and apply operations before combining to produce the
resulting `OutputVar`.
"""
struct SplitApplyVar{V <: OutputVar, S <: Tuple, A}
    "The `OutputVar` to split, apply, and combine."
    var::V

    "Split operations, each partitioning one dimension into groups of indices."
    split_ops::S

    "An apply operation to apply to each block, or `nothing` if not yet set."
    apply_op::A
end

"""
    (split_op::AbstractSplitOperation)(var::OutputVar)

Lazily applies `split_op` to `var` as part of the split-apply-combine pattern.

The result can be passed to another `AbstractSplitOperation` to group more
dimensions, or to an [`AbstractApplyOperation`](@ref) and then to
[`combine`](@ref) to produce the final `OutputVar`.
"""
function (split_op::AbstractSplitOperation)(var::OutputVar)
    return SplitApplyVar(var, (_resolve(split_op, var),), nothing)
end

"""
    (split_op::AbstractSplitOperation)(split_apply_var::SplitApplyVar)

Lazily applies `split_op` to `split_apply_var`, grouping another dimension as
part of the split-apply-combine pattern.

Grouping a dimension that is already grouped is an error.
"""
function (split_op::AbstractSplitOperation)(split_apply_var::SplitApplyVar)
    isnothing(split_apply_var.apply_op) ||
        error("An apply operation is already set")
    (; var, split_ops) = split_apply_var
    resolved_op = _resolve(split_op, var)
    any(op -> op.dim_name == resolved_op.dim_name, split_ops) &&
        error("The dimension $(resolved_op.dim_name) is already grouped")
    return SplitApplyVar(var, (split_ops..., resolved_op), nothing)
end

"""
    (apply_op::AbstractApplyOperation)(split_apply_var::SplitApplyVar)

Lazily applies `apply_op` to `split_apply_var` as part of the
split-apply-combine pattern.

The result should be passed to [`combine`](@ref) to produce the final
`OutputVar`.

You cannot apply another `AbstractApplyOperation` to a `split_apply_var` that
already has an `AbstractApplyOperation` applied on it.
"""
function (apply_op::AbstractApplyOperation)(split_apply_var::SplitApplyVar)
    isnothing(split_apply_var.apply_op) ||
        error("An apply operation is already set")
    (; var, split_ops) = split_apply_var
    return SplitApplyVar(var, split_ops, apply_op)
end

# Resolving names at pipe time moves errors to the offending pipe stage
_resolve(op::AbstractSplitOperation, var) = op
_resolve(op::GroupBy, var) =
    GroupBy(find_corresponding_dim_name_in_var(op.dim_name, var), op.by, op.on)

"""
    Base.show(io::IO, split_apply_var::SplitApplyVar)

Print the pipeline of `split_apply_var`: the short name of the variable, each
split operation with its resolved dimension name, and the apply operation if
one is set.
"""
function Base.show(io::IO, split_apply_var::SplitApplyVar)
    print(io, "SplitApplyVar(", short_name(split_apply_var.var), ")")
    for op in split_apply_var.split_ops
        print(io, " |> ")
        _show_op(io, op)
    end
    isnothing(split_apply_var.apply_op) &&
        return print(io, "  (no apply operation)")
    print(io, " |> ")
    _show_op(io, split_apply_var.apply_op)
    return print(io, "  (call combine to materialize)")
end

_show_op(io::IO, op) = print(io, nameof(typeof(op)))
_show_op(io::IO, op::GroupBy) =
    print(io, "GroupBy(\"", op.dim_name, "\", ", op.by, ")")
_show_op(io::IO, op::Reduce) = print(io, "Reduce(", op.reduction, ")")

Base.show(io::IO, f::WithCoords) = print(io, "WithCoords(", f.f, ")")

"""
    combine(split_apply_var::SplitApplyVar)

Materialize the `OutputVar` with the split and apply operations.

This splits the data of the `OutputVar` into blocks (the Cartesian product of
the groups of each grouped dimension), applies the apply operation to each
block, and combines the results together. Each grouped dimension is kept,
with one coordinate per group: the group's smallest coordinate, or its
largest if the dimension is descending. The attributes and dimension
attributes are unchanged.
"""
function combine(split_apply_var::SplitApplyVar)
    (; var, split_ops, apply_op) = split_apply_var
    isnothing(apply_op) && error("No apply operation is defined")
    groupings = map(op -> _group_indices(var, op), split_ops)
    data = _combine_blocks(var, apply_op, groupings, Val(ndims(var.data)))
    ret_dims = deepcopy(var.dims)
    for g in groupings
        ret_dims[g.name] = g.coords
    end
    return remake(var, dims = ret_dims, data = data)
end

"""
    _combine_blocks(var::OutputVar, apply_op, groupings, ::Val{N})

Apply `apply_op` to each block of the Cartesian product of `groupings` over
the data of `var` and combine the results into a single array.
"""
function _combine_blocks(
    var::OutputVar,
    apply_op,
    groupings,
    ::Val{N},
) where {N}
    slot = ntuple(d -> findfirst(g -> g.idx == d, groupings), Val(N))
    sizes = ntuple(Val(N)) do d
        isnothing(slot[d]) && return size(var.data, d)
        return length(groupings[slot[d]].groups)
    end
    dim_indices = map(g -> g.idx, groupings)
    ret = nothing
    for block_idx in CartesianIndices(map(g -> length(g.groups), groupings))
        # Ranges rather than Colon()s keep the view type concrete
        src = ntuple(Val(N)) do d
            isnothing(slot[d]) && return 1:size(var.data, d)
            return groupings[slot[d]].groups[block_idx[slot[d]]]
        end
        dst = ntuple(Val(N)) do d
            isnothing(slot[d]) && return 1:size(var.data, d)
            return block_idx[slot[d]]:block_idx[slot[d]]
        end
        group =
            Group(dim_indices, _block_coords(var, groupings, block_idx, Val(N)))
        block = _apply(apply_op, view(var.data, src...), group)
        all(d -> size(block, d) == 1, dim_indices) || error(
            "The apply operation must keep the grouped dimensions at size 1",
        )
        # similar from the block, not var.data: the reduction may change the
        # element type (e.g. mean of integers)
        isnothing(ret) && (ret = similar(block, sizes))
        ret[dst...] = block
    end
    return ret
end

"""
    _block_coords(var::OutputVar, groupings, block_idx, ::Val{N})

Return a `NamedTuple` of the coordinates of each grouped dimension in the
block `block_idx`, keyed by conventional dimension name and reshaped to
broadcast against the block.
"""
function _block_coords(
    var::OutputVar,
    groupings,
    block_idx,
    ::Val{N},
) where {N}
    return NamedTuple(
        Symbol(conventional_dim_name(g.name)) => reshape(
            var.dims[g.name][g.groups[block_idx[k]]],
            ntuple(
                d -> d == g.idx ? length(g.groups[block_idx[k]]) : 1,
                Val(N),
            ),
        ) for (k, g) in enumerate(groupings)
    )
end
