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

# A split operation must have a `dim_name` field and implement
# `_group_indices(var, split_op)` returning a `DimGrouping`
abstract type AbstractSplitOperation end

abstract type AbstractApplyOperation end

# Split operations -----------------------------------------------------------

struct GroupBy{F, C} <: AbstractSplitOperation
    "Dimension to group along"
    dim_name::String

    "Function mapping a coordinate to a label, or `nothing` to drop the index"
    by::F

    """
    Function of the coordinates of `dim_name` returning the values passed to
    `by`; when `nothing`, `by` receives the coordinates, converted to dates
    for a time dimension with a `start_date` attribute
    """
    on::C
end

GroupBy(dim_name, by; on = nothing) = GroupBy(dim_name, by, on)

# The entire dimension as a single group; returns a GroupBy with
# on = identity, so the labels are the raw coordinates
function GroupAll(dim_name)
    return GroupBy(dim_name, Returns(:all), nothing)
end

# Time split into DJF/MAM/JJA/SON per year; returns a GroupBy
function SplitSeason()
    return GroupBy("time", ClimaAnalysis.Utils.find_season_and_year, nothing)
end

struct Bins{E <: AbstractVector}
    "Bin edges; bins are [edges[i], edges[i+1]), last bin closed on the right"
    edges::E

    function Bins(edges::E) where {E <: AbstractVector}
        length(edges) >= 2 || error("At least two bin edges are needed")
        issorted(edges) || error("Bin edges must be in increasing order")
        return new{E}(edges)
    end
end

# Map x to the index of its bin, or nothing when x is in no bin
function (bins::Bins)(x)
    i = searchsortedlast(bins.edges, x)
    i == length(bins.edges) && return (x ≈ last(bins.edges) ? i - 1 : nothing)
    return i == 0 ? nothing : i
end

# Values handed to `by`: the coordinates of dim_name, converted to DateTimes
# for a time dimension with a start_date attribute
# TODO: Rename this since I don't like the name of it
function _labeling_values(var::OutputVar, dim_name)
    dim_name = find_corresponding_dim_name_in_var(dim_name, var)
    if conventional_dim_name(dim_name) == "time" &&
       haskey(var.attributes, "start_date")
        _check_time_dim(var)
        # Not dates(var): that returns the separate date dimension if one exists
        start_date = Dates.DateTime(var.attributes["start_date"])
        return time_to_date.(start_date, var.dims[dim_name])
    end
    return var.dims[dim_name]
end

"""
    DimGrouping

The grouping of one dimension of an `OutputVar`, produced by
`_group_indices` from a `AbstractSplitOperation`.
"""
struct DimGrouping{G, L, C}
    "Index of the grouped dimension in the data"
    idx::Int

    "Name of the grouped dimension in the `OutputVar`"
    name::String

    "Vector of index groups"
    groups::G

    # TODO: Might want to remove this for now?
    "Label of each group"
    labels::L

    "Representative coordinate of each group"
    coords::C
end

# Group the dimension split_op.dim_name; any AbstractSplitOperation must
# implement this and return a DimGrouping
function _group_indices(var::OutputVar, split_op::GroupBy)
    dim_name = find_corresponding_dim_name_in_var(split_op.dim_name, var)
    coord = var.dims[dim_name]
    ndims(coord) == 1 ||
        error("Cannot group along the multidimensional dimension $dim_name")
    # TODO: Improve variable naming here
    labeled =
        isnothing(split_op.on) ? _labeling_values(var, dim_name) :
        split_op.on(coord)
    length(labeled) == length(coord) || error(
        "Labeling values for $dim_name have length $(length(labeled)), " *
        "expected $(length(coord))",
    )

    labels_per_index = map(split_op.by, labeled)
    buckets = OrderedDict{eltype(labels_per_index), Vector{Int}}()
    for (i, label) in pairs(labels_per_index)
        isnothing(label) || push!(get!(() -> Int[], buckets, label), i)
    end
    isempty(buckets) && error("Grouping $dim_name produced no groups")

    # Order and represent each group by its extreme coordinate, so results do
    # not depend on the input ordering and a monotone dimension stays monotone
    rev = issorted(coord, rev = true)
    pick = rev ? maximum : minimum # improve variable name here
    labels = collect(keys(buckets))
    representatives = [pick(view(coord, buckets[k])) for k in labels]
    # MergeSort keeps ties between duplicate coordinates deterministic
    perm = sortperm(representatives, rev = rev, alg = MergeSort)
    labels, representatives = labels[perm], representatives[perm]
    groups = _as_ranges([buckets[k] for k in labels])
    length(groups) == length(labels) == length(representatives) || error(
        "The groups, labels, and coordinates must have the same length",
    )
    return DimGrouping(
        var.dim2index[dim_name],
        dim_name,
        groups,
        labels,
        representatives,
    )
end

# Convert every group to a UnitRange when all groups are consecutive
# (all-or-nothing, to keep the element type of the groups uniform)
function _as_ranges(groups)
    all(group -> group == first(group):last(group), groups) || return groups
    return [first(group):last(group) for group in groups]
end

# Apply operations ------------------------------------------------------------

struct Reduce{F} <: AbstractApplyOperation
    "Function called as `f(block; dims)`, or `f(block, group)` for WithCoords"
    reduction::F
end

struct WithCoords{F}
    "Function called as `f(block, group)`"
    f::F
end

ReduceWithCoords(f) = Reduce(WithCoords(f))

struct Group{D <: Tuple, C <: NamedTuple}
    "Indices of the dimensions being reduced"
    dim_indices::D

    "Coordinates of the grouped dimensions, reshaped to broadcast against the block"
    coords::C
end

# Resolve a dimension name to its conventional one, so that group.lat and
# group.latitude are the same coordinates
function _conventional(name::Symbol)
    error("Not implemented")
end

# group.<dim> returns the coordinates of a grouped dimension; error otherwise
function Base.getproperty(group::Group, name::Symbol)
    error("Not implemented")
end

function Base.propertynames(group::Group)
    error("Not implemented")
end

# Call f as f(block; dims) for plain functions, f(block, group) for WithCoords
function _call(f, block, group)
    error("Not implemented")
end

function _call(f::WithCoords, block, group)
    error("Not implemented")
end

function _apply(apply_op::Reduce, block, group)
    error("Not implemented")
end

# Mean over each block; weights is nothing or a function of a Group returning
# an array broadcastable against the block. Returns a Reduce.
function Average(; ignore_nan = true, weights = nothing)
    error("Not implemented")
end

# Convert the weights to the data's floating-point precision, so a Float32
# variable stays Float32
function _match_eltype(weights, data)
    error("Not implemented")
end

# cos(latitude)-weighted mean; returns a Reduce
function LatWeightedAverage(; ignore_nan = true)
    error("Not implemented")
end

# Machinery -------------------------------------------------------------------

struct SplitApplyVar{V <: OutputVar, S <: Tuple, A}
    "The `OutputVar` to split, apply, and combine."
    var::V

    "Split operations, each partitioning one dimension into groups of indices."
    split_ops::S

    "An apply operation to apply to each block, or `nothing` if not yet set."
    apply_op::A
end

function (split_op::AbstractSplitOperation)(var::OutputVar)
    error("Not implemented")
end

# Grouping a dimension that is already grouped is an error
function (split_op::AbstractSplitOperation)(split_apply_var::SplitApplyVar)
    error("Not implemented")
end

function (apply_op::AbstractApplyOperation)(split_apply_var::SplitApplyVar)
    error("Not implemented")
end

# Resolve dimension names at pipe time, so errors point at the offending stage
function _resolve(split_op::AbstractSplitOperation, var)
    error("Not implemented")
end

function _resolve(split_op::GroupBy, var)
    error("Not implemented")
end

# Print the pipeline, e.g.
# SplitApplyVar(ta) |> GroupBy("time", year) |> Reduce(nanmean)  (call combine to materialize)
function Base.show(io::IO, split_apply_var::SplitApplyVar)
    error("Not implemented")
end

function _show_op(io::IO, op)
    error("Not implemented")
end

function _show_op(io::IO, op::GroupBy)
    error("Not implemented")
end

function _show_op(io::IO, op::Reduce)
    error("Not implemented")
end

function Base.show(io::IO, f::WithCoords)
    error("Not implemented")
end

function combine(split_apply_var::SplitApplyVar)
    error("Not implemented")
end

# Apply apply_op once per block of the Cartesian product of the groupings and
# combine the results into a single array
function _combine_blocks(
    var::OutputVar,
    apply_op,
    groupings,
    ::Val{N},
) where {N}
    error("Not implemented")
end

# NamedTuple of the coordinates of each grouped dimension in the block,
# keyed by conventional dimension name and reshaped to broadcast against it
function _block_coords(
    var::OutputVar,
    groupings,
    block_idx,
    ::Val{N},
) where {N}
    error("Not implemented")
end
