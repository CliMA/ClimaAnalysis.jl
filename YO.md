# Rewrite split-apply-combine: GroupBy, block reductions, named coordinates

## Context

[src/split_apply_combine.jl](src/split_apply_combine.jl) (255 lines) implements
an experimental split-apply-combine pattern for `OutputVar`. The pipe interface
stays; the machinery behind it does not. Today it has:

- two unrelated split types (`GroupAll`, `SplitSeason`), each needing its own
  `_create_groups` method, and no way for a user to supply a grouping rule,
- a hard limit of one split dimension — a second split op errors
  ([split_apply_combine.jl:121](src/split_apply_combine.jl#L121)),
- no way for a reduction to see coordinate values, so a `cos(lat)`-weighted
  average is impossible.

Goal: one general split abstraction, several dimensions grouped at once (block
reductions), and coordinate-aware reductions that are **agnostic to dimension
order and to dimension spelling**.

The API is documented as experimental and outside semver
([docs/src/split_apply_combine.md:3-9](docs/src/split_apply_combine.md#L3-L9)),
so behavior changes are acceptable and tests get rewritten to match.

This plan has been through four adversarial reviews (Julia correctness, API
design, tests/regressions, performance). Their findings are folded in; the
"Corrections from review" section at the end records claims that were refuted.

---

## The interface today

As shipped in v0.5.23. Seven exported names; the only entry point is the pipe:

```julia
var |> ClimaAnalysis.GroupAll("time") |> ClimaAnalysis.Reduce(mean) |> ClimaAnalysis.combine
```

| Exported | Kind | Definition |
|---|---|---|
| `AbstractSplitOperation` | abstract type | Supertype of split operations |
| `AbstractApplyOperation` | abstract type | Supertype of apply operations |
| `GroupAll` | struct, field `dim_name::String` | Whole dimension as a single group |
| `SplitSeason` | singleton struct | Time split into DJF/MAM/JJA/SON per year |
| `Reduce` | struct, field `reduction::F` | Reduces each group |
| `SplitApplyVar` | struct, fields `var`, `split_op`, `apply_op` | The lazy object |
| `combine` | function | Materializes the result |

Split and apply ops are **callable structs**
([lines 111-140](src/split_apply_combine.jl#L111-L140)). There is one `split_op`
slot, so a second split op raises `"A split operation is already set"`.

`Reduce(f)` calls `f(group_view; dims = dim_idx)` with `dim_idx` an **`Int`**;
`f` must keep the reduced dimension at size 1. `combine` keeps the split
dimension, takes each group's **first** coordinate, and passes attributes
through untouched ([lines 222-255](src/split_apply_combine.jl#L222-L255)).

`SplitSeason` goes through `_check_time_dim`, so it requires time units of
`"s"` and a `start_date`; it delegates to `Utils.split_by_season_across_time`,
which errors on repeated dates.

**What it cannot do:** supply a grouping rule at all, group more than one
dimension, or let a reduction see coordinate values.

---

## The interface after this change

### Split operations

| Signature | | Meaning |
|---|---|---|
| `GroupBy(dim_name, by; on = nothing)` | **new** | Group `dim_name` by the label `by(coordinate)` |
| `GroupAll(dim_name)` | | `= GroupBy(dim_name, Returns(:all), var -> var.dims[dim_name])` |
| `SplitSeason()` | | `= GroupBy("time", find_season_and_year)` |
| `Bins(edges)` | **new** | A labeler: `GroupBy("lat", Bins(-90:30:90))` |

Indices whose labels are `isequal` form a group; a `nothing` label drops the
index. Chained `GroupBy`s group several dimensions at once and `combine` reduces
over the Cartesian product of the groups — the block reduction. Grouping the
same dimension twice errors, at the offending pipe stage.

### Apply operations

| Signature | | Reduction is called as |
|---|---|---|
| `Reduce(f)` | | `f(block; dims)` — today's contract |
| `ReduceWithCoords(f)` = `Reduce(WithCoords(f))` | **new** | `f(block, group)` |
| `Average(; ignore_nan = true, weights = nothing)` | **new** | mean, optionally weighted by `weights(group)` |
| `LatWeightedAverage(; ignore_nan = true)` | **new** | `= Average(weights = g -> cosd.(g.lat))` |

Every reduction runs **once per block**, vectorized, on a raw view.

### Terminal

`combine(split_apply_var)` — unchanged.

### What each interface can express

| Task | Today | After |
|---|---|---|
| Time mean | `GroupAll("time") \|> Reduce(mean)` | same |
| Seasonal means | `SplitSeason() \|> Reduce(mean)` | same |
| **Annual means** | not expressible | `GroupBy("time", Dates.year) \|> Average()` |
| **Monthly climatology** | not expressible | `GroupBy("time", Dates.month) \|> Average()` |
| **30° latitude bands** | not expressible | `GroupBy("lat", Bins(-90:30:90))` |
| **Tropics vs. extratropics** | not expressible | `GroupBy("lat", l -> abs(l) < 30)` |
| **30°×60° block means** | not expressible | two chained `GroupBy`s |
| **`cos(lat)`-weighted mean** | not expressible | `GroupAll("lon") \|> GroupAll("lat") \|> LatWeightedAverage()` |
| **Coordinate-dependent reduction** | not expressible | `ReduceWithCoords(f)` |

`GroupBy("time", Dates.year)` works because a time dimension with a
`start_date` is labeled on dates rather than raw seconds. `Dates.month`,
`Dates.dayofyear`, `Dates.week` and `find_season_and_year` come along for free;
`SplitSeason` is just the last of those with a name.

---

## Part 1 — split

```julia
struct GroupBy{F, C} <: AbstractSplitOperation
    "Dimension to group along"
    dim_name::String
    "Function mapping a coordinate value to a label, or `nothing` to drop it"
    by::F
    "Function `var -> vector` supplying the values to label; `nothing` uses the coordinates"
    on::C
end
GroupBy(dim_name, by; on = nothing) = GroupBy(dim_name, by, on)

GroupAll(dim_name) = GroupBy(dim_name, Returns(:all), var -> var.dims[dim_name])
SplitSeason() = GroupBy("time", find_season_and_year)
```

`GroupAll` supplies its own `on` deliberately: it never inspects the labels, and
routing through the date conversion would make it fail on an integer time
dimension, since `time_to_date` only accepts an `AbstractFloat`
([Utils.jl:532](src/Utils.jl#L532)).

`on` also lets a dimension be labeled by something other than its own
coordinates, which is the paper's "split by the values of a variable" case:

```julia
GroupBy("time", Dates.year)                                        # dates, by default
GroupBy("time", x -> fld(x, 86400), on = var -> var.dims["time"])  # raw seconds
GroupBy("time", identity, on = var -> enso_phase)                  # external labels
```

### Labeling values

```julia
function _labeling_values(var, dim_name)
    conventional_dim_name(dim_name) == "time" &&
        haskey(var.attributes, "start_date") || return var.dims[dim_name]
    _check_time_dim(var)
    # Deliberately not `dates(var)`: it returns a separate `date` dimension when
    # one exists, a different vector of the same length
    return time_to_date.(
        Dates.DateTime(var.attributes["start_date"]),
        var.dims[dim_name],
    )
end
```

Keeping `_check_time_dim` preserves the `units == "s"` check
([Var.jl:2100](src/Var.jl#L2100)); without it a file storing time in days
produces wrong seasons silently.

### Binning: explicit edges, validated

```julia
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
    return (i == 0 || i > length(b.edges) - 1) ? nothing : i
end
```

Both constructor guards are load-bearing. `searchsortedlast` assumes ascending
order, so `Bins(90.0:-30.0:-90.0)` silently labels only one of seven latitudes
and drops the rest. `Bins([0.0])` yields a bin index of `0` that becomes a real
group, and `Bins(Float64[])` throws a `BoundsError` from `last`.

Verified truth table for ascending edges: below first → `nothing`; on first edge
→ bin 1; on last edge → last bin; above last → `nothing`; `NaN` → `nothing`.

There is deliberately **no** `GroupBy(dim, width::Real)` overload —
`fld.(-90:30:90, 30)` gives seven singleton groups, and on a 2.5° grid the pole
becomes a group of one. Nor a `Dates.Period` overload:
`floor(DateTime(2010, 12, 15), Month(3))` is October, so December groups with
Oct/Nov and DJF is torn in half.

### Group construction

```julia
function _group_indices(var, op::GroupBy)
    name = op.dim_name   # resolved at pipe time
    coord = var.dims[name]
    ndims(coord) == 1 || error("Cannot group along the multidimensional dimension $name")
    labeled = isnothing(op.on) ? _labeling_values(var, name) : op.on(var)
    length(labeled) == length(coord) ||
        error("Labeling vector for $name has length $(length(labeled)), " *
              "expected $(length(coord))")

    labels_per_index = map(op.by, labeled)
    buckets = OrderedDict{eltype(labels_per_index), Vector{Int}}()
    for (i, label) in pairs(labels_per_index)
        isnothing(label) || push!(get!(() -> Int[], buckets, label), i)
    end
    isempty(buckets) && error("Grouping $name produced no groups")

    # Represent and order groups so that a monotone dimension stays monotone
    rev = issorted(coord, rev = true)
    pick = rev ? maximum : minimum
    labels = collect(keys(buckets))
    reps = [pick(view(coord, buckets[k])) for k in labels]
    perm = sortperm(reps, rev = rev, alg = MergeSort)
    labels, reps = labels[perm], reps[perm]
    groups = [_as_range(buckets[k]) for k in labels]
    return (; idx = var.dim2index[name], name, groups, labels, coords = reps)
end

_as_range(v) = v == first(v):last(v) ? (first(v):last(v)) : v
```

Three details are load-bearing:

- The bucket dict is **concretely typed** from `eltype(labels_per_index)`. An
  `OrderedDict{Any, …}` boxes every label and does a dynamic `hash`/`isequal`
  per coordinate — the one place per-coordinate dynamic dispatch would survive.
- Representatives are **precomputed** and sorted with `sortperm`. Passing
  `by = k -> pick(view(coord, buckets[k]))` to `sort!` recomputes an O(|group|)
  `minimum` plus a dict lookup on every comparison.
- `alg = MergeSort` makes ties deterministic. Ties are now reachable, because
  duplicate coordinates no longer error.

Ordering and representing groups by their extreme coordinate (rather than by
first occurrence) is what makes results independent of input ordering — required
by the shuffled-dates test at
[test_split_apply_combine.jl:122-136](test/test_split_apply_combine.jl#L122-L136).
`_as_range` matters because seasonal, monthly and binned groups on a sorted
dimension are always consecutive, and a `UnitRange` view is affine where a
`Vector{Int}` view is a gather.

---

## Part 2 — apply

```julia
struct Reduce{F} <: AbstractApplyOperation
    "Function called as `f(block; dims)`"
    reduction::F
end

struct WithCoords{F}
    "Function called as `f(block, group)`"
    f::F
end

_call(f, block, group) = f(block; dims = group.dim_indices)
_call(f::WithCoords, block, group) = f.f(block, group)

_apply(op::Reduce, block, group) = _call(op.reduction, block, group)

ReduceWithCoords(f) = Reduce(WithCoords(f))
```

`Reduce(mean)`, `Reduce(nanmean)`, `Reduce(std)`, `Reduce(minimum)` all work
verbatim — this is today's contract, `f(block; dims)` returning an array with
the grouped axes at size 1. The only change is that `dims` is **always a
`Tuple`**, since several dimensions can be grouped.

Marking the *function* rather than the operation is what makes this scale: the
paper's other two apply categories are future work, and they get coordinate
support for free — `Transform(WithCoords(f))` — instead of needing a
`TransformWithCoords` and a `FilterWithCoords`.

### What `group` carries

```julia
struct Group{D <: Tuple, C <: NamedTuple}
    "Indices of the axes being reduced"
    dim_indices::D
    "Coordinates of those axes, reshaped to broadcast against the block"
    coords::C
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

const _CONVENTIONAL = NamedTuple(
    Symbol(name) => Symbol(conventional_dim_name(name)) for name in Iterators.flatten((
        LONGITUDE_NAMES, LATITUDE_NAMES, TIME_NAMES, DATE_NAMES,
        ALTITUDE_NAMES, PRESSURE_NAMES,
    ))
)
_conventional(name::Symbol) =
    hasfield(typeof(_CONVENTIONAL), name) ? getfield(_CONVENTIONAL, name) : name
```

The `fieldnames` guard rather than `name === :dim_indices` is required: with only
the latter, `g.coords` routes through `_conventional`, misses, and errors — which
would break `values(g.coords)`, destructuring, and `propertynames`.

`_CONVENTIONAL` is a `const` **NamedTuple**, not a `Dict`. A `Dict` reads mutable
fields, so `get` is not `:foldable` and the lookup survives into the generated
code; `getfield` on a `const` immutable with a literal symbol folds away. Built
from the existing name constants
([outvar_dimensions.jl:2-7](src/outvar_dimensions.jl#L2-L7)) so the lists are not
duplicated.

Only grouped dimensions appear in `coords`; asking for anything else errors and
lists what is available. Weighting a dimension that is not being reduced cancels
out of a normalized mean, so it is almost always a bug.

### Agnosticism

**Dimension order.** Coordinates are reshaped to broadcast against the block, so
`cosd.(g.lat)` lands on the right axis wherever latitude sits. `vec(g.lat)`
recovers the plain vector.

**Dimension name.** `g.lat` and `g.latitude` are one thing — but **only those
two**: `LATITUDE_NAMES == ["lat", "latitude"]`
([outvar_dimensions.jl:3](src/outvar_dimensions.jl#L3)), so `g.lati` does *not*
resolve. On the split side, `GroupBy("lat", …)` resolves through
`find_corresponding_dim_name_in_var`
([outvar_dimensions.jl:298](src/outvar_dimensions.jl#L298)) to whatever the file
calls it.

**Naming a dimension at all.** `values(g.coords)` is ordered to match the grouped
dimensions:

```julia
ReduceWithCoords() do A, g
    w = prod(_cell_width.(values(g.coords)))
    sum(A .* w; dims = g.dim_indices) ./ sum(w; dims = g.dim_indices)
end
```

### A worked example, with shapes

For a `(time, lat, lon)` variable of size `(120, 181, 361)`:

```julia
var |> GroupAll("lat") |> GroupAll("lon") |> LatWeightedAverage() |> combine
```

| | Value | Size |
|---|---|---|
| `A` | the block — the whole array, since `GroupAll` makes one group per dimension | `(120, 181, 361)` |
| `g.dim_indices` | `(2, 3)` — found by name, never written by the user | |
| `g.lat` | latitudes, reshaped to broadcast along axis 2 | `(1, 181, 1)` |
| `w` | `cosd.(g.lat)` | `(1, 181, 1)` |
| `nanmean(A, w; dims = g.dim_indices)` | grouped axes kept at size 1 | `(120, 1, 1)` |

Both dimensions reduce in one call, so this is the joint weighted average
matching `weighted_average_lonlat`, not an average of averages. On a
`(lat, lon, time)` variable only `g.dim_indices` and the reshape of `g.lat`
change; the reduction's source is identical.

### The batteries

```julia
function Average(; ignore_nan = true, weights = nothing)
    if isnothing(weights)
        return ignore_nan ? Reduce(nanmean) : Reduce(mean)
    end
    ignore_nan && return ReduceWithCoords(
        (A, g) ->
            nanmean(A, _match_eltype(weights(g), A); dims = g.dim_indices),
    )
    # Mirrors average_lonlat: normalize elementwise, then sum
    return ReduceWithCoords() do A, g
        w = _match_eltype(weights(g), A)
        norm = mapslices(A; dims = g.dim_indices) do slice
            sum(ifelse.(isnan.(slice), NaN, 1.0) .* w)
        end
        return sum((A .* w) ./ norm; dims = g.dim_indices)
    end
end

# Weights follow the data's precision, so a Float32 variable stays Float32
_match_eltype(w, ::AbstractArray{T}) where {T <: AbstractFloat} =
    convert.(T, w)
_match_eltype(w, ::AbstractArray) = w

LatWeightedAverage(; ignore_nan = true) =
    Average(; ignore_nan, weights = g -> cosd.(g.lat))
```

Four things here are not stylistic:

- `ignore_nan ? Reduce(nanmean) : Reduce(mean)` rather than
  `Reduce(ignore_nan ? nanmean : mean)`. The latter is a runtime branch over two
  distinct function types, so `F` is not inferred and `Reduce{F}` is abstract.
- Written as a `function`, not a ternary whose branch is a `do` block. It parses,
  but does not round-trip through the pinned `JuliaFormatter =2.10.1`.
- The `ignore_nan = false` branch mirrors `average_lonlat`
  ([Var.jl:975-987](src/Var.jl#L975-L987)) elementwise so the regression test can
  assert `==`. The obvious `sum(A .* w; dims) ./ sum(w; dims)` is **wrong**: `w`
  has the block's shape only along the weighted axis, so the denominator is short
  by the product of the sizes of the other grouped dimensions.
- `_match_eltype` narrows the weights to the data's precision. Dimensions read
  from NetCDF are `Float64` even when the data is `Float32`, so `cosd.(g.lat)` is
  `Float64` and the plain `nanmean(A, W)` promotes the whole result — a `Float32`
  variable comes back `Float64` at twice the memory. Narrowing only ever applies
  to floating-point data: the `AbstractArray` fallback leaves `Float64` weights
  alone on integer data, where converting them would be catastrophic.

`weights` is a function of the group returning anything broadcastable against
the block:

```julia
Average()
Average(weights = g -> cosd.(g.lat) .^ 2)
LatWeightedAverage()
```

`nanmean(A, W; dims)` is native to NaNStatistics and is the exact call the repo
already makes for `weighted_average_lonlat`
([Var.jl:973](src/Var.jl#L973)), so on the `ignore_nan = true` path
`LatWeightedAverage()` reproduces existing semantics bitwise **for `Float64`
data** — the overwhelmingly common case, and the one the regression test should
assert with `==`.

For `Float32` data the two now differ in the last bits, deliberately:
`weighted_average_lonlat` returns `Float64`, `LatWeightedAverage()` returns
`Float32`. That test needs `≈`, and it is the one place where "matches the
existing function exactly" and "does not silently double the memory of a `Float32`
variable" cannot both hold. Record it in NEWS.

Two limits on the narrowing:

- It covers only `Average`'s `weights` keyword. A hand-written
  `ReduceWithCoords((A, g) -> nanmean(A, cosd.(g.lat); dims = g.dim_indices))`
  promotes exactly as before; say so in the `ReduceWithCoords` docstring.
- The `ignore_nan = false` branch still promotes through its mask, whose `NaN` and
  `1.0` literals are `Float64`. Writing them as `T(NaN)` and `one(T)` belongs with
  the `dropdims` fix that branch needs anyway (appendix A.2), not here.

> `LatWeightedAverage` matches `weighted_average_lonlat`, **not**
> `weighted_average_lat`, which builds a per-column NaN-masked normalization by
> hand ([Var.jl:828-841](src/Var.jl#L828-L841)). Say so in the docstring.
>
> `average_lonlat` also warns when latitudes look like radians
> ([Var.jl:963-964](src/Var.jl#L963-L964)). Carry that check over.

---

## Part 3 — the rest of the machinery

```julia
struct SplitApplyVar{V <: OutputVar, S <: Tuple, A}
    var::V
    split_ops::S
    apply_op::A
end

(op::AbstractSplitOperation)(var::OutputVar) =
    SplitApplyVar(var, (_resolve(op, var),), nothing)

function (op::AbstractSplitOperation)(sav::SplitApplyVar)
    isnothing(sav.apply_op) || error("An apply operation is already set")
    op = _resolve(op, sav.var)
    any(o -> o.dim_name == op.dim_name, sav.split_ops) &&
        error("The dimension $(op.dim_name) is already grouped")
    return SplitApplyVar(sav.var, (sav.split_ops..., op), nothing)
end

function (op::AbstractApplyOperation)(sav::SplitApplyVar)
    isnothing(sav.apply_op) || error("An apply operation is already set")
    return SplitApplyVar(sav.var, sav.split_ops, op)
end

# Resolving names eagerly moves errors from `combine` to the offending pipe stage
_resolve(op::AbstractSplitOperation, var) = op
_resolve(op::GroupBy, var) =
    GroupBy(find_corresponding_dim_name_in_var(op.dim_name, var), op.by, op.on)
```

The `AbstractSplitOperation` fallback for `_resolve` is needed or any
user-defined subtype hits a `MethodError`.

A `Tuple` rather than a `Vector` for `split_ops`: a
`Vector{AbstractSplitOperation}` boxes the operations and dispatches `by`
dynamically once per coordinate.

### Printing the pipeline

`SplitApplyVar` has no `show` today, so a forgotten `combine` prints the default
struct dump — which recurses into `Base.show(io::IO, ::OutputVar)`
([Var.jl:2919](src/Var.jl#L2919)) and spills every attribute and dimension
attribute of the variable, wrapped in the struct's type parameters. Ten lines fix
it:

```julia
function Base.show(io::IO, sav::SplitApplyVar)
    print(io, "SplitApplyVar(", short_name(sav.var), ")")
    for op in sav.split_ops
        print(io, " |> ")
        _show_op(io, op)
    end
    isnothing(sav.apply_op) && return print(io, "  (no apply operation)")
    print(io, " |> ")
    _show_op(io, sav.apply_op)
    return print(io, "  (call combine to materialize)")
end

_show_op(io::IO, op) = print(io, nameof(typeof(op)))
_show_op(io::IO, op::GroupBy) =
    print(io, "GroupBy(\"", op.dim_name, "\", ", op.by, ")")
_show_op(io::IO, op::Reduce) = print(io, "Reduce(", op.reduction, ")")
Base.show(io::IO, f::WithCoords) = print(io, "WithCoords(", f.f, ")")
```

```
julia> var |> GroupBy("t", Dates.year)
SplitApplyVar(ta) |> GroupBy("time", year)  (no apply operation)

julia> var |> GroupBy("t", Dates.year) |> Average()
SplitApplyVar(ta) |> GroupBy("time", year) |> Reduce(nanmean)  (call combine to materialize)
```

Four notes:

- Only the two-argument `show`, matching `Base.show(io::IO, var::OutputVar)`;
  `display` falls back to it.
- It deliberately does **not** compute groups. Name errors are already caught at
  pipe time by `_resolve`, but `"Grouping $name produced no groups"` is not, and a
  `show` that throws makes the variable unprintable. Group counts stay unavailable
  until `combine`.
- Because `_resolve` runs at pipe time, `op.dim_name` is the **resolved** name:
  the user typed `"t"` and the pipeline reports `"time"`. That is the cheapest
  possible confirmation that name resolution did what they expected.
- `Average()` prints as `Reduce(nanmean)`, since `Average` is a constructor rather
  than a type. Mildly surprising, but it discloses the NaN policy — the one place
  where `Average()` and `Reduce(mean)` quietly disagree.

`short_name` returns `""` for a variable without one
([Var.jl:488-490](src/Var.jl#L488-L490)), so nothing here can throw. Anonymous
`by` and `weights` functions print as `#7`; unavoidable.

### combine

```julia
function combine(sav::SplitApplyVar)
    (; var, split_ops, apply_op) = sav
    isnothing(apply_op) && error("No apply operation is defined")
    gs = map(op -> _group_indices(var, op), split_ops)
    data = _combine_blocks(var, apply_op, gs, Val(ndims(var.data)))
    ret_dims = deepcopy(var.dims)
    for g in gs
        ret_dims[g.name] = g.coords
    end
    # Attributes and dimension attributes are passed through, as they are today
    return remake(var, dims = ret_dims, data = data)
end

function _combine_blocks(var, apply_op, gs, ::Val{N}) where {N}
    slot = ntuple(d -> findfirst(g -> g.idx == d, gs), Val(N))
    sizes = ntuple(d -> isnothing(slot[d]) ? size(var.data, d) : length(gs[slot[d]].groups), Val(N))
    dim_indices = map(g -> g.idx, gs)
    ret = nothing
    for I in CartesianIndices(map(g -> length(g.groups), gs))
        src = ntuple(d -> isnothing(slot[d]) ? Colon() : gs[slot[d]].groups[I[slot[d]]], Val(N))
        dst = ntuple(d -> isnothing(slot[d]) ? Colon() : (I[slot[d]]:I[slot[d]]), Val(N))
        block = _apply(apply_op, view(var.data, src...), Group(dim_indices, _block_coords(var, gs, I, Val(N))))
        all(d -> size(block, d) == 1, dim_indices) ||
            error("The apply operation must keep the grouped dimensions at size 1")
        isnothing(ret) && (ret = similar(block, sizes))
        ret[dst...] = block
    end
    return ret
end

function _block_coords(var, gs, I, ::Val{N}) where {N}
    return NamedTuple(
        Symbol(conventional_dim_name(g.name)) => reshape(
            var.dims[g.name][g.groups[I[k]]],
            ntuple(d -> d == g.idx ? length(g.groups[I[k]]) : 1, Val(N)),
        ) for (k, g) in enumerate(gs)
    )
end
```

Details that are load-bearing:

- `_apply` takes **three** arguments; `dims` lives on the `Group`.
- `dim_indices` is always a `Tuple`, never `only(idxs)` — `Group{D <: Tuple}`
  cannot hold a scalar, and NaNStatistics normalizes `Int` to `Tuple` anyway
  (`_nanmean(A, dims::Int, st) = _nanmean(A, (dims,), st)`), so the `Union` buys
  nothing and costs a dynamic reduction call per block.
- `similar(block, sizes)`, **not** `similar(var.data, sizes)` — the latter throws
  `InexactError` for `Reduce(mean)` on integer data.
- The shape check is explicit. `ret[:, 1:1, :] = block` silently accepts a block
  that dropped the reduced dimension, so `setindex!` does not enforce the
  contract on its own.
- `_block_coords` **materializes** `var.dims[g.name][g.groups[…]]` before
  reshaping. Reshaping a view of non-consecutive indices does not error, but
  gives a `ReshapedArray` over an `IndexCartesian` parent; materializing puts a
  concrete `Array` in the NamedTuple and keeps `Group`'s type parameter clean.
- The keys are `Symbol(conventional_dim_name(g.name))`, not `Symbol(g.name)`, or
  `_conventional` lookup misses.
- `Val(N)` rather than a runtime `N`, so the index tuples' element types are
  inferred. `N = ndims(var.data)` is a type-domain constant, so this is free.
- The block loop lives in its own function, but **that is not enough to make the
  `Group`'s `C` parameter concrete** — measured, `_combine_blocks` as written above
  infers `ret::Any`. See "Type stability — measured" under Performance; the three
  causes and their cost are quantified there.
- `deepcopy(var.dims)`, matching
  [line 250](src/split_apply_combine.jl#L250): the result must own its dims.
  `remake` deep-copies only *omitted* keyword arguments, so a passed `dims` is
  used verbatim.

---

## Coverage of Wickham (2011)

Checked against ["The Split-Apply-Combine Strategy for Data
Analysis"](https://www.jstatsoft.org/article/view/v040i01), cited by the docs.

### Split

| Paper | Plan |
|---|---|
| Split an array by one or more **margins** | ✅ Chained `GroupBy`s. The framings are duals: plyr's margins are the dimensions you *keep*, a `GroupBy` dimension is one you *reduce*. |
| Split by the values of one or more **variables** | ✅ Coordinates are the grouping variables and `by` is the rule; several variables on one axis is a tuple-valued label. `on` covers grouping by a variable that is not the coordinate. |
| Combinations across axes | ✅ Chained `GroupBy`s form the Cartesian product. |
| `.drop` | ➖ Moot: labels only exist where data exists, and a product of non-empty groups is non-empty. `nothing` from `by` is the deliberate way to discard. |
| Split a list | ➖ No analogue for `OutputVar`. |

### Split — overlapping windows (future work, but constrains the design now)

Moving averages need element *i* to belong to *k* different windows. `GroupBy`
maps one coordinate to **one** label and `_group_indices` pushes each index into
exactly one bucket, so its groups are a partition by construction. No `by` can
produce overlapping groups — this is an impossibility, not an awkward spelling.

Note that **non**-overlapping windows are already expressible, so only rolling
operations are affected:

```julia
# 3-month block means / coarsening / resample-to-annual
var |> GroupBy("time", i -> fld(i - 1, 3), on = var -> eachindex(var.dims["time"])) |>
    Average() |> combine
```

This is a capability the shipped `_create_groups` and the `group_and_reduce_by`
prototype in `src/temp.jl` both had — the latter's docstring says explicitly
"`group_by` does not need to partition the values of the dimension" — and which
the labeler-based `GroupBy` gives up. Worth recording as a deliberate trade.

**`combine` already supports overlapping groups.** In `_combine_blocks`, each
group index `I` reads from `src` and writes to `dst = I:I`. Sources may overlap
freely; only destinations must be disjoint, and they are by construction. Nothing
in the block loop assumes a partition.

So the only missing piece is a second split operation that yields index groups
directly rather than through labels:

```julia
struct Rolling <: AbstractSplitOperation
    "Dimension to roll along"
    dim_name::String
    "Number of coordinates per window"
    width::Int
    "Step between successive windows"
    stride::Int
end
Rolling(dim_name, width; stride = 1) = Rolling(dim_name, width, stride)
```

Prefer `Rolling` over `Window`/`Windows`: the package already has
`window(var, dim_name; left, right)`
([outvar_selectors.jl:201](src/outvar_selectors.jl#L201)) for *subsetting* a
dimension, and reusing the word for a rolling reduction would collide. `Rolling`
also matches the vocabulary this audience knows from pandas and xarray.

**What this requires of the design as written — four things to preserve:**

1. **Keep `AbstractSplitOperation` as a real extension point.** Rolling windows
   are the concrete second split operation with genuinely different semantics, so
   the abstract type is not vestigial. It needs a documented contract:
   `_group_indices(var, op)` must return `(; idx, name, groups, labels, coords)`.
2. **Keep `_group_indices` dispatching on the operation type.** Do not inline the
   labeler-bucketing logic into `combine`; `Rolling` supplies its own method.
3. **Do not assume groups partition anywhere in `combine`.** They do not today,
   and nothing added later should start to.
4. **Let each split operation choose its own representative coordinate.** With
   data-derived coordinates a centred rolling mean would otherwise be labelled by
   the window's *first* coordinate. `Rolling`'s `_group_indices` should set
   `coords` to the window centre directly — a per-operation decision needing no
   global change.

Given those, `Rolling` slots in with **zero changes** to `combine`, `Reduce`,
`Average`, `LatWeightedAverage` or `Group`, and `Reduce(mean)` /
`LatWeightedAverage()` work over windows unmodified.

Two caveats to settle when it lands: edge handling (`valid` gives `n - k + 1`
windows, `same` pads to `n`, and only the former needs no policy), and cost —
each window reduces `k` elements independently, so this is `O(n·k)` rather than
the `O(n)` a cumulative-sum implementation would achieve. Fine for a 12-month
window over a decade; wasteful for a wide window over a long record.

A rolling reduction is **not** the deferred `Transform` case below; it is a
reduction over overlapping groups. Both pandas and xarray keep `rolling` separate
from `groupby` for the same reason.

### Apply

| Paper | Plan |
|---|---|
| **Reductions / summaries** | ✅ `Reduce`, `ReduceWithCoords`, `Average`, `LatWeightedAverage` |
| **Transformations** (size preserving), e.g. per-season z-score | ❌ Future work |
| **Filters** (top *k* per group) | ❌ Future work |

Leave the seam: make the output size a function of the apply operation.

```julia
_out_size(::Reduce, var, gs, d, slot) = length(gs[slot].groups)
_out_size(::Transform, var, gs, d, slot) = size(var.data, d)
```

Coordinate support comes free — `WithCoords` marks the function, so
`Transform(WithCoords(f))` needs no new type.

### Combine

Concatenation ordered by group: ✅, and independent of input ordering. Retaining
group labels in the output, which the paper treats as part of combining: ⚠️
deferred, see below.

### Not adopted

`.progress`, `.parallel`, `each()`, `colwise`. plyr also accommodates apply
functions returning different shapes per group; this plan requires a uniform
shape and errors otherwise, which is what makes a single preallocated output
possible.

---

## Group metadata — DEFERRED

> **Not part of this implementation.** The shape is still being iterated on.
> Nothing else depends on it: one field or one type on the split side and one
> line in `combine`, layerable later.

**Decided:** all four kinds are wanted eventually — the label, the group size,
the coordinate range, and `long_name` provenance — primarily for debuggability.
**Where** they live is open.

The label comes from `by`; the name is only the attribute key. Grouping leaves
exactly one label per group in output order (`g.labels`), so recording them is
`ret.attributes[name] = g.labels`.

```julia
seasonal.attributes["season"]  # [("DJF", 2010), ("MAM", 2010), ...]
```

The payoff: [`average_season_across_time`](src/Var.jl#L2123-L2159) — 37 lines of
hand-rolled `cat` and attribute plumbing — collapses into this API.

Three shapes, none committed to:

```julia
GroupBy("time", find_season_and_year, label_name = "season")   # keyword on GroupBy
Labeled(GroupBy("time", find_season_and_year), "season")        # wrapper split op
combine(sav, labels = "season")                                 # keyword on the terminal
```

### The wider design space

| Available per group | Why it matters |
|---|---|
| **Label** | Which season/year/band this slice is |
| **Group size** | Distinguishes a complete DJF from one with two months — currently invisible, and a silent bias when averaging across seasons |
| **Coordinate range** | The span each slice covers, not just where it starts |
| **Grouping rule** | Provenance |

Where it goes: `attributes` (matches `average_season_across_time`),
`dim_attributes[dim_name]` (structurally honest — these are per-element metadata
of the split dimension), or `long_name` provenance via
`_update_long_name_generic!` ([Var.jl:1208-1252](src/Var.jl#L1208-L1252)), the
repo's existing idiom for recording that a reduction happened.

A typing constraint: `OutputVar.attributes` is `Dict{String, B}`, so writing a
`Vector{Tuple{String, Int}}` into all-`String` attributes needs the dict rebuilt
with a promoted value type. `average_season_across_time` already works around
this ([Var.jl:2146-2149](src/Var.jl#L2146-L2149)).

### Coordinate bounds (future work, but decide the shape now)

The "coordinate range" row above is not just a debugging nicety — it is the CF
`bounds` convention, and it is the difference between a coarsened variable that
downstream code can use correctly and one it silently misreads.

After grouping, **one coordinate value stands for an interval**. Nothing in the
output records that. A `GroupBy("lat", Bins(-90:30:90))` result is six latitudes
with no indication that each spans 30°, so `integrate_lat`, `weighted_average_lat`,
area weighting and every plotting recipe treat them as point samples and infer
spacing from the gaps between representatives. CF solves this with a `bounds`
attribute on the coordinate variable naming an `(n, 2)` companion variable —
`time_bnds`, `lat_bnds`. `OutputVar` has nowhere to put one.

**It is nearly free to compute.** `_group_indices` already scans each group's
coordinates to pick the representative; `extrema(view(coord, buckets[k]))` returns
`(lo, hi)` from the same scan. Ask for both in that one pass rather than
retrofitting a second traversal later.

Where it could live, in increasing order of invasiveness:
`dim_attributes[dim_name]["bounds"]` as an `(n, 2)` array (closest to CF, and it
hits the same value-type promotion problem as the label channel, since
`dim_attributes` is `OrderedDict{String, C}` with `C <: AbstractDict`); a
`bounds` field on `OutputVar` (touches every constructor, `remake`,
`arecompatible`, `flatten`); or a companion dimension.

**Five things to preserve now so this stays layerable:**

1. `_group_indices` must keep each group's index set in hand while it computes the
   representative. It does — do not optimize that away.
2. **Do not derive the representative coordinate from the label.** With
   data-derived coordinates the representative and the bounds come from the same
   scan and are consistent by construction. This is a second, independent argument
   for the decision already taken in the appendix.
3. Whatever holds bounds must survive `remake`, since `combine` builds its result
   that way.
4. `Rolling` needs this more than anything else does: a centred window's
   representative coordinate loses the width entirely, so a rolling mean without
   bounds is indistinguishable from a point series.
5. Decide what `bounds` *means* before writing any of it. For
   `GroupBy("lat", Bins(-90:30:90))` the honest data bounds are the extreme
   latitudes actually present (say −89 and −61), while the bin edges are −90 and
   −60. Both are useful and they are different numbers. For
   `GroupBy("time", Dates.year)` there are no edges at all, only data. So bounds
   are per-split-operation, the same way the representative coordinate is —
   which is a third reason to let each split operation decide, rather than
   computing either one centrally.

One caveat that argues for keeping this modest: exact bounds need the *cell*
extent, not the extreme cell centres. Recovering −90 from a centre at −89 requires
assuming uniform spacing, and the appendix already notes the package is
inconsistent about that (`integrate_lonlat` uses `Δφ·cos φ` while
`average_lonlat(weighted = true)` uses bare `cos φ`). Record data bounds, which are
cheap and exact; do not synthesize cell bounds.

---

## Concerns from the original PR

| Concern | Status |
|---|---|
| 1. Should this be in a module? | Unchanged — stays in `Var`, exported |
| 2. Docstrings for structs used as functions | Unchanged, and slightly worse: `Bins` and `WithCoords` add two more callable structs. Consider making them closures instead, which is how `Template` avoids the problem ([Template.jl:139-141](src/Template.jl#L139-L141)) |
| 3. No reductions that need dimension values | **Resolved** — `ReduceWithCoords` |
| 4. Errors at materialization, not construction | **Resolved** for name errors and duplicate grouping, via `_resolve` at pipe time. Shape mismatches can only be caught in `combine` |
| 5. New coordinate is always the group's first element | **Changed** — now the extreme coordinate, so results do not depend on input ordering. Identical for ascending dimensions |
| 6. Attributes not updated | Unchanged; the label channel is deferred |
| Future: arbitrary splits | **Now the core** (`GroupBy`) |
| Future: transformations and filters | Still future; `_out_size` and `WithCoords` are the seams |

The PR's motivation — keeping the reduced dimension so downstream calibration
code such as `ObservationRecipe` still sees a time dimension — is preserved:
`combine` always keeps every grouped dimension.

---

## Performance

**Every reduction is vectorized**, called once per block on a raw view, so
NaNStatistics' kernels are reached normally. There are no `Array`-only fast paths
on the `nanmean` path that a `SubArray` would miss.

**Allocations are unchanged, not reduced.** `similar(block, sizes)` plus
`ret[dst...] = block` *is* the `cat` copy under another name, and each block's
reduction still allocates its own result. What genuinely improves is peak live
memory — today all G block results are alive at once *plus* the `cat` output —
and `cat`'s overhead: `cat(reduced..., dims)` splats a runtime-length `Vector`
through `Core._apply_iterate` with no specialization, and past 32 groups hits
`Base.Any32` paths. For a 120-group monthly climatology that is the dominant
cost.

**`_as_range` keeps the reduction's inner loop affine.** Today's
`view(A, eachindex(dim), :, :)` and the plan's `view(A, 1:n, :, :)` are
equivalent (both `IndexCartesian`); the win is against `Vector{Int}` groups,
where the parent access becomes a gather. NaNStatistics attaches `@simd ivdep`
to the *lowest-numbered* reduced dimension, so this matters most when the
grouped dimension is the fastest reduced axis.

**`getproperty` costs a folded `getfield`** with the `const` NamedTuple, and at
worst a field-name scan. Coordinates are built once per block and read once per
reduction, never per element. Do not call `g.lat` inside an inner loop.

**`deepcopy(var.dims)`** is `O(sum of dimension lengths)` — negligible for large
variables, but note `remake` deep-copies `attributes` and `dim_attributes` on top,
which dominates for a small variable such as a 1×1×120 time series.

**The weighted path is the most expensive operation in the design.**
`nanmean(A, W; dims)` expands to `sum(A.*W.*mask) ./ sum(W.*mask)` with a `Bool`
mask — three full-block temporaries, ~2.1× the block plus the mask. This is
pre-existing (`average_lonlat` does the same), but it means `LatWeightedAverage`
is not cheap, and a fine block grid multiplies the constant.

### Type stability — measured

Measured on Julia 1.12.6 with a standalone prototype of `_combine_blocks` exactly
as written above, reducing with `Reduce(mean)`. `Base.return_types` on the loop as
written gives **`Any`**, so the design has a real per-block dynamic dispatch.

Three *independent* causes, each verified by changing one thing at a time:

| Cause | Inferred type | Fix |
|---|---|---|
| `_block_coords` keys are runtime `Symbol`s | `_block_coords` → `NamedTuple`, so `Group{D, C}` is not concrete | store the names in a `Tuple{String,…}` **field** and have `getproperty` scan it |
| `Colon()` mixed with ranges in `src`/`dst` | `SubArray{Float64, 3, Array{Float64, 3}}` — parameters 4 and 5 unknown | `1:size(data, d)` instead of `Colon()` |
| `_as_range` applied per group | `groups` is `Vector{AbstractVector{Int64}}` as soon as one group has a gap | apply `_as_range` all-or-nothing per split operation |

With all three fixed the loop infers `Union{Nothing, Array{Float64, 3}}` — concrete
apart from the `ret = nothing` sentinel.

Two of those corrections matter for the plan as written:

- **`Colon()` is the culprit, not the `Nothing` sentinel.** `slot`'s
  `Union{Nothing, Int}` elements are harmless — inference union-splits them, and
  `1:size(data, d)` with `isnothing(slot[d])` still infers a fully concrete
  `SubArray`. `axes(var.data, d)` is **not** sufficient either: it yields
  `Base.OneTo{Int}`, a different type from the `UnitRange{Int}` of a group, so the
  eight branch combinations again exceed union-splitting. Only a literal
  `1:size(data, d)` makes both branches the same type.
- **A single non-consecutive group defeats everything else.** With
  `groups::Vector{AbstractVector{Int64}}` the loop infers `Any` even with the other
  two fixes applied. `_as_range`'s value is not the affine inner loop the plan
  already claims — it is this.

Cost, `plan` vs. all three fixed, best of several runs:

| Case | as written | fixed | ratio | allocated |
|---|---|---|---|---|
| one block, 360×180×12 | 92 µs | 79 µs | 1.16× | |
| 120 blocks of 360×180×1 (monthly climatology) | 54.0 ms | 46.9 ms | 1.15× | |
| 180 blocks of 360×1×12 | 3.06 ms | 2.31 ms | 1.32× | |
| 180 blocks of 360×1×120 | 45.8 ms | 48.6 ms | 0.94× | |
| **64800 blocks of 4×4×12** (0.25°→1° coarsening) | 326 ms | 47 ms | **6.9×** | 245 → 55 MiB |
| 64800 blocks of 1×1×12 | 271 ms | 32 ms | **8.4×** | 241 → 51 MiB |

So the instability is genuinely negligible when blocks are large — including the
120-group monthly climatology the plan names as its worst case — and costs 7-8× in
wall time and 4.4× in allocations when blocks are small and numerous. Conservative
coarsening by block-averaging is a real workflow that lands squarely in that
regime.

**This also qualifies "allocations are unchanged".** That holds for the
large-block cases. In the many-small-blocks regime the boxing from the dynamic
dispatch is itself 4.4× the necessary allocation, which no amount of `cat`
avoidance recovers.

**The two internal fixes buy 1.5× of the 8.4×; the rest needs the API change.**
Measured separately on the 64800×(1×1×12) case: as written 271 ms / 241 MiB;
ranges instead of `Colon()` 182 ms / 173 MiB; ranges *and* names-in-a-field
32 ms / 51 MiB. So replacing `Colon()` and fixing `_as_range` are free wins with
no visible consequence, while the remaining 5.7× requires giving up the
`coords` NamedTuple for a `names::Tuple` field plus a `coords::Tuple`.

That third one is smaller than it sounds. `g.lat === g.latitude` still works —
`getproperty` scans a tuple of at most `ndims` names instead of doing a NamedTuple
lookup, once per reduction rather than per element — and `values(g.coords)` becomes
just `g.coords`, which shortens the "naming a dimension at all" example above to
`prod(_cell_width.(g.coords))`. What is actually lost is `keys(g.coords)` returning
`Symbol`s, replaced by `g.names` returning `String`s, and NamedTuple destructuring
of `g.coords`. It still changes an exported type's shape and its docstring, so it
is a decision rather than a bug fix — just a cheaper one than the ratio suggests.

**Future work:** an `_apply!(op, out_view, block, group)` seam would let
`nanmean!`/`sum!` write straight into `view(ret, dst...)`, removing the last
per-block temporary. It is the only way the "allocations go down" claim becomes
true. Deferred — revisit if profiling justifies it.

---

## Files to touch

| File | Change |
|---|---|
| [src/split_apply_combine.jl](src/split_apply_combine.jl) | Full rewrite |
| [test/test_split_apply_combine.jl](test/test_split_apply_combine.jl) | Two amendments plus ~10 new testsets |
| [docs/src/split_apply_combine.md](docs/src/split_apply_combine.md) | New sections, **and** fix the existing prose at lines 95-98, which says the coordinate is "the first element of that group" |
| [docs/src/api.md:143-153](docs/src/api.md#L143-L153) | Add every new export, plus `Base.show(io::IO, sav::SplitApplyVar)` alongside the existing `Base.show(io::IO, var::OutputVar)` entry |
| [NEWS.md:3](NEWS.md#L3) | Entry under `main` |

Export list: `AbstractSplitOperation`, `AbstractApplyOperation`,
`SplitApplyVar`, `Group`, `GroupBy`, `GroupAll`, `SplitSeason`, `Bins`,
`Reduce`, `WithCoords`, `ReduceWithCoords`, `Average`, `LatWeightedAverage`,
`combine`.

`Group` must be exported, not merely documented — `@docs ClimaAnalysis.Group`
fails for an unexported name, and it is the type every user lambda receives.

`src/temp.jl` and `thoughts.md` (untracked prototypes) are superseded.

### Infrastructure traps

- `checkdocs = :exports` ([docs/make.jl:32](docs/make.jl#L32)) means every export
  needs a docstring **and** an `@docs` entry.
- `GroupAll`/`SplitSeason` become functions, so their existing docstrings must
  move onto the function definitions or they vanish, and `api.md:149-150` breaks.
- `DocTestSetup` does not import `Statistics`, so a `jldoctest` writing
  `Reduce(mean)` fails with `UndefVarError`. The existing docstrings use
  non-executing ```julia fences; keep that.
- Drop the `ClimaAnalysis.` self-qualification used at
  [lines 205](src/split_apply_combine.jl#L205) and
  [254](src/split_apply_combine.jl#L254) — `test/aqua.jl` asserts
  `check_no_self_qualified_accesses`.
- `Average` is a generic name in a reexported namespace; note the potential
  clash with other statistics packages in NEWS.

---

## Behavior changes to record in NEWS

1. Repeated dates no longer error under `SplitSeason`. The `allunique` check
   ([Utils.jl:372](src/Utils.jl#L372)) existed because the old code used
   `indexin`, which collapses duplicates; positional labeling is immune.
2. Group ordering, and each group's coordinate, come from the extreme coordinate
   rather than the first index. Identical for ascending dimensions; on a
   descending dimension `GroupAll("lat")` now returns `90.0`, not `-90.0`.
3. Chaining a second split operation composes instead of erroring.
4. `dims` reaches the reduction as a `Tuple`, not an `Int`. Only a hand-written
   reduction branching on `dims == 1` needs updating.
5. `GroupAll` on a zero-length dimension now errors instead of producing one
   empty group and a `NaN`.
6. `GroupAll`/`SplitSeason` are functions rather than types, so dispatching on
   them no longer works.
7. `LatWeightedAverage()` preserves the data's floating-point precision, so on a
   `Float32` variable it returns `Float32` where `weighted_average_lonlat`
   returns `Float64`. Also note `Average` is a generic name in a reexported
   namespace, and so are `Bins` (clashes with `DimensionalData.Bins`) and the
   already-shipped `combine` (clashes with `DataFrames.combine`).

---

## Verification

Per the repo workflow, use a tmux-hosted REPL, never a one-off process:

```
tmux new-session -d -s julia 'julia --project=. --history-file=no'
```

then in the REPL:

```julia
include("test/test_split_apply_combine.jl")
include("test/format.jl")
include("test/aqua.jl")
include("test/runtests.jl")
```

### Existing tests: two break

- **[lines 147-150](test/test_split_apply_combine.jl#L147-L150)** repeated-dates
  `@test_throws` — nothing throws now. Replace with a positive assertion that the
  duplicates group together.
- **[lines 173-179](test/test_split_apply_combine.jl#L173-L179)**
  `"A split operation is already set"` — chaining composes. Replace with
  `@test_throws r"already grouped"` on the *same* dimension, and add a positive
  test that two different dimensions compose.

Everything else passes. Three assertions change meaning silently and should be
tightened: [line 44](test/test_split_apply_combine.jl#L44) (`first` vs
`minimum` — identical on the ascending test dimensions),
[lines 50-53](test/test_split_apply_combine.jl#L50-L53) and
[159-162](test/test_split_apply_combine.jl#L159-L162) (errors now raised at the
pipe stage rather than in `combine`).

> **Floating-point caveat.** [Line 136](test/test_split_apply_combine.jl#L136)
> asserts `==` between a shuffled and a sorted variable. The plan buckets in
> positional order, so summation order within a group differs. It survives only
> because `Template.initialize` defaults to `reshape(1:n, …)` — exactly
> representable integers. Any rewrite that switches to non-integral data must use
> `≈`.

### New testsets

Cover, at minimum: `GroupBy` labeling and `nothing`-dropping; `Dates.year` and
`Dates.month`; `Bins` edge cases including the two constructor guards; block
reduction over two dimensions **and** that chaining order does not matter;
`LatWeightedAverage` vs `weighted_average_lonlat` with NaNs and with
`ignore_nan = false`; `SplitSeason() |> Average()` vs
`average_season_across_time` on data and dims; descending dimensions staying
descending; degenerate variables (1-D, single-element groups where `std` returns
`NaN` rather than throwing, duplicate coordinates, zero-length dimensions);
element types (`Int` stays `Int` under `sum`, becomes `Float64` under `mean`;
`Float32` data and coordinates, **and** the mixed `Float32` data with `Float64`
dimensions case, which is what NetCDF actually produces — assert the result is
`Float32`); the `Group` contract (`g.lat === g.latitude`,
`values(g.coords)` reachable, `g.lon` errors when longitude was not grouped,
`size(g.lat)` reshaped correctly, `dims` arriving as a `Tuple`); and `show` on a
`SplitApplyVar` both with and without an apply operation, asserting the
**resolved** dimension name appears when the user typed an alias.

Also assert dimension-order independence via `permutedims` and dimension-name
independence via a `latitude`/`longitude` variable.

### Performance check

`@allocated` on `GroupAll("time") |> Reduce(mean)` will show **no change** and is
the wrong acceptance criterion. Instead measure wall time and GC time on a
120-group monthly climatology (exposes `cat`'s splat) and on a fine block grid
such as 0.25°→1° coarsening (exposes the per-block constant, where the measured
penalty is 7×).

Do **not** assert that `combine` infers a concrete `ret`. As written it infers
`Any`, and the three causes are catalogued under "Type stability — measured";
only the third of them is fixable without changing `Group`. If the two internal
fixes are taken, the useful regression test is
`only(Base.return_types(_combine_blocks, …)) !== Any` — which still fails until
the NamedTuple goes, so gate it on that decision rather than writing it now.

---

## Corrections from review

Claims in earlier drafts that were checked and refuted, recorded so they are not
reintroduced:

- "Allocations go down." Total bytes are identical; only peak memory and `cat`
  overhead improve.
- "`getproperty` constant-folds with a literal symbol." Not with a `Dict`; hence
  the `const` NamedTuple.
- "`_make_regridder` rejects non-increasing dimensions." It `@warn`s and returns
  `nothing` ([Var.jl:146-151](src/Var.jl#L146-L151)), and bails for length-1
  dimensions anyway. The ordering-independence argument for the extreme-coordinate
  rule stands on its own; the monotonicity argument does not.
- "`g.lati` works." It does not; `LATITUDE_NAMES` has only two entries.
- "`reshape` on a non-contiguous view errors." It does not — 1.05× a materialized
  array. Materializing is still marginally better.
- "`combine` errors when a reduction drops the reduced dimension." `setindex!`
  silently accepts it; hence the explicit shape check.
- "The function barrier makes the `Group`'s `C` parameter concrete." Measured on
  Julia 1.12.6: `_block_coords` infers `NamedTuple` and `_combine_blocks` infers
  `Any`. A barrier fixes its *arguments*, but `_block_coords` computes its keys
  from runtime strings inside the barrier. See "Type stability — measured".
- "`slot`'s `Union{Nothing, Int}` costs inference." It does not — inference
  union-splits it. The `Colon()`/range mixing in `src` is what widens the view
  type, and `axes(var.data, d)` does not fix it either.

---
---

# APPENDIX — review suggestions, NOT incorporated

Everything below this line is **raw reviewer output for later evaluation**. None
of it is folded into the plan above. Reviews 1-4 (Julia correctness, API design,
tests/regressions, performance) *are* already incorporated in the body; reviews
5-7 below are not. Two further reviews (downstream integration, docs/migration)
were stopped before reporting.

## Decisions already taken during review

| Question | Decision |
|---|---|
| Simplicity vs. capability | **Keep the features, fix the bugs.** Do not adopt the ~305-line minimal formulation |
| Group coordinates from label or data | **Keep data-derived coordinates.** Document the caveat; add a test pinning the behavior |
| Group labels in output | **Deferred**, shape still open |
| `on` keyword | **Keep** |
| `g.dim_indices` naming | **Keep** (not `g.dims`, which collides with `OutputVar.dims`) |
| `ReduceWithCoords` + `WithCoords` | **Keep both** — alias plus wrapper |
| `dims` to the reduction | **Always a `Tuple`** |
| In-place `_apply!` seam | **Skip now**, noted as future work |
| `ignore_nan = false` weighted formula | **Match `average_lonlat` exactly** |
| `Average` name | **Keep**, note the collision risk in NEWS |
| What a coordinate-aware reduction receives | **Unresolved** — review stopped here |

## A. Bugs that must be fixed before implementing

These are defects in the plan body, not suggestions. Reported by the simplicity
and usability reviewers after the body was written.

1. **`on` breaks `GroupAll` dimension-name agnosticism.**
   `GroupAll(dim) = GroupBy(dim, Returns(:all), var -> var.dims[dim])` captures
   the *unresolved* name. `_resolve` rewrites `dim_name` but carries `op.on`
   through verbatim, so `GroupAll("t")` on a variable whose dimension is `"time"`
   throws a `KeyError` inside `combine`. Fix: have `on` receive the resolved name,
   or rebuild the closure in `_resolve`, or guard the date conversion on
   `eltype(coord) <: AbstractFloat` instead of using `on` at all.

2. **`Average(ignore_nan = false, weights = ...)` does not run.**
   `mapslices(f, A; dims)` passes a slice with the non-reduced dimensions
   **dropped** — `(181, 361)` for a `(120, 181, 361)` block reduced over `(2, 3)`
   — while `w` is `(1, 181, 1)`, so the broadcast throws. `average_lonlat` works
   only because it calls `dropdims(lat_weights; dims = dims_to_drop)` first
   ([Var.jl:975-987](src/Var.jl#L975-L987)). That step was lost in compression.

3. **`_out_size(::Transform, ...)` references a type that does not exist.**
   Written literally it does not compile. It is illustrative only — mark it as
   such or delete it.

4. **`_cell_width` does not exist.** Zero occurrences in the repository. The
   "Naming a dimension at all" example is unrunnable. If the intent was
   `Numerics._integration_weights_generic_left`, that is private and appends a
   **trailing zero**, so used as an averaging weight it silently discards the
   last latitude band.

5. **Dead clause in `Bins`.** After the `i == length(edges)` line, and given
   `searchsortedlast` returns `0..length(edges)`, the `i > length(b.edges) - 1`
   half of the next test is unreachable.

6. **The radians warning must not live inside the reduction.**
   `average_lonlat` warns once per variable ([Var.jl:963-964](src/Var.jl#L963-L964));
   inside `LatWeightedAverage` it would run per block and warn spuriously for a
   thin equatorial `Bins` band. Hoist it to pipe time.

## B. Simplicity review

**Line count**, formatted at margin 80, one docstring per export:

| Component | code | doc | total |
|---|---:|---:|---:|
| exports (14 names) | 13 | 0 | 13 |
| two abstract types | 2 | 12 | 14 |
| `GroupBy` + kw constructor | 9 | 24 | 33 |
| `GroupAll`, `SplitSeason` | 2 | 22 | 24 |
| `_labeling_values` | 11 | 6 | 17 |
| `Bins` + call | 15 | 18 | 33 |
| `_group_indices` + `_as_range` | 32 | 12 | 44 |
| `Reduce`, `WithCoords`, `_call`×2, `_apply`, `ReduceWithCoords` | 16 | 46 | 62 |
| `Group`, `getproperty`, `propertynames`, `_CONVENTIONAL`, `_conventional` | 33 | 28 | 61 |
| `Average`, `LatWeightedAverage` | 23 | 34 | 57 |
| `SplitApplyVar` + pipes + `_resolve`×2 | 27 | 26 | 53 |
| `combine`, `_combine_blocks`, `_block_coords` | 47 | 22 | 69 |
| **Total** | **230** | **250** | **~480** |

Against 255 today (111 code / 144 doc): **1.9× the file, 2.1× the code.**
The dominant term is the export count — 14 exports × ~18 doc lines is 52% of the
file. Fewer lines means cutting names, not statements.

**Cut list, ranked by value-per-line removed.** Not adopted; recorded for later.

1. `Group` + `getproperty` + `propertynames` + `_CONVENTIONAL` + `_conventional`
   — **~61 lines**. `_block_coords` already keys by
   `Symbol(conventional_dim_name(...))`, so the alias apparatus buys exactly one
   thing: typing `g.lat` instead of `g.latitude`. A plain NamedTuple
   `(; dims = (2, 3), latitude = ..., longitude = ...)` gives `g.dims`,
   `g.latitude`, `values` and destructuring with zero code.
2. `WithCoords` + `_call`×2 + the alias — **~22 lines, one export**. Justified by
   `Transform(WithCoords(f))` for a type the plan explicitly defers. Two concrete
   `_apply` methods are five lines.
3. The `on` keyword — **~14 lines**, and cutting it removes bug A.1. Its only
   caller is the `GroupAll` workaround; the ENSO use case has no caller anywhere.
4. `Average`'s `weights` keyword and the `ignore_nan = false` branch —
   **~14 lines**, and cutting it removes bug A.2. One caller
   (`LatWeightedAverage`); a user wanting custom weights writes one line via
   `ReduceWithCoords`.
5. `_out_size` seam — **2 lines**, does not compile (bug A.3).
6. `Bins` as a struct — **~10 lines**. The plan concedes this worsens PR concern
   #2 and points at [Template.jl:139-141](src/Template.jl#L139-L141) as the
   repo's answer: return a closure. Both constructor guards stay either way.
7. Both abstract types — **~14 lines, two exports** (medium confidence).
   `GroupBy` becomes the only split operation.
8. `Average` itself — marginal. Real value: `nanmean` is imported into `Var` but
   not re-exported, so `Reduce(nanmean)` needs a user-side `import NaNStatistics`.
9. `rev`/`pick` orientation preservation — 2 lines, marginal either way.

**Do not cut**: the explicit shape check (best value-per-line in the plan);
`sortperm` over precomputed representatives (labels are not orderable —
`("DJF", 2010)` sorts DJF < JJA < MAM < SON, so sorting by representative
coordinate is the *only* way to satisfy the shuffled-dates test); the concretely
typed bucket dict; `Val(N)` and the `_combine_blocks` barrier; `_as_range`;
`similar(block, sizes)`; the date conversion in `_labeling_values`.

**Minimal formulation**: ~305 total (150 code / 155 doc), 9-11 exports. Gives up
`g.lat`, the `on` escape hatch, `LatWeightedAverage(ignore_nan = false)`, and
`Bins` as a dispatchable type. The reviewer's blunt verdict: 50% is not
reachable; ~305 is the floor while keeping multi-dimensional grouping,
coordinate-aware reductions, and one docstring per export.

**Comment density.** One outright violation in the plan's code:
`# Attributes and dimension attributes are passed through, as they are today` is
diff-narration describing code that is not there. The larger risk is the prose —
"Three details are load-bearing", "Details that are load-bearing", "Three things
here are not stylistic" — which an implementer will transcribe into `src/` as 60
lines of justification comments. **That rationale belongs in the PR description,
not in the source.**

## C. Climate-scientist usability review

Workflows in the proposed API, `var` with dims `(lon, lat, time)`:

| Workflow | Verdict |
|---|---|
| Climatological monthly means | Works mechanically; output silently misordered and unlabeled |
| Annual means | Works for daily data; subtly wrong for monthly data (no day-length weighting) |
| DJF-vs-JJA composite | Works, and is the best thing in the plan — but carries no labels |
| Global area-weighted mean series | Works, bitwise-equal to `weighted_average_lonlat`; more typing than the existing one-liner, loses the `long_name` update |
| Zonal mean | Works; identical to `average_lon`, no advantage |
| Lat/lon box mean | Works and is elegant; no user will discover it; output coordinate is the box corner |
| Land-only / ocean-only mean | Works and is correct. Fractional land-area weighting impossible |
| Anomaly vs. climatology | **Impossible** — needs `Transform`, and no post-hoc workaround because `arecompatible` requires equal dim arrays ([Var.jl:1429-1447](src/Var.jl#L1429-L1447)) |
| ENSO composite | Works; `on` is the right seam |
| Percentiles per season | Works, but NaNStatistics is not re-exported |
| Diurnal cycle | Works, same ordering issue; local-solar-time impossible |
| monthly→annual, daily→monthly | Works; undiscoverable — nobody searches docs for "GroupBy", they search for "resample" |

**The ordering finding** (highest severity in that review, superseded by the
decision to keep data-derived coordinates, recorded anyway): groups are ordered
by the earliest coordinate in each group. For a run starting in **July** —
spin-up trimmed, water-year, restart — a monthly climatology comes out ordered
7,8,9,10,11,12,1,2,3,4,5,6 with a time axis reading Jul 2010 … Jun 2011, and is
bit-for-bit indistinguishable from a raw 12-month window. The reviewer notes the
shuffled-dates test only requires coordinate ordering for *coordinate-valued*
groupings; for label-valued groupings (month, hour, season, bin index) sorting by
label satisfies it equally well.

**The `g.time` asymmetry.** `_labeling_values` converts time to `DateTime` before
handing it to `by`, but `_block_coords` builds `g.time` from `var.dims[...]` —
**raw seconds**. So `by` sees dates and `weights` sees seconds, for the same
dimension in the same pipeline. This is why month-length weighting takes three
lines of plumbing instead of one keyword.

**xarray comparison** — where the plan loses: `groupby("time.month")` renames the
output dimension to `month = 1…12`; `groupby_bins` labels bins as intervals;
`resample` is the word this audience types; `.weighted(w).mean(dim=("lat","lon"))`
needs one concept where the plan needs three (you must *group* a dimension to
reduce it, weights are a function of a `Group`, weighting an ungrouped dimension
errors). Where the plan wins: per-year seasons via `SplitSeason` (fiddly in
xarray), arbitrary `by` lambdas beat the `"time.xxx"` string DSL, and the joint
weighted lon-lat reduction is genuinely better than xarray's average-of-averages.

**Gratuitous divergences**: `combine` as a mandatory fourth stage; one `GroupAll`
per dimension where xarray writes `mean(dim=("lat","lon"))` — fixable with
varargs `GroupAll("lon", "lat")`; and `Average()` ignoring NaN while
`Reduce(mean)` does not, with the tutorial teaching the NaN-unsafe spelling
([docs/src/split_apply_combine.md:64-68](docs/src/split_apply_combine.md#L64-L68)).

**On `cosd` weighting**: it is exactly right, not an approximation, for uniform
latitude spacing — band area ∝ `2 sin(Δ/2) cos φ`. It degrades on Gaussian grids,
on `resampled_as` output onto non-uniform targets, and where poles are cell
centres (`cosd(±90) = 0` zeroes a half-cell). Note the package already disagrees
with itself: `integrate_lonlat` uses `Δφ·cos φ`
([Numerics.jl:125-132](src/Numerics.jl#L125-L132)) while
`average_lonlat(weighted = true)` uses bare `cos φ`, and neither docstring says
"assumes equispaced latitudes".

**Ranked gaps**: group labels in output; order by label; `Transform`; varargs
`GroupAll`; `weights` accepting an array or `OutputVar` aligned by dimension name;
dates on the group; named groupers (`ByYear`, `ByMonth`, `MonthlyClimatology`);
`combine` updating `long_name` via `_update_long_name_generic!`; `Bins` output
coordinate as midpoint; documented recipes for `by → nothing` as regional
selection and `on = eachindex` as coarsening; auto-drop or document `dropdims`;
rolling windows.

**Discoverability walkthrough** (new user attempts a lat-weighted seasonal mean):
the docs page teaches only `GroupAll` + `Reduce`; the obvious
`SplitSeason() |> LatWeightedAverage()` errors; they must learn "to weight a
dimension you must group it", which exists nowhere in the xarray/pandas
vocabulary; `GroupAll("lat")` alone silently gives a different (zonal-band)
answer; the result is `(1, 1, n)` and will not plot without `dropdims`; and they
cannot tell which season is which. Estimated 20-40 minutes with the API index
open, or a fallback to `weighted_average_lat(average_season_across_time(var))`,
which is an average of averages and slightly wrong.

## D. Alternative-design review

**Concern 1 — data-dependent coordinates break the comparison workflow.**
`arecompatible` ends in exact equality of dim arrays
([Var.jl:1446](src/Var.jl#L1446)) and gates every binary op
([outvar_operators.jl:15](src/outvar_operators.jl#L15)). Bin a model variable and
an obs variable with the same `Bins(-90:30:90)`; if one grid starts at −90.0 and
the other at −89.5, `sim - obs`, `bias` and `global_rmse` error out. Same for
`GroupBy("time", Dates.year)` on runs starting Jan 1 vs Jan 15. **Decision taken:
keep data-derived coordinates.** Recorded because it is the reviewer's
highest-severity finding and the mitigation is localised — a
`_representative(labeler, label, coords)` hook, ~6 lines, with today's `minimum`
as the generic fallback.

**Concern 2 — `weights = g -> ...` cannot express what matters.** `Group` carries
no attributes, so month-length weighting is impossible (you cannot get
`start_date` inside the lambda), as are area weighting from `areacella` and land-
fraction weighting from `masks.jl`. And `weights(g)` runs once per block for a
quantity that usually depends only on the variable.

**Concern 3 — `SplitApplyVar` is lazy in name only.** `combine` calls
`_group_indices` fresh for every op; nothing is cached or deferred. It costs two
abstract types, three callable methods, two "already set" error paths, `_resolve`
plus its fallback, four exports, four `@docs` entries and two testsets, and buys
one capability nothing uses.

**Concern 4** — 14 exports for a feature expressible in 4-5.
**Concern 5** — `Bins(edges)` collides with `DimensionalData.Bins(f, bins)`; a
Rasters.jl user writing `Bins(month, 12)` gets a `MethodError`.
**Concern 6** — the plan's stated payoff is false as scoped: with labels deferred
and attributes untouched, `SplitSeason() |> Average() |> combine` is strictly
*less* informative than `average_season_across_time`, which writes `"season"`,
`"year"` and calls `_update_long_name_generic!`.
**Concern 7** — `f(block; dims)` excludes `Reduce(x -> quantile(x, 0.9))` and any
`f(::AbstractVector)::Number`, failing with a confusing `dims` `MethodError`.
Suggests a three-line `ReduceSlices(f)` on its own merits.
**Concern 8** — the date-vs-seconds switch is a silent heuristic keyed on whether
`start_date` happens to be present.

### Alternative A — one function, `dim => labeler` pairs, no pipeline objects

```julia
groupreduce(f, var::OutputVar, groups::Pair...)
groupreduce(f, groups::Pair...) = var -> groupreduce(f, var, groups...)  # for |>
```

| plan | Alternative A |
|---|---|
| `var \|> GroupAll("time") \|> Reduce(mean) \|> combine` | `var \|> groupreduce(mean, "time" => All())` |
| `var \|> SplitSeason() \|> Average() \|> combine` | `var \|> groupreduce(nanmean, "time" => season)` |
| `GroupBy("time", identity, on = var -> enso)` | `"time" => enso` (a plain vector) |

Keeps `_group_indices`, `_combine_blocks`, `Bins` and the reduction contract
verbatim. Deletes `SplitApplyVar`, both abstract types, all three callable
methods, `_resolve` and its fallback, both "already set" errors, `combine`, and
`on`. Accounting: **−70 code lines, −60 doc lines, −8 exports, −8 `@docs`
entries, −2 testsets.** Costs: no partial application; no user-defined split
subtypes; `groupreduce` collides with `SplitApplyCombine.groupreduce`.

### Alternative B — inspectable groups (the DimensionalData / YAXArrays shape)

`groups(var, "time" => season)` returns a `GroupedVar <: AbstractArray{OutputVar}`
you can `show`, index and plot **before** applying anything, then
`combine(nanmean, gv)`. Biggest win: the split becomes inspectable, which is the
single largest debuggability gain available and which the plan gives up entirely.
Also dissolves the group-metadata question structurally. Not recommended as
primary because it calls the reduction once per group rather than once per block,
losing the plan's vectorization property. Steal the `show` and labels-as-coords.

### Alternative C — pass a sub-`OutputVar` instead of a `Group`

```julia
function broadcast_dim(var::OutputVar, dim_name, values)
    idx = var.dim2index[find_corresponding_dim_name_in_var(dim_name, var)]
    return reshape(values, ntuple(d -> d == idx ? length(values) : 1, ndims(var.data)))
end

LatWeightedAverage(; ignore_nan = true) = Reduce(WithVar() do v; dims
    nanmean(v.data, broadcast_dim(v, "lat", cosd.(latitudes(v))); dims = dims)
end)
```

Deletes `Group`, its `getproperty` and `propertynames`, `_CONVENTIONAL`,
`_conventional`, the conventional-key rule, the reshape in `_block_coords`, one
export, one `@docs` entry and one docstring — **~45 lines and one public type** —
in favour of a 4-line helper that is independently useful
([Var.jl:965-971](src/Var.jl#L965-L971) inlines exactly this reshape today).
Reductions can then reach `dates(v)`, `v.attributes`, units and masks, which is
what makes concern 2 solvable. Costs: two small dicts allocated per block; loses
the "not a grouped dimension" guard. **This was the open question when the review
was stopped.**

### Alternative D — weights resolved once, from the variable

`weights :: Nothing | AbstractArray | OutputVar | Function(var) -> array`,
evaluated **once** before splitting rather than per block. Makes cell-area and
land-fraction weighting one-liners and month-length weighting expressible.
Composes with Alternative C rather than competing.

### Ecosystem verdict: do not depend, do borrow

DimensionalData ships this feature (`groupby(A, Ti => month)`, `Bins(f, bins)`,
`DimGroupByArray`), but `OutputVar` is not an `AbstractDimArray`, so using it
means round-tripping and losing `attributes`/`dim_attributes`, and DD is 0.x with
frequent breaking releases. This repo has nine dependencies and pins NaNStatistics
to a hand-curated allowlist of patch versions — a team that has priced dependency
churn. Taking DD for a 250-line feature is a bad trade. SplitApplyCombine.jl is
the wrong shape entirely (no concept of reducing along an axis while keeping it).
Borrow for free: the `dim => labeler` pair syntax, the `groupby`/`groupreduce`
vocabulary, `Bins`'s function-first signature or a non-colliding name.

A weak-dep `ClimaAnalysisDimensionalDataExt` providing `DimArray(var)` /
`OutputVar(::AbstractDimArray)` is where a DD dependency should be spent if ever
— it opens the Rasters/YAXArrays ecosystem without owning any of DD's API.

### What that review would reopen

Drop `combine` (make apply eager); delete `GroupAll`/`SplitSeason` now that
`GroupBy` exists; bring labels and `long_name` into scope; reconsider per-slice
reductions as `ReduceSlices(f)`; rename `Bins` or adopt DD's signature; decide
`drop`/`dropdims` now rather than adding a second API later
([thoughts.md:6](thoughts.md#L6) wanted it, and every other `average_*` squeezes).

**Smallest change capturing most of the value, if only one thing moves:** make
the output coordinate a function of the label rather than the data. Second:
Alternative C. Third: Alternative A. All three are independent.

## E. Reviews not completed

- **Downstream integration** — masks, `FlatVar`, Leaderboard/calibration, units,
  Makie recipes, `reverse_dim`/`permutedims`/`resampled_as`, and which existing
  functions could be reimplemented on the new API. Stopped before reporting. The
  one question worth answering from it: after a time average the grouped dimension
  stays at size 1, so **can the result be plotted directly**, or is a `dropdims`
  step needed? The usability review suggests the latter.
- **Docs and migration** — actual docstrings for every export, the actual docs
  page, the NEWS entry, an old→new migration table with a column for where
  semantics differ, and a concrete answer to PR concern #2 (callable-struct
  docstrings). Stopped before reporting.
