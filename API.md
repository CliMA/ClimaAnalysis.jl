# Grid/Dim API reference

## Types

```julia
struct Dim{A <: AbstractArray, B, T}
    name::String                      # as in the source file, e.g. "lon"
    values::A                         # live; usually 1D, N-D legal (topography z)
    units::Base.RefValue{String}      # "" = missing; promoted first-tier
    bounds::B                         # cell bounds (CF-style n×2) or nothing; promoted first-tier
    attributes::Dict{String, T}       # everything else (Any-valued OK)
end

abstract type AbstractGrid end

struct Grid{D <: Tuple{Vararg{Dim}}} <: AbstractGrid
    dims::D                      # order == data axis order; names unique (checked)
end
```

`units` and `bounds` are promoted out of `attributes` as first-tier
fields. `units` is held in a `Ref` so the struct stays immutable while
`set_dim_units!` (public, in-place) keeps working — the same
interior-mutability pattern as the live `values` array and `attributes`
dict (design decision 4). `bounds` is a plain field (`nothing` when
absent): no public in-place setter exists, so filling bounds later is a
functional rebuild (new `Dim` via `replace_dim`), consistent with the
other structural edits. Constructors extract `"units"`/`"bounds"` from
attribute dicts (NetCDF reads land there) into the fields. Two
consequences: `Base.copy(grid)`/`copy(d::Dim)` must allocate a fresh
`units` `Ref` (sharing it would make `set_dim_units!` on a copy mutate
the original); and materialized `var.dim_attributes` merges the fields
back in as entries, so *reads* keep working, but *writing*
`var.dim_attributes[name]["units"]` through the materialized dict is
unsupported (documented; use `set_dim_units!`) — same class of break as
dict insertion.

Fields are otherwise Private: all consuming code goes through the
accessor/verb API (never `grid.dims` directly) so a Vector-backed
`AbstractGrid` sibling can drop in later.

Complete function list for the `Grid`/`Dim` abstraction designed in
[GRID_DESIGN.md](GRID_DESIGN.md), with signatures. Markers: **exported**,
**Base** (overload, no export needed), **internal**.

The "Used by / absorbs" column records the existing code each function
replaces or backs — i.e. why it earns its place. Entries marked
**nice to have** have no current absorber: they exist for API
completeness/ergonomics and could be dropped from the initial scope
without leaving a hole in Var.jl.

The "Visibility" column is the commitment level, kept to a bare minimum
because the design may still change: **Public** = documented, semver-stable
from v0.6.0; **Private** = implementation detail or experimental —
docstring marked internal, free to change/remove, promotable to Public
later once proven. The principle: Public covers only *construct, read,
subset, compare, display* on grids, plus the OutputVar-facing surface;
all verbs, predicates, and introspection start Private since their only
required consumers are internal (Var.jl etc.). This narrows the design
doc's export list — only Public non-Base names should be exported
(`Grid`, `Dim`, `AbstractGrid`, `dim_names`); Base overloads need no
export and are Public purely by documentation. (`public` keyword
unavailable while supporting Julia 1.10; publicness = docs inclusion.)

## Constructors

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `Dim(name::String, values::AbstractArray; units = "", bounds = nothing, attributes = Dict{String,Any}())` | exported | Public | every `Grid` construction site; `read_var`'s dim-building; extracts `"units"`/`"bounds"` if passed inside `attributes` |
| `Grid(dims::Dim...)` | exported | Public | primary form; `Grid()` empty grid needed by 0-dim reductions and `_split_along_dim` |
| `Grid(pairs::Pair{String,<:AbstractArray}...)` | exported | Private | nice to have — test/user ergonomics |
| `Grid(dims_dict::AbstractDict, dim_attributes::AbstractDict = ...)` | exported | Private | back-compat bridge: keeps 4-arg `OutputVar` constructor (Var.jl:241-271) and `remake` working unchanged |
| `check_grid_consistency(grid, data)` | internal | Private | length-vs-size check in `OutputVar` constructor (Var.jl:246-261) |

## Accessors (lookup)

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `grid[name::AbstractString]` | Base | Public | every `var.dims[name]` read; backs `times`/`longitudes`/... delegators |
| `dim(grid, name::AbstractString)::Dim` | exported | Private | verbs and interpolant internals that need values + attributes together |
| `dim(grid, i::Integer)::Dim` | exported | Private | positional access in indexing/`permutedims` internals |
| `dim_names(grid)` | exported | Public | discovery entry point for `grid[name]`; `keys(var.dims)` sites; name checks in `cat`, compatibility predicates |
| `dim_index(grid, name)` | internal | Private | absorbs `dim2index` (~20 sites, incl. both Makie extensions) |
| `Base.keys(grid)`, `Base.iterate(grid)`, `Base.length(grid)` | Base | Private | `getproperty` materialization of `var.dims`; iteration protocol (dims-vs-points ambiguity — keep changeable) |
| `Base.haskey(grid, name)` | Base | Private | exact-name checks (e.g. `read_var` guards) |
| `hasdim(grid, name)` | exported | Private | conventional-name checks; backs `has_time`/`has_longitude`/... (var versions stay the public surface) |
| `Base.size(grid)` | Base | Public | data-shape source of truth: constructor check, `cat`, `ones`, flat.jl |
| `Base.ndims(grid)`, `Base.isempty(grid)` | Base | Public | `Base.isempty(var)` (Var.jl:771-778); dimension-count checks |
| `is_cartesian(grid)` | exported | Private | generalizes `is_z_1D` (Atmos, interpolation guards); concept still experimental |
| `has_time(grid)`, `has_longitude(grid)`, ... | exported | Private | var versions in outvar_dimensions.jl become one-line delegators and remain the public API |
| `time_name(grid)`, `longitude_name(grid)`, ... | exported | Private | same |
| `times(grid)`, `longitudes(grid)`, `latitudes(grid)`, ... | exported | Private | same (`dates(var)` stays var-level, needs `start_date`) |

## Dim introspection

Var-level `dim_units`/`set_dim_units!`/`range_dim` remain the public API;
everything here is grid-internal plumbing until `Dim` proves stable.

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `units(d::Dim)` | exported | Private | `d.units[]`; backs `dim_units` (accessors hide the `Ref`) |
| `set_units!(d::Dim, u)` | exported | Private | `d.units[] = u`; backs `set_dim_units!` (Var.jl:710-722, the dict-insertion hazard) |
| `bounds(d::Dim)` | exported | Private | field accessor; `nothing` when absent; consumed by `cell_widths` when present (no setter — rebuild the `Dim` to change bounds) |
| `convert_units(d::Dim, new_units; conversion_function)` | exported (new method of existing var-level export) | Private | the deepcopy-dims-and-refill dance in `convert_dim_units` (Var.jl:657-668) becomes `replace_dim(grid, name => convert_units(...))` |
| `Base.extrema(d)` | Base | Private | backs `range_dim` |
| `Base.length(d)`, `Base.eltype(d)`, `Base.ndims(d)` | Base | Private | interpolant construction (`_make_interpolant`, Var.jl:139-150) |
| `Base.first(d)`, `Base.last(d)`, `Base.issorted(d)` | Base | Private | `_find_extp_bound_cond` (Var.jl:216-239); sorted-dim guard |
| `isequispaced(d)` | exported | Private | wraps `Utils._isequispaced`; interpolant + `resampled_as` paths |
| `Base.step(d)` | Base | Private | equispaced dims only; extra-longitude-point logic (Var.jl:178-209) |

## Collection semantics

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `Base.copy(grid)` | Base | Private | `Base.copy(var)` (Var.jl:780-785) |
| `Base.deepcopy(grid)` | Base | Private | callers' `deepcopy` habits via `remake` |
| `Base.empty(grid)` | Base | Private | `empty(var.dims)` idiom (Var.jl:1674, 1705, 2157) |
| `Base.map(f, grid)` | Base | Private | nice to have — rebuild-all-dims idiom; generic `map` + splat covers it (see design tension below) |
| `Base.filter(pred, grid)` | Base | Private | subset-dim comparison in `check_dims_consistent`; internal building block for indexing/`dropdims`/`spatial_dims` |
| `coordinates(grid)` | exported | Private | nice to have — explicit Cartesian-product-of-points iteration (`Iterators.product`) |
| `flattened_length(grid)` | internal | Private | hand-rolled `prod(length.(values(dims)))` (flat.jl:246, 253, 312) |

## Structural indexing

Integer drops the dim, range/vector keeps it. Subset arrays are copies;
N-D dims accept only `:`. Public — this is user requirement #4 and a
headline feature, with semantics locked to Base arrays.

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `grid[i, j, ...]` | Base | Public | `window`, `_select` (outvar_selectors.jl:450-459), `_split_along_dim` |
| `grid[time = 1:10, lat = 3]` | Base | Public | named form of the same; unmentioned dims = `:` |
| `Base.view(grid, inds...)` | Base | Private | nice to have — no current absorber; symmetry with array `view`, shared arrays instead of copies |
| `Base.view(grid; kwargs...)` | Base | Private | nice to have — named form of the above |
| `Base.maybeview(grid, args...; kwargs...)` | Base | Private | nice to have — only so `@views` blocks produce grid views instead of silent copies |

## Operation verbs (functional — each returns a new `Grid`)

All Private initially: their required consumers are the Phase 2 internal
migrations. Promote (and export) individually once the design settles.

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `Base.dropdims(grid; dims)` | Base | Private | `_reduce_over` pop-loop (Var.jl:836-842) → every reduction |
| `Base.dropdims(grid)` | Base | Private | drop all singleton dims (xarray `squeeze`); grid companion to `Utils.squeeze`, closes the TODO at Var.jl:830 |
| `Base.reverse(grid; dims)` | Base | Private | functional counterpart to `reverse_dim!`; today `reverse_dim` copies the whole var for non-mutating semantics; pairs with backlog `sort` |
| `replace_dim(grid, name => new_values)` | exported | Private | resample refills (Var.jl:1674-1715), `transform_dates`, `_dates_to_seconds`, `shift_longitude`, `average_season_across_time` |
| `replace_dim(grid, old_name => new_dim::Dim)` | exported | Private | `Atmos.to_pressure_coordinates` type-widening dance (Atmos.jl:88-102) |
| `rename_dim(grid, old => new)` | exported | Private | sugar over `replace_dim`; Atmos rename sites |
| `Base.permutedims(grid, perm)` | Base | Private | ~20-line dim_attributes reorder dance (Var.jl:1549-1564), `reordered_as`, flatten's perm |
| `permutation(grid, perm)::Vector{Int}` | internal | Private | data-reorder companion to `permutedims` |
| `Base.cat(grids...; dims::String)` | Base | Private | dim-checking/concat halves of `Base.cat(vars...)` (Var.jl:2770-2853) |
| `Base.ones([T,] grid)` | Base | Private | 3 identical `ones_var` blocks (Var.jl:2251, 2337, Atmos.jl:195) |
| `ones_like(var)` | exported | Private | var-level companion to `ones` |
| `reverse_dim!(grid, name)` | exported | Private | existing `reverse_dim`/`reverse_dim!` (var-level stays public); returns axis index |

## Geometry / weights

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `cell_widths(d::Dim)` | internal | Private | spacing re-derivation in integration (Numerics.jl via `integrate_lon/lat`); center-vs-edge longitude logic (Var.jl:178-209, 224-232) |
| `latitude_weights(grid)` | exported | Private | hand-rolled weights in `weighted_average_lat`/`weighted_average_lonlat` (Var.jl:866-898) |
| `coordinate_field(grid, name)` | exported | Private | reshape gymnastics in `masks._generate_binary_mask` (masks.jl:167-178) |

## Value-based selection

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `get_index(grid, name, value, selector::AbstractSelector)` | exported | Private | var-level `get_index` methods delegate and remain the public selector extension API (outvar_selectors.jl:253-296) |
| `window_indices(grid, name, left, right)` | internal | Private | `_get_window_indices` (outvar_selectors.jl:306) |

## Compatibility predicates

Var-level `arecompatible` remains the public API; grid-level predicates
start Private and get promoted if users need grid-to-grid comparison.

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `matching_dim_names(a::Grid, b::Grid; order)` | exported | Private | one of three consolidations of `arecompatible` (Var.jl:1415, fixes vacuous check at :1419) |
| `matching_dim_values(a::Grid, b::Grid; approx)` | exported | Private | same consolidation |
| `matching_dim_units(a::Grid, b::Grid; on_missing)` | exported | Private | same consolidation |
| `iscompatible(a::Grid, b::Grid)` | exported | Private | composite; replaces `arecompatible` |
| `check_dims_consistent(a::Grid, b::Grid; dims)` | exported | Private | throwing composite; replaces `_check_dims_consistent` (Var.jl:1446-1522) and flatten's checks (flat.jl:146-192) |

## Interpolation

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `make_interpolant(grid, data)` | internal | Private | Var.jl:131-239; callers at Var.jl:1398/1669/1717, masks.jl:308 |
| `extrapolation_bc(d::Dim)` | internal | Private | periodic-longitude / flat-latitude detection (same block) |

## Equality / hashing / display

| Signature | Marker | Visibility | Used by / absorbs |
|---|---|---|---|
| `Base.:(==)(a::Dim, b::Dim)` | Base | Public | name + values + units + bounds + attributes (`Ref`s compared by contents, not identity) |
| `Base.:(==)(a::Grid, b::Grid)` | Base | Public | pairwise, in order; `grid == OrderedDict` NOT defined (~50 test assertions go through materialized `var.dims`) |
| `Base.hash(::Dim)`, `Base.hash(::Grid)` | Base | Public | consistency requirement with `==` |
| `Base.isapprox(a::Grid, b::Grid; kwargs...)` | Base | Private | conventional-name pairing is a design choice that may change; absorbs `≈` checks in `cat` (Var.jl:2770-2853) and `flatten` (flat.jl:168) |
| `Base.show(io, ::Dim)`, `Base.show(io, ::Grid)` | Base | Public | extracted from `show(::OutputVar)` (Var.jl:2875-2937), which reuses them |
| `Base.summary(grid)`, `Base.summary(d::Dim)` | Base | Public | compact one-liner (`"Grid(time × lon × lat, 120×360×180)"`) for `show(::OutputVar)` headers and error messages |

## OutputVar-side additions

The actual user-facing surface — all Public.

| Signature | Visibility | Used by / absorbs |
|---|---|---|
| `OutputVar(attributes, grid::AbstractGrid, data)` | Public | new primary constructor; dict-based constructors kept, build a `Grid` internally |
| `var.grid` | Public | property accessor (no `grid(var)` function) |
| `var.dims`, `var.dim_attributes` | Public | back-compat: materialize genuine `OrderedDict`s via `getproperty` (live arrays); keeps ~39 rebuild sites compiling through Phase 1 |
| `remake(var; grid, ...)` | Public | gains `grid` kwarg; `dims`/`dim_attributes` kwargs stay accepted |

## Design tension noted

`Base.map(f, grid)::Grid` reads ambiguously — "map over dims" (the
dict-flavored protocol view, consistent with `iterate`) vs "map over the
Cartesian product of points" (the array-intuition view, spelled
`coordinates(grid)`). Base precedent: dicts have a dict-returning
`filter` but no dict-returning `map`. Left as designed for now; if it
confuses in practice, drop the overload (generic `map` + `Grid(...)`
splat covers the use case) or rename to `map_dims`. Kept Private partly
for this reason.

## Maybe add

Candidates with a plausible but not-yet-forced case — decide when a
caller appears. All would start Private.

| Signature | Rationale |
|---|---|
| `Base.issubset(a::Grid, b::Grid)` | "is `a` a subgrid of `b`" (names and coordinate ranges contained); the right guard before `resampled_as` extrapolation errors; sits between backlog's point-level `in` and grid-level `intersect` |
| `eachslice(grid, data; dims = "time")` | yields `(coordinate, data_slice)` pairs — the loop inside `_split_along_dim` and every "for each time step" workflow; subsumed by backlog `groupby` (it's `groupby` with `identity`), so only add if `groupby` doesn't land |
| `guess_bounds(d::Dim)` | Iris's name for estimating cell bounds from centers; returns a new `Dim` with the `bounds` field set; front door for the `cell_widths` heuristics now that `bounds` is first-tier |

## Open items (stubs in GRID_DESIGN.md, not yet specified)

- Grid-level `remake`
- `unflatten`

### N-D dims (topography z): semantics and containment

An N-D `Dim` is a category error kept for compatibility: a 3D
topography-following `z` is not a dimension but a *coordinate* laid over
several axes. Every N-D headache is that error surfacing — `length(d)`
means nothing useful, `first`/`issorted`/`step`/`isequispaced` are
undefined, indexing accepts only `:`, `size(grid)` requires
`is_cartesian`, and the constructor's shape check has an N-D escape
hatch (Var.jl:249).

Prior-art consensus (no mature library models N-D as a dim):

- **CF conventions**: dims are 1D coordinate variables; anything N-D is
  an *auxiliary coordinate variable* linked via the `coordinates`
  attribute. Terrain-following altitude specifically is a *parametric
  vertical coordinate* — 1D level dim + `formula_terms` referencing 2D
  surface fields; the 3D array is computed, never stored as a dim.
- **xarray**: N-D "non-dimension coordinates" are carried through ops
  and displayed but can never be label-indexed (`.sel` errors).
- **Iris**: `DimCoord` (1D, monotonic, one per axis) vs `AuxCoord`
  (any-D, with an explicit map of which data axes it spans); derived
  3D altitude built lazily via `HybridHeightFactory`.
- **DimensionalData.jl/Rasters.jl**: lookups are strictly 1D by type;
  bare axes use `NoLookup`; curvilinear coordinates are derived
  (`mappedcrs`), not stored N-D.

Decisions for this design:

1. **N-D `Dim` is officially a stopgap** encoding of an auxiliary
   coordinate. No new API is designed around it (no `spans` field, no
   clever N-D `length`); anything N-D-shaped waits for the Phase 5
   axes-vs-coordinates split, which matches the ecosystem consensus.
2. **Interim semantics copy xarray**: N-D dims are carried and
   displayed, never indexed or introspected.
3. **Quarantine by dispatch, not runtime checks**: define
   `const AxisDim = Dim{<:AbstractVector}` and write every 1D-only
   method (`first`, `last`, `issorted`, `isequispaced`, `step`,
   `extrema`, `window_indices`, `get_index`, non-`:` indexing) against
   `AxisDim`, so misuse fails as a `MethodError` at the call boundary
   instead of scattered `ndims(d) != 1` guards. `is_cartesian(grid)` =
   all dims are `AxisDim`. Base-semantics functions (`length`, `eltype`,
   `ndims`) forward from `values` unchanged for any dimensionality.
4. The back-compat bridge must still store/return the 3D array
   (`var.dims["z"]`), so `Dim` cannot restrict `values` to
   `AbstractVector` outright. Because the fields are Private and the
   public surface is minimal, migrating N-D `z` to a proper
   auxiliary-coordinate slot later is non-breaking.
5. Longer-term option worth stealing (CF/Iris): if ClimaCore output ever
   carries formula terms, represent `z` as 1D levels + 2D surface field
   + formula and stop materializing the 3D array entirely.

### Multiple independent columns (user requirement #1) — to think through

The other grid kind `AbstractGrid` was reserved for: data on scattered
columns (site list, unstructured sampling) rather than a lon × lat
product. Worth designing on paper now because it likely *clarifies* the
N-D question — the two problems are the same problem:

- A columns grid is naturally **one "points" axis** with *several*
  coordinates laid over it (lon(points), lat(points), maybe elevation).
  That is: multiple coordinates sharing one axis.
- Topography `z` is **one coordinate** laid over *several* axes.
- Both violate the current `Dim` assumption (exactly one coordinate ↔
  one axis, 1:1). The axes-vs-coordinates split resolves both at once:
  bare named axes + coordinate variables that each record which axes
  they span. Solving columns without that split just re-creates the
  N-D hack in transposed form.

Prior art: CF's UGRID / discrete-sampling-geometries (coordinates as
1D arrays over a station/points dim); xarray `stack`/`MultiIndex`
(tuple-valued coordinate over the stacked dim); ClimaCore column
spaces; flat.jl's `flatten` is already the ad-hoc version of this.

Candidate shapes to evaluate (not decided):

1. **Phase 5 `stack`/`unstack`** (`stack(grid, ("lon", "lat") =>
   "points")`): a `Dim` whose values are coordinate *tuples* — cheap,
   stays inside the current `Dim` model, but tuple-valued coordinates
   are awkward for units/bounds (per-component) and selectors.
2. **Several 1D `Dim`s sharing one data axis**: `Grid` drops the
   one-dim-per-axis invariant; `dim_index` becomes many-to-one. Closer
   to CF/UGRID, but breaks `size(grid)`-from-dims and structural
   indexing assumptions — effectively forces the axes/coords split.
3. **A dedicated `ColumnsGrid <: AbstractGrid`**: keeps `Grid` clean;
   the discipline rule (all code goes through accessors/verbs, never
   `grid.dims`) exists precisely so this can drop in. Risk: verb surface
   duplicated per grid kind.

Open questions: which verbs are even meaningful on columns
(`permutedims`? `cat` along points?); what `is_cartesian` returns;
whether interpolation/`resampled_as` degrade gracefully (nearest-
neighbor over points vs error); how `flatten`/`unflatten` unify with
`stack`/`unstack`.

### Breaking vs non-breaking — leaning breaking

GRID_DESIGN.md decision 6 prefers non-breaking, carried by heavy
back-compat machinery: `getproperty` materialization of
`var.dims`/`var.dim_attributes` as live-array `OrderedDict`s, the
dict-wrapping `Grid` constructor bridge, and derived
`dim2index`/`index2dim`. Current lean: **make it honestly breaking and
have users migrate**. Note the project's norm is that minor bumps do
*not* break users (GRID_DESIGN.md decision 6 misstates this), so a
breaking v0.6.0 departs from project practice and must be loudly
communicated — though Julia's pre-1.0 compat semantics (`"0.5"` won't
resolve to 0.6.x) do protect downstream packages from auto-upgrading
into it.

The case for breaking has gotten stronger as the design evolved:

- The bridge was already leaky (dict *insertion* silently unsupported);
  promoting `units`/`bounds` to fields made it leakier (*writing*
  `var.dim_attributes[name]["units"]` through the materialized dict is
  now unsupported too). Each hazard is a silent no-op or desync, not an
  error — the worst failure mode for users.
- The materialization machinery is a permanent complexity tax
  (`getproperty` forwarding, live-array aliasing semantics, per-access
  dict construction cost that Phase 2 then works to avoid) in service
  of an idiom (`var.dims["lon"]`) the new API deliberately replaces
  with a better one (`var.grid["lon"]`).
- The minimal-Public surface chosen for the grid API means the
  committed new interface is small — a migration guide covers it in a
  page (`var.dims[name]` → `var.grid[name]`, `var.dim_attributes` →
  `units`/`attributes` accessors, dict constructors → `Grid`/`Dim`).

What breaking would change in the plan: Phase 1 replaces the
`getproperty` sugar with **deprecation errors or depwarns that name the
replacement** (a hard `error("var.dims is removed, use var.grid[...]")`
is kinder than a silent desync); the ~50 `var.dims == OrderedDict(...)`
test assertions get rewritten against `Grid` equality instead of passing
via materialization; the ~39 internal rebuild sites migrate in Phase 1-2
rather than compiling through sugar; NEWS.md gains a full was/now
migration table. Open sub-questions: keep the dict-accepting
`OutputVar`/`Grid` constructors (cheap, harmless, ease migration) even
in the breaking version? One transitional minor version with read-only
depwarned materialization before removal, or cut directly at v0.6.0?

Tier A–C backlog items (`indices`, `groupby`, `sort`, `shift_longitude`,
`Grid(::NCDataset, dim_names)`, `spatial_dims`/`temporal_dims`,
`bounding_box`, `intersect`/`common_grid`, `coarsen`/`refine`,
`pairs`/`NamedTuple`, `axes`/`eachindex`, `validate`, `mesh`, `unique_dim`,
`in`, DimensionalData interop, `dims_mismatch_message`) are
optional, not scope commitments — see the Extended API backlog in
GRID_DESIGN.md.

## Other thoughts
- `subgrid` function? I am not sure what this would be but I got reminded of it
  by looking at `subset`
