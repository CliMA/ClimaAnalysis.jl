# Replace `OutputVar.dims` with a `Grid`/`Dim` abstraction

## Context

`OutputVar` stores its grid as four co-varying fields (`dims::OrderedDict{String,T}`, `dim_attributes::OrderedDict{String,C}`, plus derived caches `dim2index`/`index2dim`), duplicated again in `Metadata` and faked in `FlatVar`. Every operation that transforms a variable edits these in lockstep: ~39 rebuild sites, triplicated grid-comparison logic (with one live bug), repeated `deepcopy`/`empty`-then-refill boilerplate, and no home for per-dimension metadata like time bounds. The goal: a first-class grid type owning dimension names + coordinate arrays + per-dim attributes together, so grid manipulations become a small verb vocabulary, future grid kinds (independent columns) and per-dim metadata (bounds, reference date) have a home, and `Var.jl` shrinks.

## Decisions locked in with the user

1. **One general `Grid` type** now (not `CartesianGrid` — coordinate arrays may be N-D, e.g. 3D topography `z`); `AbstractGrid` reserved for future kinds. `is_cartesian(grid)` predicate generalizes `is_z_1D`.
2. **`Dim` struct** bundles name (runtime `String`, not type-level) + live coordinate array + attribute dict.
3. **`dim2index`/`index2dim` deleted** — derived on the fly from ordered dim names.
4. **Interior mutability kept at the content level**: `Dim` is an immutable struct holding live mutable arrays/dicts, so `transform_dates!`, `set_reference_date!`, `reverse_dim!` keep working via `.=`/`reverse!`. (Structural in-place edits — swapping a whole `Dim` — are out, see container decision below; the one affected site, `center_longitude!`, switches to `.=`.)
5. **Back-compat via `Base.getproperty`**: `var.dims` / `var.dim_attributes` materialize genuine `OrderedDict`s on demand whose values ARE the live arrays. Element mutation keeps working; dict-level insertion becomes unsupported (documented; only 2 such sites exist and both get fixed).
6. Non-breaking preferred; the field swap itself ships as v0.6.0 (project treats minor bumps as breaking).
7. Do-notation rebuild function: out of scope (withdrawn).
8. **Type stability** everywhere except attribute dicts (`Dict{String,Any}` accepted there).
   **Container decision (user's call): tuple-backed Grid, Vector door kept open.** `dims::D where D <: Tuple{Vararg{Dim}}` gives concrete per-dim types even for mixed grids (`Vector{DateTime}` time + `Float32` lon — which today's `OrderedDict{String,T}` collapses to abstract `T`), and inferable structural indexing. Accepted costs: string lookup `grid["time"]` is type-unstable (mitigated by function barriers in hot paths — already the pattern — and const-prop of literal names); `OutputVar{G}` specializes per file layout (compile-time, not runtime cost); no in-place `replace_dim!`; the empty-grid typed doctest (docs/src/var.md:188) must be updated since an empty tuple carries no eltype memory. **Discipline rule**: all `OutputVar`-facing code uses only the accessor/verb API (`dim`, `dim_names`, `dim_index`, `dropdims`, `getindex`, `replace_dim`, `permutedims`, `cat`, `ones`) — never `grid.dims` directly — so a Vector-backed `AbstractGrid` sibling can drop in later without touching `Var.jl`.
9. **Performance**: soft goal — avoid regressions; internal hot helpers use direct grid access, materialization only at user-facing back-compat sites.

## API design

New file `src/grid.jl`, `include`d in module `Var` before the `OutputVar` struct ([src/Var.jl:85](src/Var.jl#L85)). The conventional-name constants and pure string helpers move there from `outvar_dimensions.jl` (which keeps var-level delegators).

### Types

```julia
struct Dim{A <: AbstractArray, T}
    name::String                 # as in the source file, e.g. "lon"
    values::A                    # live; usually 1D, N-D legal (topography z)
    attributes::Dict{String, T}  # units, bounds, ... (Any-valued OK per decision 8)
end

abstract type AbstractGrid end

struct Grid{D <: Tuple{Vararg{Dim}}} <: AbstractGrid
    dims::D                      # order == data axis order; names unique (checked)
end
```

- The tuple IS the order: absorbs `index2dim`; `dim_index(grid, name)` = `findfirst` over names (≤~4 dims, O(n) fine) absorbs `dim2index`.
- **Type stability**: each `Dim` keeps its concrete array type even in mixed grids (`Vector{DateTime}` time + `Vector{Float32}` lon) — strictly better than today's `OrderedDict{String,T}`, which promotes mixed `T` to an abstract type. Structural indexing `grid[1, :, 1:10]` is fully inferable. String lookup `grid["time"]` is type-unstable by nature (runtime name → heterogeneous return); hot paths extract arrays once behind function barriers (already the pattern in interpolant construction and reductions), and literal names benefit from constant propagation.
- Empty grid `Grid(())`: `isempty` true, `size(grid) == ()` — required by 0-dim reductions and `_split_along_dim`. An empty tuple carries no eltype memory, so the materialized empty `var.dims` loses its element type — the docs/src/var.md:188 jldoctest output must be updated in Phase 1 (accepted small break).
- Mutation discipline: element-level mutation of `values`/`attributes` contents is supported (live arrays/dicts); structural edits (add/remove/swap a `Dim`) are always functional — the grid is fully immutable. `center_longitude!` switches from key assignment to `.=` (shape-preserving `circshift`).
- **Discipline rule for future flexibility**: `Var.jl` and friends never touch `grid.dims` directly — only accessors/verbs — so a Vector-backed `AbstractGrid` sibling can be added later if compile-time specialization churn (each file layout specializes `OutputVar{G}` methods) proves costly in practice.

### Constructors

- `Dim(name, values; attributes = Dict{String,Any}())`
- `Grid(dims::Dim...)`, `Grid()` empty, `Grid("time" => t, "lon" => lon)` pairs form
- **Back-compat bridge**: `Grid(dims_dict::AbstractDict, dim_attributes::AbstractDict = ...)` — wraps the same arrays/dicts (no copy, preserves current aliasing semantics); missing attribute entries become empty dicts.
- `check_grid_consistency(grid, data)`: the length-vs-size check from [src/Var.jl:246-261](src/Var.jl#L246-L261), called from the `OutputVar` constructor (grid alone doesn't know data); preserves the current skip-if-empty/N-D/length-1 escape hatch as one named, tightenable function.

### Accessors

| Method | Notes |
|---|---|
| `grid[name::AbstractString]` | live coordinate array; exact-then-conventional name resolution (user req. #3) |
| `dim(grid, name)::Dim`, `dim(grid, i)` | full Dim object |
| `dim_names(grid)`, `keys`, `iterate`, `length`, `haskey` (exact) / `hasdim` (conventional) | dict-flavored surface |
| `dim_index(grid, name)` | absorbs `dim2index` (~20 sites) |
| `size(grid)` (requires `is_cartesian`), `ndims`, `isempty` | data-shape source of truth |
| `units(d::Dim)`, `set_units!(d::Dim, u)`, `extrema(d)` | back `dim_units`, `set_dim_units!`, `range_dim` |
| `has_time/has_longitude/...`, `time_name/...`, `times/latitudes/...` at grid level | var versions become one-line delegators; `dates(var)` stays var-level (needs `start_date`) |
| **Dim introspection**: `length(d)`, `eltype(d)`, `ndims(d)`, `first(d)`/`last(d)`, `issorted(d)`, `isequispaced(d)` (wraps `Utils._isequispaced`), `step(d)` for equispaced | implicit computations in `_make_interpolant` (Var.jl:139-150) and `_find_extp_bound_cond` (Var.jl:216-239) become named predicates |
| **Collection semantics**: `Base.copy(grid)` (new Dims, shared arrays) vs `deepcopy`; `Base.empty(grid)`; `Base.map(f, grid)` / `Base.filter(pred, grid)` over Dims; `coordinates(grid)` = `Iterators.product` of coord arrays; `flattened_length(grid)` = `prod(size(grid))` | `empty(var.dims)` idiom (Var.jl:1674, 1705, 2157); dim filtering in `_check_dims_consistent` (Var.jl:1467-1470) and `_select`; hand-rolled `prod(length.(values(dims)))` (flat.jl:246, 253, 312) |
| **Geometry/weights**: `cell_widths(d::Dim)` (cell-edge/center spacing helpers); `latitude_weights(grid)` (cosd/area weights); `coordinate_field(grid, name)` (coordinate broadcast to full data shape) | spacing re-derivation in integration (Numerics.jl via `integrate_lon/lat`) and the interpolant's center-vs-edge longitude logic (Var.jl:178-209, 224-232); hand-rolled weights in `weighted_average_lat`/`weighted_average_lonlat` (Var.jl:866-898); reshape gymnastics in `masks._generate_binary_mask` (masks.jl:167-178) |
| **Value-based selection**: `get_index(grid, name, value, selector)` and `window_indices(grid, name, left, right)` | `get_index` var methods delegate (public `AbstractSelector` extension API preserved, outvar_selectors.jl:253-296); `_get_window_indices` (outvar_selectors.jl:306). Selector-valued structural indexing (`grid[time = NearestValue(5.0)]`) deferred to future hooks. |

### Structural indexing (user req. #4)

Base-array semantics — **integer drops the dim, range/vector keeps it** (matches `_select` at [src/outvar_selectors.jl:453](src/outvar_selectors.jl#L453), so `slice`/`select` port loss-free):

- `grid[1, :, 1:10]` positional; `grid[time = 1:10, lat = 3]` named (conventional resolution, unmentioned dims = `:`).
- Subset arrays are copies (consistent with `select`); N-D dims accept only `:` (error otherwise, same as today — future `spans` info would enable more).
- No collision with `grid[name::String]` (dispatch on element type).

### Operation verbs (each returns a new `Grid`; untouched dims share arrays, callers keep their `copy`/`deepcopy` habits via `remake`)

| Verb | Absorbs |
|---|---|
| `Base.dropdims(grid; dims)` | `_reduce_over` pop-loop (Var.jl:836-842) → every reduction |
| `grid[...]` | `window`, `_select`, `_split_along_dim` |
| `replace_dim(grid, name => new_values)` / `replace_dim(grid, old => ::Dim)` / `rename_dim` sugar (all functional) | resample refills (Var.jl:1674-1715), `transform_dates`, `_dates_to_seconds`, `shift_longitude`, `average_season_across_time`, `Atmos.to_pressure_coordinates` type-widening dance (Atmos.jl:88-102). (`center_longitude!` keeps in-place semantics via `.=` on the live array instead.) |
| `Base.permutedims(grid, perm)` (+ `permutation(grid, perm)::Vector{Int}` for data) | the ~20-line dim_attributes reorder dance (Var.jl:1549-1564), `reordered_as`, flatten's perm |
| `Base.cat(grids...; dims::String)` | dim-checking/concat parts of `Base.cat(vars...)` (Var.jl:2770-2853) |
| `Base.ones([T,] grid)` + `ones_like(var)` | 3 identical `ones_var` blocks (Var.jl:2251, 2337, Atmos.jl:195) |
| `reverse_dim!(grid, name)` (returns axis index) | `reverse_dim`/`reverse_dim!` |

### Compatibility predicates (consolidate the triplication)

`matching_dim_names(a, b; order)`, `matching_dim_values(a, b; approx)`, `matching_dim_units(a, b; on_missing)` + composites `iscompatible(a, b)` and throwing `check_dims_consistent(a, b; dims)`. Replaces `arecompatible` (Var.jl:1415, **fixing the vacuous-check bug at :1419** — y's names computed from x's), `_check_dims_consistent` (Var.jl:1446-1522), and the checks in `flatten(var, metadata)` (flat.jl:146-192; its dates-aware time case stays local).

### Interpolation

`make_interpolant(grid, data)` and `extrapolation_bc(d::Dim)` (longitude-periodic/latitude-flat detection, dispatching on `conventional_name(d)`) absorb Var.jl:131-239; callers at Var.jl:1398/1669/1717, masks.jl:308.

### Equality / show / exports

- `Dim == Dim`: name + values + attributes. `Grid == Grid`: pairwise, in order. Matching `hash` methods. `grid == OrderedDict` NOT defined — the ~50 test assertions pass because `var.dims` materializes a genuine OrderedDict.
- `Base.isapprox(a::Grid, b::Grid; kwargs...)`: `≈` on values pairing dims by conventional name — names four strictness levels together with `==`, `iscompatible`, and `matching_dim_*`; absorbs the hand-rolled `≈` checks in `Base.cat`'s off-axis comparison (Var.jl:2770-2853) and `flatten(var, metadata)` (flat.jl:168).
- `show(::Dim)` one-line, `show(::Grid)` — extracted from `show(::OutputVar)` (Var.jl:2875-2937), which reuses them.
- Export: `Grid`, `Dim`, `AbstractGrid`, `dim`, `hasdim`, `is_cartesian`, `replace_dim`, `rename_dim`, `iscompatible`, `matching_dim_*`, `coordinates`, `coordinate_field`, `latitude_weights`, `isequispaced`. Base overloads (`==`, `isapprox`, `copy`, `empty`, `map`, `filter`, indexing, `dropdims`, `permutedims`, `cat`, `ones`, introspection) need no export. No `grid(var)` function — `var.grid` is the accessor. Internal: `make_interpolant`, `check_grid_consistency`, `dim_index`, `permutation`, `cell_widths`, `window_indices`, `flattened_length`, `*_NAMES`.

### Embedding

```julia
struct OutputVar{G <: AbstractGrid, A <: AbstractArray, B}
    attributes::Dict{String, B}
    grid::G
    data::A
end
```

(No code anywhere dispatches on the old `OutputVar{T,A,B,C}` params — grep-verified, safe to change.) `Metadata` likewise gets `grid::G`; `FlatVar.getproperty` (flat.jl:325-335) forwards. `HasDimAndAttribs` survives with contract "has `attributes` and `grid` properties".

**getproperty forwarding** on `OutputVar`/`Metadata` (+ `propertynames`): `:dims` → `OrderedDict(d.name => d.values ...)` (live arrays; value type is the promotion of the per-dim array types — identical to today for homogeneous grids, and the empty grid materializes untyped), `:dim_attributes` similarly, `:dim2index`/`:index2dim` derived (kept, depwarn later). Constructors `OutputVar(attribs, dims_dict, dim_attribs, data)` and `OutputVar(dims, data)` stay (build a `Grid` internally); new primary `OutputVar(attributes, grid, data)`. `remake` gains a `grid` kwarg; `dims`/`dim_attributes` kwargs stay accepted.

## Migration plan (each phase independently mergeable and green)

### Phase 0 — Add `Grid`/`Dim` (pure addition, v0.5.x)
New `src/grid.jl` + `test/test_grid.jl` + one `include` line + runtests entry. Nothing existing changes. Docstrings marked internal/experimental so the API can still move. (~250 LOC src + ~250 LOC tests.)

### Phase 1 — Field swap + back-compat sugar (**v0.6.0**, the load-bearing phase)
Struct becomes 3 fields; add `getproperty`/`propertynames`. Sites that MUST change (field-sensitive or dict-insertion hazards):
1. 4-arg constructor (Var.jl:241-271) — builds `Grid`, keeps accepting any dict-likes.
2. `Base.copy` (Var.jl:780-785) — currently splats all 6 fields positionally; rewrite explicitly.
3. `Base.isempty` (Var.jl:771-778) — delete dead `:interpolant` filter; `isempty(attributes) && isempty(grid) && isempty(data)`.
4. `set_dim_units!` (Var.jl:710-722) — **dict insertion**, silent no-op under materialization; reimplement on grid. Exported function, highest-priority hazard.
5. `center_longitude!` (Var.jl:1267) — key replacement → `.=` (shape-preserving `circshift`).
6. `_select` (outvar_selectors.jl:450-459) — drop `typeof(var.dims)(...)` reflection.
7. `Atmos.to_pressure_coordinates` (Atmos.jl:88-102) — drop `keytype`/`valtype` widening.
8. `remake` — re-route defaults through `var.grid` (avoids double materialization).
9. Tests: rewrite the one dict-insertion line (test_Var.jl:181-184, `remake_var.dims["z"] = ...` — already corrupts dim2index today, never supported); add forwarding/live-array/empty-type tests. The ~50 `var.dims == OrderedDict(...)` assertions pass unchanged.
10. Docs: update the empty-dims jldoctest output at docs/src/var.md:188 (empty tuple grid materializes an untyped empty OrderedDict — accepted break); add `var.grid` docs + view-semantics note. NEWS.md v0.6.0 headline with was/now table.

Everything else (all ~39 rebuild sites, all `dim2index` readers, both Makie exts, all of flat.jl, `transform_dates!`/`set_reference_date!`/`reverse_dim!`) keeps compiling via the sugar. (~200 LOC changed.)

### Phase 2 — Migrate internals to grid verbs (v0.6.x, several small PRs, net −150 to −250 LOC)
In payoff order: (1) index lookups incl. both Makie exts, (2) `_reduce_over` → `dropdims`, (3) `window`/`_select` → indexing, (4) `permutedims` dance, (5) resample refills → `replace_dim`, (6) Atmos rename, (7) `ones_like`, (8) interpolant absorption, (9) `cat`/`_split_along_dim`/`shift_longitude`/masks as touched, (10) geometry absorption: `latitude_weights` into `weighted_average_lat`/`weighted_average_lonlat`, `coordinate_field` into `masks._generate_binary_mask`, `cell_widths` into the integration/extra-lon-point paths, `isapprox`/`flattened_length` into `cat` and flat.jl. Also migrate high-frequency helpers (`outvar_dimensions.jl` `has_*`/`*_name`, `dim_units`) to direct grid access — removes the per-access materialization cost.

### Phase 3 — Unify `Metadata`/`FlatVar` (v0.6.x)
`Metadata` gets `grid::G` + getproperty; `unflatten`'s hand-rolled dim2index (flat.jl:257-263) → `dim_index`. flat.jl only (~70 LOC).

### Phase 4 — Consolidate comparison logic + fix `arecompatible` bug (v0.6.x)
Grid predicates replace the triplication; regression test for the previously-vacuous name check. Watch `@test_throws` message assertions. NEWS "Bug fixes" entry. (net −80 LOC.)

### Phase 5 — Later / future hooks (design room reserved, not built)
- Depwarn `:dim2index`/`:index2dim` only (never `:dims`/`:dim_attributes` — documented public idiom, kept indefinitely as materialized views).
- `AbstractGrid` subtypes when a concrete need appears (a Vector-backed grid if tuple specialization churn hurts; specialized lon-lat grids).
- `reference_date`/`bounds` in the time `Dim`'s attributes + `bounds(d::Dim)`/cell-bounds accessors (purely additive; Var.jl:2613-2616 already anticipates the rename; `cell_widths` is the natural consumer).
- `stack(grid, ("lon","lat") => "points")` / `unstack` (xarray-style): collapse several dims into one whose coordinate is tuples of points — this is the natural representation for the multiple-independent-columns grid (user requirement #1) and the principled generalization of flat.jl's flatten.
- `insert_dim(grid, dim; at)` singleton/unsqueeze (ensemble stacking).
- Selector-valued structural indexing: `grid[time = NearestValue(5.0)]`.
- **Cell-sampling metadata** (from CF/DimensionalData): a dim-attribute convention `"sampling" => "center"|"edge"` consulted by `extrapolation_bc` and `cell_widths` when present, falling back to today's `(max−min)+Δ ≈ 360°` heuristics (Var.jl:224-232). Removes guessing at zero cost; purely additive.
- **Time bounds in CF layout**: when bounds land, store as CF-style `n×2` arrays in the time Dim's attributes — free interop with cf-xarray/Iris conventions.
- **Axes-vs-coordinates split** (from xarray): the principled resolution of the N-D dim question. N-D `z` is really a *coordinate* laid over several *axes*, not a dim; a future `AbstractGrid` refinement separating bare named axes from coordinate variables (each recording which axes it spans) is the correct fix for N-D structural indexing — preferred over bolting a `spans` field onto `Dim`. Nothing in the current design blocks this.

### Design tensions recorded (do not resolve now)
- **Cached dim traits vs interior mutability**: DimensionalData caches order/span/sampling on the lookup at construction; we recompute `issorted`/`isequispaced` per use because in-place mutation (`reverse_dim!`, `var.dims["time"] .= x`) would silently stale any cache. If the package ever moves to fully immutable grids (Phase 5+ option), adopt DD-style traits then — not before.
- **Interpolant caching keyed on grid identity**: rejected for now for the same reason — mutation through live arrays cannot invalidate a grid-keyed cache. Interpolants stay built-on-demand (`make_interpolant(grid, data)`).

### Dict-insertion hazard: resolution
Plain `OrderedDict` materialization (user's choice), insertion documented as unsupported. A `DimsView <: AbstractDict` forwarding `setindex!` is the prepared fallback if users hit it — rejected for now because it changes the printed type in doctests/REPL and insertion is semantically incoherent anyway (data can't grow an axis; today it silently desyncs `dim2index`).

## Verification

- Per phase: full suite `julia --project=test` via Pkg.test conventions (runtests includes Aqua, format check, and `test/doctest.jl` which runs Documenter doctests — the empty-dims doctest at docs/src/var.md:188 is updated deliberately in Phase 1; every other doctest must pass unchanged).
- Phase 1 acceptance: entire existing suite green with **no test deletions** — only the one documented dict-insertion rewrite; new tests for live-array mutation through `var.dims`, `propertynames`, materialized empty-dict type.
- Phase 2+: behavior-preserving refactors under the existing suite; Makie/GeoMakie ext tests cover the `index2dim` migration.
- `arecompatible` fix: new regression test for mismatched-conventional-name pairs.

## Prior-art comparison: how other libraries model the grid

### xarray (Python)
- **Model**: `DataArray` = data + `dims` (a tuple of *names only*) + `coords` (a dict of coordinate variables — each itself a mini-array with its own dims and attrs) + `attrs`. There is no single "grid" object; the dims/coords/indexes triple plays that role (a recognized design wart that the later flexible-indexes work tried to patch).
- **Direct analogs**: an xarray *coordinate variable* ≈ our `Dim` (values + attributes bundle); "dimension coordinates" (1D, same name as a dim) ≈ our 1D dims; **"non-dimension coordinates" (N-D auxiliary coords) ≈ our topography 3D `z`** — xarray legitimizes exactly the non-Cartesian case that forced our `Grid`-not-`CartesianGrid` decision.
- **Indexing**: `.isel(time=0)` positional vs `.sel(time=val, method="nearest")` label-based — our `grid[time = 1]` vs selector layer (`NearestValue`/`MatchValue`) mirrors this split; integer drops the dim in both.
- **What we're not adopting**: automatic alignment/broadcasting by dim name on binary ops (ClimaAnalysis deliberately requires `arecompatible` — explicit over implicit); pandas-index machinery.
- **cf-xarray** is the analog of our conventional-name layer (`LONGITUDE_NAMES` etc.) and exposes CF cell `bounds` — validation that per-dim bounds metadata (our future hook) is the standard home.

### DimensionalData.jl
- **Model**: `DimArray(data, (X(lookup), Ti(lookup)))` — dimension *identity in the type domain* (`X`, `Ti`, or `Dim{:name}`), dims in a **Tuple** (the container we chose), each wrapping a *lookup* carrying values, metadata, and traits: order (`ForwardOrdered`/…), span (`Regular(step)`/`Irregular(bounds)`/`Explicit(bounds_matrix)`), sampling (`Points`/`Intervals`).
- **Direct analogs**: selectors `At`/`Near`/`Between`/`Contains` ≈ our `MatchValue`/`NearestValue`/`window_indices`; `rebuild(...)` is their universal `remake`; the tuple-of-dims container validates our Phase-0 choice at ecosystem scale.
- **Where we deliberately diverge**: type-level dim names don't fit runtime NetCDF names (our `name::String` field), and their lookup-trait lattice is heavier than we need now — but `Regular/Irregular/Explicit` span is precisely the mature version of our `isequispaced`/`cell_widths`/future-`bounds` trio, and the natural blueprint if `Dim` ever grows typed structure.

### YAXArrays.jl
- **Model**: datacubes built *on* DimensionalData dims, oriented at out-of-core chunked data (Zarr/NetCDF). Its distinctive idea is `mapCube` with `InDims`/`OutDims`: the user declares which dims a kernel consumes and produces, and the framework loops over all remaining dims.
- **Relevance**: that declaration style is the scalable answer to requirement #5 ("separate operations on data from operations on the grid") — every ClimaAnalysis op's "grid signature" from the exploration (drop/subset/replace/permute) is an `InDims/OutDims` declaration in disguise. Worth revisiting if a general mapping facility is ever reintroduced (the withdrawn do-notation idea).

### Rasters.jl, Iris (Python), lightweight Julia options
- **Rasters.jl** (DD-based): adds geospatial semantics as lookup metadata (`crs`, `mappedcrs`, `missingval`) — the pattern for attaching domain semantics (our conventional lon/lat handling, mask metadata) to dims rather than to the array.
- **Iris** is the closest philosophical match: cubes carry explicit coordinate objects (`DimCoord`/`AuxCoord`) with **bounds arrays** and *cell measures* (precomputed area weights) — direct prior art for our `latitude_weights`/`cell_widths` and time-bounds hooks.
- **AxisKeys.jl / NamedDims.jl**: name-only or keys-only wrappers; they show names and lookup values are separable concerns — our `Dim` fuses them deliberately because ClimaAnalysis always has both.

### Takeaways adopted into this design
1. Coordinate-as-object (`Dim`) with attributes attached — xarray/Iris consensus.
2. Tuple container for dims — DimensionalData-validated, fits the type-stability requirement.
3. Positional-vs-label indexing as two explicit layers — xarray's `isel`/`sel` split, already latent in ClimaAnalysis selectors.
4. `remake`/`rebuild` as the single reconstruction seam — DD's `rebuild` pattern.
5. Per-dim `bounds`/span metadata as the future home for cell geometry and time bounds — CF/Iris/DD agreement.
6. Rejected with precedent: automatic dim-name alignment (xarray) and type-level names (DD), both misfits for explicit, NetCDF-runtime-named climate output.

## Extended API backlog (extensive sweep — optional, not scope commitments)

Ideas vetted against actual ClimaAnalysis workflows and prior art (xarray, DimensionalData.jl, Rasters.jl). Tiered by pull-in priority; anything here can be added later without design changes because it builds on the core types/verbs.

### Tier A — earns its place during Phase 2 (absorbs existing code or unlocks a current pain point)

- `indices(pred, d::Dim)` — predicate-based index selection **with do-notation support** (`indices(d) do v; a <= v <= b; end`); `window_indices` becomes the two-sided special case. Generalizes `_get_window_indices` and gives `window`/`select` an escape hatch for arbitrary conditions.
- `groupby(d::Dim, f)` → ordered `label => indices` groups. `_split_along_dim` (Var.jl:2141-2165), `split_by_season(_across_time)`, `split_by_month`, and `average_season_across_time` are all "group time indices by a function of the value, then subset/reduce" — this is their common core.
- `Base.sort(grid; dims)` → `(sorted_grid, permutations)` and `Base.issorted(grid)` — interpolation silently requires sorted dims today (`_make_interpolant` warns and bails, Var.jl:145-150); a first-class sort pairs with `reverse_dim` and removes a user footgun.
- `shift_longitude(grid, lon_0)` grid-level (returns grid + data permutation/circshift recipe) — absorbs the grid half of `shift_longitude`/`center_longitude!` (Var.jl:1251-1340), including the duplicate-endpoint `pop!` subtlety.
- `Grid(nc::NCDataset, dim_names)` constructor — absorbs the dim-building halves of `read_var(::String)` (Var.jl:361-381) and `read_var(::Vector{String})` (Var.jl:411-442).
- `Grid(template::TemplateVar)` / `Template.initialize_grid` — Template.jl:94-130 already builds exactly a grid lazily; useful standalone for tests.
- `spatial_dims(grid)` / `temporal_dims(grid)` — grouped conventional queries; several functions test `has_longitude && has_latitude` in sequence.
- `bounding_box(grid)` — lon/lat extents as a NamedTuple; masks and plotting extensions both re-derive extents.

### Tier B — genuinely useful, no current absorber (add when a caller appears)

- `intersect(a::Grid, b::Grid)` / `common_grid(a, b; by = :coarsest)` — compute the overlap/target grid for obs-vs-sim comparison before `resampled_as`; today users construct target grids by hand (docs/src/howdoi.md resampling recipes). The single strongest candidate for a new user-facing workflow win.
- `coarsen(grid; factor or spacing)` / `refine` — target-grid builders that pair with `resampled_as(var; kwargs)` (Var.jl:1741-1774, which currently fabricates a throwaway dest var).
- `Base.pairs(grid)` (`name => values` iteration) and `NamedTuple(grid)` (when names are valid identifiers) — ergonomic destructuring: `(; lon, lat) = NamedTuple(grid)`.
- `Base.axes(grid)`, `CartesianIndices(grid)`, `eachindex(grid)` — completes the array-protocol story next to `size`/`coordinates`.
- `validate(grid; strict)` — one diagnostic pass: duplicate names, unsorted dims, NaN/duplicate coordinates (the `shift_longitude` endpoint-dedup case), missing units. Today these are scattered warnings.
- `mesh(grid)` — tuple of `coordinate_field`s (meshgrid); companion to `coordinate_field` for vectorized geometry computations.
- `unique_dim(d)` / dedup helper — duplicate-coordinate repair (wrap-around longitude point, repeated time stamps from restarted simulations — `Base.cat` would benefit when concatenating overlapping runs).
- `Base.in((; lon, lat), grid)` point-membership / bounds check — guards extrapolation errors before calling `var(target_coord)`.

### Tier C — interop and display polish

- DimensionalData.jl interop extension (`DimArray(var)` / `OutputVar(::AbstractDimArray)`): a package extension mapping `Dim` ↔ `DimensionalData.Dim{name}` — opens the whole DD/Rasters ecosystem (plotting, zonal stats) without coupling the core.
- `Base.summary(grid)` compact one-liner (`"Grid(time × lon × lat, 120×360×180)"`) for use in `show(::OutputVar)` headers and error messages.
- Error-message helpers: `dims_mismatch_message(a, b)` producing the aligned name/units/extents diff — the throwing predicates (`check_dims_consistent`) and `cat` all want the same readable diff.

### Deliberately rejected

- `push!`/`delete!` on a grid embedded in an OutputVar — desyncs data by construction; structural edits stay functional.
- `grid == OrderedDict` cross-type equality — would paper over the type change and create asymmetric `==`.
- Type-level dim names (DimensionalData-style `X`/`Ti` types) — names are runtime NetCDF strings here; revisit only if profiling shows string-lookup dispatch cost that constant propagation doesn't already eliminate.
