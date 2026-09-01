# Auxiliary coordinates across ecosystems: a deep survey

How mature gridded-data software models, discovers, propagates, and selects on
coordinates that are not simple 1D per-axis vectors — CF auxiliary coordinates
(2D curvilinear lat/lon, 3D topography-warped altitude), parametric vertical
coordinates, bare index dimensions, and discrete-sampling-geometry point axes.
Written to inform the ClimaAnalysis.jl Grid/Dim design (see GRID_DESIGN.md,
API.md; prior shallow survey in aux_coord_examples/HANDOFF.md — this document
goes deeper on exact semantics).

Behaviors marked **[verified]** were reproduced locally (xarray 2026.7.0,
iris 3.16.0, scripts in the session scratchpad) against the two example files
in `aux_coord_examples/` (`rasm.nc` — curvilinear `Tair(time,y,x)` with
`coordinates = "yc xc"`; `ta_with_topography.nc` — ClimaDiagnostics output with
`z_physical(z_reference,lat,lon)` and **no** `coordinates` attribute).

---

## 1. Executive summary

Every mature system converges on the same data model: **a coordinate is a
named array + attributes, attached to a variable by a mapping onto a subset of
the variable's dimensions; a "dimension coordinate" is the special case of a
1D coordinate over its own dimension.** That is xarray's dimension/
non-dimension coordinate split, iris's `DimCoord`/`AuxCoord` + `coord_dims`,
netCDF-Java's `CoordinateAxis` within a `CoordinateSystem`, the CF data
model's `DimensionCoordinate`/`AuxiliaryCoordinate` constructs, and CF-netCDF
chapter 5 itself. Nobody models "the grid" as a first-class object *inside*
the data variable (uxarray does bolt a sidecar `Grid` object onto xarray — the
exception that proves the rule, and exactly the design the CF data model paper
rejected).

The systems differ sharply in three places — **propagation under reduction,
conflict handling, and selection** — and those differences are the real design
decisions:

| System | Model | Discovery | Slice/index aux | Reduce over spanned dim | Select by aux values |
|---|---|---|---|---|---|
| xarray | coords = DataArrays, dims ⊆ var dims | `coordinates` attr only (`decode_coords`); no heuristics | sliced; int index → scalar coord kept | **coord dropped** silently | opt-in index (`set_xindex`, `NDPointIndex` nearest-only, 2025+) |
| iris | `AuxCoord` + explicit `coord_dims` | CF attr + std_name/units identification | sliced; scalar coords kept (points+bounds) | **coord collapsed too** (mean of points, bounds=[min,max]) | refused: `CoordinateMultiDimError`; regrid/UnstructuredNearest instead |
| netCDF-Java | `CoordinateAxis` (may be N-D) in `CoordinateSystem` | layered: coord vars → CF attrs → `_Coordinate*` attrs → units/name heuristics | read-only layer; N/A | N/A (not an analysis lib) | `findXYindexFromLatLon` index search on 2D axes |
| cfdm/cf-python | full CF data model incl. `DomainAncillary` | CF attrs (strict) | subspace propagates all constructs | 1D coords collapsed w/ bounds; **≥2D aux coords removed** | yes: `indices` conditions on N-D coords (mask-based, compress/envelope/full) |
| CDO | curvilinear/unstructured = first-class grid types | grid module (per-file griddes) | operators grid-aware | result carries the reduced grid | `remapnn`/remap operators, not label indexing |
| Julia (DD/Rasters/YAX) | strictly 1D lookup per dim | dim-name match only | 1D lookups only | lookup reduced to length-1 | none (`MergedLookup` scan; `Near` throws) |
| ClimaAnalysis (adopted) | `coord = (name, spans ⊆ dims, values, attrs)` | chain: CF attr → heuristics → name conventions | slice along spans | drop-with-provenance recommended (§11) | explicitly not supported (error + pointer to resample) |

Headline findings, expanded in the body:

1. **xarray drops, iris reduces.** `mean` over a spanned dim silently deletes
   a non-dimension coordinate in xarray [verified]; iris collapses the
   coordinate alongside the data into a scalar cell with `bounds=[min,max]`
   [verified]. This is the single biggest semantic fork in the survey.
2. **Nobody does label/orthogonal indexing on N-D coordinates.** The universal
   punt: masks (`where`), nearest-neighbor trees (xoak, `NDPointIndex`,
   pyresample, CDO `remapnn`), or "regrid first". Iris raises a dedicated
   exception type for it.
3. **Discovery is layered everywhere except xarray.** xarray reads only the
   `coordinates` attribute; iris, netCDF-Java, cf_xarray, and
   CommonDataModel.jl all add standard_name/units/axis heuristics because
   real files (MPAS, ClimaDiagnostics; WRF's attribute is present but
   incomplete) omit or under-fill the attribute.
4. **Scalar coordinates are load-bearing.** Both xarray and iris keep a 0-d
   coordinate when indexing drops all spanned dims [verified]; iris's
   merge rebuilds dimensions from them. They are the memory of where a slice
   came from.
5. **No Julia package has any of this.** DimensionalData/Rasters/YAXArrays are
   strictly one-1D-lookup-per-dimension; the only multidim machinery
   (`MergedLookup`, `Transformed`) solves different problems. GitHub issue
   traffic confirms the gap is known and unaddressed.

---

## 2. xarray — exact semantics

### 2.1 The model

xarray distinguishes **dimension coordinates** ("one dimensional coordinates
with a name equal to their sole dimension", index-backed) from
**non-dimension coordinates**, which "can be multidimensional … and there is
no relationship between the name of a non-dimension coordinate and the
name(s) of its dimension(s)". The docs are explicit about their second-class
status: non-dimension coordinates "can be useful for indexing or plotting;
otherwise, xarray does not make any direct use of the values associated with
them. They are not used for alignment or automatic indexing, nor are they
required to match when doing arithmetic"
(https://docs.xarray.dev/en/stable/user-guide/data-structures.html).

### 2.2 Propagation rules [verified on rasm.nc, xarray 2026.7.0]

For aux coords `xc(y,x)`, `yc(y,x)` on `Tair(time,y,x)`:

| Operation | Result for `xc(y,x)` |
|---|---|
| `isel(x=slice(0,10))` | sliced along x → `xc(y, 10)` |
| `isel(x=5)` (dim dropped) | loses the dim → `xc(y)` |
| `isel(x=5, y=7)` | **0-d scalar coordinate retained** (value kept, `ndim==0`) |
| `isel(x=5, drop=True)` | `xc(y)` still kept — `drop` only discards coords that would *become scalar* ("drop coordinates variables indexed by integers instead of making them scalar", isel API docs) |
| `mean("x")` | **`xc`, `yc` dropped entirely** — even though they still span `y` |
| `mean("time")` (unspanned) | kept unchanged |
| `weighted(w).mean("lat")` | same drop rule as unweighted [verified on topo file] |
| `transpose("x","y","time")` | coord **array physically transposed** → `xc(x,y)` (`transpose_coords=True` default) |
| `concat` along `time` (unspanned), coords equal | kept once (compared under default `coords="different"`, `compat="equals"`) |
| `concat` along `time`, coords **conflict** | `xc` **gains the concat dim** → `xc(time,y,x)` — the notorious `coords="different"` behavior |
| `concat(..., coords="minimal", compat="override")` | first object's `xc(y,x)` wins, no comparison |
| `concat` along `x` (spanned) | concatenated along x (spanning coords are always concatenated, all `coords` modes) |
| `groupby("time.season").mean()` | survives (reduction dim `time` unspanned) → present on `(season,y,x)` result |
| `groupby_bins("xc", bins=4).mean()` | grouping **by** the 2D coord works: y,x consumed, result `(time, xc_bins)` |
| `d1 + d2` with conflicting `xc` | **`xc` silently dropped** from the result; matching `yc` kept |
| `sel(xc=..., method="nearest")` | `ValueError: Could not automatically create PandasIndex for coord 'xc' with 2 dimensions. Please explicitly set the index using set_xindex` |
| `swap_dims({"x": "xc"})` | `ValueError: replacement dimension 'xc' is not a 1D variable along the old dimension 'x'` |
| `squeeze()` | squeezed dims' coords → scalar coords (`drop=True` to discard), squeeze API docs |

Notes on the two sharp edges:

- **The reduction drop is barely documented.** The user guide's Coordinates
  section never states it; it is confirmed by issue traffic:
  "Applying reduction operations like .mean() to a Dataset or DataArray will
  remove all coordinates that have the reduced dimension"
  (https://github.com/pydata/xarray/issues/8317, which asks for an opt-in
  `reduce_coords`; also #9168, #3510). Documented workarounds: `reset_coords()`
  to demote coords to data vars before reducing, or `keepdims=True` (#2170).
- **concat's `coords` argument** (https://docs.xarray.dev/en/stable/generated/xarray.concat.html):
  `"minimal"` = "Only coordinates in which the dimension already appears are
  included"; `"different"` (default) = "Coordinates which are not equal
  (ignoring attributes) across all datasets are also concatenated"; `"all"` =
  everything. Non-concatenated variables are checked per `compat`
  (identical/equals/broadcast_equals/no_conflicts/override). The default is
  both surprising (a static coord can sprout a time dim) and slow (it eagerly
  loads and compares coord values across files); the io guide recommends
  `data_vars="minimal", coords="minimal", compat="override"` for
  many-file opens.
- Binary-op coordinate merging uses `compat="minimal"` internally, whose
  implementation `coord_names.intersection_update(variables)` is what silently
  drops conflicting non-index coords
  (xarray/structure/merge.py, `merge_coordinates_without_align`).

### 2.3 Discovery: `decode_coords`

Semantics of `xr.open_dataset(..., decode_coords=...)`
(https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html):

- `True`/`"coordinates"` (default): promote variables named in each
  variable's (or the dataset's) `coordinates` attribute. "Only existing
  variables can be set as coordinates. Missing variables will be silently
  ignored."
- `"all"`: additionally promote variables referenced by the attributes in
  `CF_RELATED_DATA = ("bounds", "grid_mapping", "climatology", "geometry",
  "node_coordinates", "node_count", "part_node_count", "interior_ring",
  "cell_measures", "formula_terms")` (xarray/conventions.py; the last three of
  `grid_mapping`/`cell_measures`/`formula_terms` need `key: value` parsing).
  **`ancillary_variables` is not on the list.**
- `False`: nothing promoted; attribute stays in `attrs`.

There are **no heuristics**: if `coordinates` is absent, the only automatic
promotion is structural — a variable whose name equals a dimension name
becomes a (dimension) coordinate (`merge_core` in xarray/structure/merge.py).
[verified: opening `ta_with_topography.nc` with `decode_coords="all"` promotes
`time_bnds`/`date_bnds` (via `bounds`) but leaves `z_physical` and `date` as
data variables — nothing references them.]

Round-trip: on decode the `coordinates` attribute is *moved* from `attrs` to
`encoding` [verified: `Tair.encoding["coordinates"] == "yc xc"`]. On save,
xarray auto-writes `coordinates` as "a space-delimited list of names of
coordinate variables that share dimensions with the DataArray being written";
orphan coords go to a **global** `coordinates` attribute — "This is not
CF-compliant but again facilitates roundtripping"
(https://docs.xarray.dev/en/stable/user-guide/io.html). `coordinates` present
in both attrs and encoding raises `ValueError`.

### 2.4 Selection on N-D coords: the explicit-indexes endgame

- Historical stance: `.sel` requires an index; only 1D dimension coordinates
  got one automatically. The official recipe for 2D coords was `where(...,
  drop=True)` masks and `groupby_bins`
  (https://docs.xarray.dev/en/stable/examples/multidimensional-coords.html).
  [verified: `where((xc>200)&(xc<220), drop=True)` trims to the bounding box
  of the mask, `(time, y:103, x:107)`.]
- The **flexible indexes** project
  (https://github.com/pydata/xarray/blob/main/design_notes/flexible_indexes_notes.md)
  turned the dimension↔index relationship many-to-many: any set of coords can
  back a custom `Index` (KDTree, BallTree, S2 …), set via
  `set_xindex` (added v2022.09.0; the explicit-indexes refactor itself shipped
  v2022.06.0). Design notes explicitly scope *out*: chunked/lazy coordinate
  indexing (delegated to index implementations), mutable indexes, and
  automatic resolution when multiple indexes could satisfy a `.sel` (error
  instead). Prior art credited: **xoak**, "point-wise selection of irregular,
  n-dimensional data encoded in coordinates with an arbitrary number of
  dimensions" via scipy cKDTree/sklearn BallTree/S2 adapters
  (https://xoak.readthedocs.io/en/latest/).
- **`xarray.indexes.NDPointIndex`** (v2025.07.1, PR #10478): KDTree over a set
  of same-shaped N-D coords. [verified:
  `da.set_xindex(("yc","xc"), NDPointIndex)` then
  `sel(yc=57.2, xc=118.3, method="nearest")` returns the nearest cell;
  exact `sel` raises "NDPointIndex only supports selection with
  method='nearest'"; and the index **does not survive `isel` slicing** — after
  `isel(x=slice(0,50))` only the time PandasIndex remains.]
- Still unsupported, deliberately: label/orthogonal **slice** indexing on N-D
  coords (`sel(lat=slice(40,60))` on 2D lat) — open feature request
  (https://github.com/pydata/xarray/issues/10572); no rationale documented
  beyond per-index capability scoping.

### 2.5 cf_xarray

Criteria table (from cf_xarray/criteria.py and
https://cf-xarray.readthedocs.io/en/latest/coord_axes.html): identification by
`standard_name` (latitude/longitude/vertical set), `units`
(degrees_north family etc.), `axis` (X/Y/Z/T), `positive` (up/down ⇒
vertical), `_CoordinateAxisType` (netCDF-Java's convention!), `cartesian_axis`
(GFDL), `grads_dim`, plus **name regexes** as a last resort (e.g. latitude
`y?(nav_lat|lat|gphi)[a-z0-9]*`, Z
`(z|nav_lev|gdep|lv_|[o]*lev|bottom_top|sigma|h(ei)?ght|altitude|depth|isobaric|pres|isotherm)[a-z_]*[0-9]*`).
`decode_vertical_coords()` evaluates CF `formula_terms` (chosen by the
parametric `standard_name`) into a new physical-height variable —
the direct analog of materializing `z_physical`
(https://cf-xarray.readthedocs.io/en/latest/parametricz.html).
[verified on `ta_with_topography.nc`: `ds.cf.axes == {X: [lon], Y: [lat],
Z: [z_reference], T: [time]}` from `axis` attrs alone, but `ds.cf["vertical"]`
raises KeyError and `z_physical` is invisible — heuristics identify axes, they
do not attach unreferenced aux coords.]

### 2.6 Warts to learn from

- **Index/coordinate conflation** took a multi-year refactor to undo
  (pre-2022: index ⟺ dimension coordinate; MultiIndex levels were "virtual
  coordinates"; MultiIndex couldn't serialize, #1077). Post-refactor edges
  remain (#7695 alignment with multiple indexes per dim, #8056 stale indexes
  after `assign_coords`).
- **Silent drops** (reductions, arithmetic conflicts) are the top recurring
  complaint: #8317, #9168, #3510.
- **concat defaults** are a performance and correctness trap (#2975 concat of
  conflicting dim coords "pretty weird/counterintuitive"; discussion #9058 on
  alignment costs; io guide's recommended non-default flags).

---

## 3. Iris — exact semantics

*(Complements §2; iris is the most CF-faithful in-memory model. Doc-level
claims below cross-checked empirically where marked.)*

### 3.1 The model

A `Cube` = data array + list of coordinate objects, each attached with an
explicit `coord_dims` tuple mapping coord axes to cube axes. `DimCoord` is 1D,
**numeric and strictly monotonic**, at most one per cube dim; `AuxCoord` is
any-rank, attached to any tuple of dims (including `()` — scalar coords).
`CellMeasure` (area/volume weights) and `AncillaryVariable` attach the same
way. Cubes hold at most one standard_name per coord, and coordinate identity
is by rich metadata comparison, not object id.
(https://scitools-iris.readthedocs.io/en/stable/userguide/iris_cubes.html)

### 3.2 Propagation rules [verified on rasm.nc, iris 3.16.0]

Loading rasm.nc: `yc`/`xc` become `AuxCoord`s with inferred
`standard_name="latitude"/"longitude"` (file has only units + long_name —
iris identifies by units), `coord_dims=(1,2)`.

| Operation | Result for `AuxCoord latitude(y,x)` |
|---|---|
| `cube[:, :, 0:10]` | sliced → shape (205,10), coord_dims (1,2) |
| `cube[:, :, 5]` | loses dim → shape (205,), coord_dims (1,) |
| `cube[:, 7, 5]` | **scalar coord kept**: shape (1,), `coord_dims=()`, value retained |
| `collapsed("time", MEAN)` | untouched; `time` itself becomes a **scalar coord** (per `Coord.collapsed`, points→aggregate and bounds→extremes; bounds verified below for the lat/lon collapse) |
| `collapsed(["latitude","longitude"], MEAN)` | **coords collapsed too**: scalar cells with `points=mean(values)`, `bounds=[[min, max]]` (lat: point 53.16, bounds [16.53, 89.79]) |
| `collapsed("latitude", MEAN)` (one 2D coord) | collapses **every dim the coord spans** — result shape (36,), i.e. y *and* x gone |
| `transpose([2,1,0])` | coord array **not** touched; `coord_dims` remapped (1,2)→(1,0) |
| `concatenate_cube` along time | requires aux coords equal; conflicting values → hard `ConcatenateError: "Auxiliary coordinates are unequal for phenomenon …"` |
| `merge_cube` of two scalar-time cubes | **rebuilds the time dimension** from scalar coords → (2,205,275) |
| `Constraint(latitude=lambda cell: 60<cell<65)` | `iris.exceptions.CoordinateMultiDimError: Cannot apply constraints to multidimensional coordinates` (iris/_constraints.py:338) |
| `interpolate([("latitude",[60])], Linear())` | `ValueError: Coordinates repeat a data dimension - the interpolation would be over-specified.` |

The collapse behavior is documented at `Cube.collapsed`: coordinates spanning
the collapsed dimension are themselves collapsed (points become the cell
midpoint/aggregate, bounds the extremes), preserving provenance — the exact
opposite philosophy to xarray's drop. Iris likewise *keeps* scalar coords on
slicing so that `merge` can later reconstruct dimensions from them
(the merge/concatenate pair:
https://scitools-iris.readthedocs.io/en/stable/userguide/merge_and_concat.html).

Doc-level confirmation of the table (all
https://scitools-iris.readthedocs.io):

- Indexing: "Cube indexing … implemented at the data level … All metadata
  will be subsequently indexed appropriately" — every dim/aux coord, cell
  measure, and ancillary variable is sliced by the keys for its dims; the
  cube source moves a coord to the aux container "If the associated dimension
  has been sliced so the coord is a scalar" (lib/iris/cube.py).
- `Coord.collapsed` "returns a copy … with points & bounds replaced by a
  simple bounded region" (iris.coords API); a `CellMethod` (e.g.
  `model_level_number: mean`) is added; the cube-statistics guide shows
  collapsed hybrid terms like `level_height 696.67 m, bound=(0.0, 1393.33)`.
- `transpose` renumbers the dim-mapping tuples of all attached constructs
  (`remap_cube_metadata` in cube.py) — coordinate arrays untouched.
- `CubeList.merge` "combines multiple input cubes into a single resultant
  cube with new dimensions created from the *scalar coordinate values* of the
  input cubes" (multiple new dims if the scalar sequences form an orthogonal
  basis); `concatenate` extends existing dims "by joining together sequential
  dimension coordinates", requiring shape/metadata/coords consistency —
  helpers `equalise_attributes()` / `unify_time_units()` exist precisely
  because mismatches are the classic blockers (merge & concatenate guide).
- `CellMeasure` (area/volume weights) and `AncillaryVariable` attach via
  their own dims-mappings (`cell_measure_dims`, `ancillary_variable_dims`),
  parallel to `coord_dims` — the CF construct taxonomy realized in one class
  hierarchy (iris.coords API).

### 3.3 Derived coordinates (`AuxCoordFactory`)

An `AuxCoordFactory` "can manufacture an additional auxiliary coordinate on
demand, by combining the values of other coordinates" — the formula terms are
its dependencies. Concrete factories mirror CF Appendix D:
`AtmosphereSigmaFactory`, `HybridHeightFactory`, `HybridPressureFactory`
(`p = ap + b*ps`), `OceanSigmaFactory`, `OceanSigmaZFactory`, `OceanSFactory`,
`OceanSg1Factory`, `OceanSg2Factory`
(https://scitools-iris.readthedocs.io/en/latest/generated/api/iris.aux_factory.html).
Two design points worth copying:

- **Laziness**: the derived coordinate has lazy points/bounds whenever any
  dependency does (real & lazy data guide) — a 4D ROMS depth is never
  materialized unless asked for.
- **Slicing survival**: `Cube.__getitem__` records "a mapping from old
  coordinate IDs to new coordinates, for subsequent use in creating updated
  aux_factories" — the factory is re-pointed at the *sliced* dependencies and
  recomputes on demand (cube.py). Derived coords therefore behave exactly
  like stored aux coords under subsetting, at zero storage cost.

### 3.4 Selection stance, and the cost of the model

`iris.Constraint(coord_values=...)` works on any *named* coordinate, dim or
aux — "Constraint filtering is performed at the cell level" (point + bounds
comparison) — but multidimensional coords are a hard error [verified §3.2].
The sanctioned value-space paths for 2D coords are regridding:
`iris.analysis.UnstructuredNearest` ("a nearest-neighbour regridding scheme
for regridding data whose horizontal (X- and Y-axis) coordinates are mapped
to the *same* dimensions, rather than being orthogonal on independent
dimensions"), `PointInCell`, and `iris.analysis.trajectory.interpolate`
(`"nearest"` required for multi-dim coords). `cube.interpolate` with
`Linear`/`Nearest` "requires 1D numeric, monotonic, coordinates";
`AreaWeighted` requires "monotonic, bounded, 1D spatial coordinates"
(interpolation & regridding guide). So: no label-based extraction over 2D
lat/lon, ever — only explicit regrid/nearest machinery.

Pain points on record: the strict CF model makes multi-source combination
hard — users must "level the playing field" of metadata before merging
(https://github.com/SciTools/iris/issues/4446 "Lenient Cube Merge"; #3234
"Unify merge and concatenate"; #2592 `var_name` mismatch blocking merge;
#3084 merge blocked by non-scalar AuxCoords with mismatched points). Iris's
own comparison page and the xarray FAQ frame the trade: "Iris strictly
interprets CF conventions" vs xarray's "metadata should not be allowed to get
in the way". The model's overhead also showed up unprompted in this survey: a
plain `cube[...]` slice crashed on iris 3.16 + xxhash 3.6 (TypeError inside
coordinate-metadata hashing during `add_aux_coord`) — every getitem runs
metadata-equality machinery over all coords.

---

## 4. netCDF-Java CDM

The reference architecture for "decode step produces an interpreted grid".
The CDM "architecture has three layers": data access ("handles data reading
and writing"), a **coordinate system layer** ("identifies the coordinates of
the data arrays"), and a scientific feature-types layer (grids, radial,
point data)
(https://docs.unidata.ucar.edu/netcdf-java/current/userguide/common_data_model_overview.html).

- **Model**: `CoordinateAxis` is a Variable subtype, optionally typed by
  `AxisType` (Lat, Lon, Height, Pressure, Time, GeoX/GeoY/GeoZ, RunTime,
  Ensemble, …), **may be multidimensional**, with the same subset rule as CF:
  "All dimensions used by a Coordinate Axis must be shared with the data
  variable." A `CoordinateSystem` is "one or more CoordinateAxis, and zero or
  more CoordinateTransforms"; a data variable can have 0..n coordinate
  systems; a `CoordinateTransform` "currently is either a Projection or a
  Vertical Transform".
- **Discovery** is the deepest fallback chain surveyed. `NetcdfDataset`
  selects a convention-specific `CoordSysBuilder` from the global
  `Conventions` attribute; axes are then identified from (1) classic
  coordinate variables (name == dimension name), (2) the CF `coordinates`
  attribute, (3) the internal `_Coordinate*` attribute convention
  (`_CoordinateAxisType`, `_CoordinateAxes` — a space-separated axis list on
  the data variable, `_CoordinateSystems`, `_CoordinateZisPositive`,
  `_CoordinateSystemFor`, `_CoordinateAliases`), and (4) when no convention
  matches, units/name heuristics: `degrees_east`→Lon, `degrees_north`→Lat,
  udunits date units→Time, names like "level"/"sigma_level"→vertical.
  Implicit coordinate systems are assembled from coordinate variables
  matching the data variable's dimensions
  (https://docs.unidata.ucar.edu/netcdf-java/5.8/userguide/coord_attr_conv.html,
  https://docs.unidata.ucar.edu/netcdf-java/5.5/userguide/coord_system_builder.html).
- **Vertical transforms**: registered per CF Appendix D
  (`atmosphere_hybrid_sigma_pressure_coordinate`, `ocean_s_coordinate_g1/g2`,
  `atmosphere_sigma_coordinate`, …, plus non-CF `explicit_field`), recognized
  by `standard_name` + `formula_terms`. Evaluation is deliberately deferred:
  "For performance, the actual work is not done until you call …
  makeVerticalTransform(), and then VerticalTransform.getCoordinateArray()"
  — computed per time step
  (https://docs.unidata.ucar.edu/netcdf-java/dev/userguide/std_vertical_coord_transforms.html).
- **2D coords / selection**: `GridCoordSystem` allows "x and y are 1 or 2
  dimensional" (2D = curvilinear). There is no label indexing; the offered
  primitive is coordinate→index *search*: `findXYindexFromCoord` /
  `findXYindexFromLatLon` ("Given a lat,lon point, find the x,y index of the
  containing grid point"), returning −1/empty when outside; the modern grid
  rewrite keeps this shape with a dedicated `HorizCoordSys2D` for curvilinear
  search (https://docs.unidata.ucar.edu/netcdf-java/4.6/javadocAll/ucar/nc2/dt/GridCoordSystem.html,
  https://github.com/Unidata/netcdf-java/discussions/712). Subsetting is by
  index ranges / bounding box through the feature type — never by generalized
  aux-coord constraints.

---

## 5. cfdm / cf-python and the CF data model

### 5.1 The construct taxonomy (Hassell et al. 2017)

Source: "A data model of the Climate and Forecast metadata conventions",
GMD 10:4619, https://gmd.copernicus.org/articles/10/4619/2017/. Key
discriminations, verbatim:

- "a dimension coordinate construct provides monotonic numeric coordinates
  for a single domain axis, and an auxiliary coordinate construct provides
  any type of coordinate information for one or more of the domain axes";
  "a domain axis can be associated with at most one dimension coordinate
  construct, whose data array values must all be non-missing and strictly
  monotonically increasing or decreasing".
- When aux is *required*: "Auxiliary coordinate constructs have to be used,
  instead of dimension coordinate constructs, when a single domain axis
  requires more than one set of coordinate values, when coordinate values are
  not numeric, strictly monotonic, or contain missing values, or when they
  vary along more than one domain axis construct simultaneously."
- Bare axes are legal: "If a domain axis construct does not correspond to a
  continuous physical quantity, then it is not necessary for it to be
  associated with a dimension coordinate construct" — the CF blessing for
  rasm.nc's index-only `x`/`y`.
- **Formula terms are NOT coordinates**: "CF-netCDF variables named by the
  formula_terms attribute of a CF-netCDF coordinate variable correspond to
  domain ancillary constructs"; a domain ancillary "provides information
  which is needed for computing the location of cells in an alternative
  coordinate system". The `CoordinateReference` construct (grid_mapping or
  formula_terms) holds datum + conversion; "the coordinate values are not
  relevant to the coordinate reference construct, only their properties."
- **No grid object**: coordinate systems are emergent — "The domain may
  contain various coordinate systems, each of which is constructed from a
  subset of the dimension and auxiliary coordinate constructs"; the Domain
  itself "is not a construct of the data model … but an abstract concept"
  (cfdm later reifies it so a domain "may exist independently of a field
  construct", https://ncas-cms.github.io/cfdm/cf_data_model.html). Design
  principles: "a minimal set of elements … sufficient for accommodating all
  aspects of the CF conventions", "independent of the encoding".

So the full attachment taxonomy is: dimension coordinate (1D monotonic
numeric, ≤1 per axis) / auxiliary coordinate (anything else, physical
locations) / domain ancillary (inputs to coordinate *computations*) / cell
measure (area/volume weights) / field ancillary (per-cell metadata about the
*data*, e.g. quality flags — attached like coords but semantically about the
field, not the domain).

### 5.2 cf-python operational semantics

- **Subspace**: `f[...]`/`f.subspace` "also subspaces any metadata constructs
  of the field construct … which span any of the domain axis constructs that
  are affected" (https://ncas-cms.github.io/cf-python/tutorial.html).
- **Selection on N-D coords — the one system that does it**: `f.indices`
  accepts conditions on multidimensional metadata constructs
  ("Conditions may also be applied to multi-dimensional metadata
  constructs"), with three modes: `'compress'` (default; keep cells matching,
  compressing to the bounding subarray), `'envelope'`, `'full'` — mask-based,
  not label-lookup
  (https://ncas-cms.github.io/cf-python/method/cf.Field.indices.html). This is
  formalized `where`-style selection, not an index structure.
- **Collapse**: the collapsed axis's dimension coordinate becomes one cell
  with bounds = original extremes (iris-style), but the source **removes all
  ≥2-D auxiliary coordinates spanning a collapse axis** (field.py comment:
  "REMOVE all 2+ dimensional auxiliary coordinates which span this axis"),
  removes varying 1-D aux coords (keeps constant ones), and deletes
  coordinate references whose terms span collapsed axes
  (https://ncas-cms.github.io/cf-python/method/cf.Field.collapse.html). Note
  well: even the CF reference implementation *drops* multidimensional aux
  coords on reduction — it just does it deliberately, with the 1D collapse
  story intact.
- **Regridding**: ESMF-based `f.regrids` accepts source/destination grids
  whose lat/lon "may be 1-d dimension coordinates or 2-d auxiliary
  coordinates", plus UGRID meshes and DSGs
  (https://ncas-cms.github.io/cf-python/method/cf.Field.regrids.html) — 2D
  coords' primary *consumer* is the regridder.

---

## 6. CF conventions: the normative details

### 6.1 Chapter 5 — coordinate systems

[Quotes verified against ch05.adoc (main) and the released texts
https://cfconventions.org/Data/cf-conventions/cf-conventions-1.12/cf-conventions.html
and -1.11; wording identical unless noted. The unversioned
cfconventions.org/cf-conventions URL serves the working draft.]

- The `coordinates` attribute: "Any longitude, latitude, vertical or time
  coordinate which depends on more than one spatiotemporal dimension **must**
  be identified by the `coordinates` attribute"; its value "is a blank
  separated list of the names of auxiliary coordinate variables", and "There
  is no restriction on the order in which the auxiliary coordinate variables
  appear in the coordinates attribute string." Optionally, dimension
  coordinates too: "it is permissible, but optional, to list coordinate
  variables as well".
- The inverse rule, easy to miss: a multi-valued coordinate that varies in
  only one dimension **and varies independently of other spatiotemporal
  coordinates** "is not permitted to be stored as an auxiliary coordinate
  variable" — it must be the dimension's coordinate variable. The
  "independently" qualifier is what legitimizes DSG's covarying
  `lat(station)`/`lon(station)` (§6.4); a second free-standing 1D labeling of
  an already-labeled axis (ClimaDiagnostics's `date(time)` next to
  `time(time)`) sits in a gray zone CF never really sanctions.
- The subset rule: "The dimensions of an auxiliary coordinate variable must be
  a subset of the dimensions of the variable with which the coordinate is
  associated, with three exceptions" — (1) char-typed string-label coords
  carry a string-length dim; (2) compression by gathering; (3) DSG ragged
  arrays. **Order is not constrained** — only the dimension *set* matters.
- Naming: "It is recommended that the name of a multidimensional coordinate
  variable should not match the name of any of its dimensions because that
  precludes supplying a coordinate variable for the dimension" (rasm.nc
  respects this: `xc(y,x)`, not `x(y,x)`).
- CF's own sanctioned discovery order: "An application that is trying to find
  the latitude coordinate of a variable should always look first to see if any
  of the variable's dimensions correspond to a latitude coordinate variable"
  — dimension coordinates first, then the `coordinates` attribute.
- `axis` may appear on aux coords, but a data variable must not end up with
  two coordinates claiming the same axis letter.
- Missing values: "Missing data is not allowed in coordinate variables";
  "Missing data is allowed in data variables and auxiliary coordinate
  variables" (§2.5.1; the aux allowance dates to CF-1.7). "Generic
  applications should treat the data as missing where any auxiliary
  coordinate variables have missing values."
- When the attribute is absent, CF offers nothing further: the sanctioned
  lookup is dimension→coordinate-variable, then the `coordinates` list
  (quoted above). An unlisted 2D lat/lon is simply not CF-attached; any
  dimension-matching association a tool performs is a heuristic outside the
  standard.

### 6.2 Chapter 4 + Appendix D — parametric vertical coordinates

A parametric vertical coordinate is a (usually dimensionless) 1D coordinate
whose `standard_name` selects a formula from Appendix D and whose
`formula_terms` attribute maps formula terms to file variables
("term: varname …" pairs); the computed quantity's identity can be declared
with `computed_standard_name` (§4.3.3: "To maintain backwards compatibility
with COARDS the use of these attributes is not required, but is strongly
recommended"). If a parametric coordinate has `bounds`, "its boundary
variable must have a formula_terms attribute too" (§7.1). Appendix D rule:
"A term that is omitted from the formula_terms attribute should be assumed
to be zero." The full Appendix D
catalog (https://raw.githubusercontent.com/cf-convention/cf-conventions/main/appd.adoc),
with one-line formulas — note the computed coordinate's rank is the union of
the term ranks, so 1D×(2-3D terms) ⇒ 3-4D result:

| standard_name | formula (computed: air_pressure unless noted) |
|---|---|
| `atmosphere_ln_pressure_coordinate` | `p(k) = p0 * exp(-lev(k))` |
| `atmosphere_sigma_coordinate` | `p = ptop + sigma(k)*(ps(n,j,i)-ptop)` |
| `atmosphere_hybrid_sigma_pressure_coordinate` | `p = a(k)*p0 + b(k)*ps(n,j,i)` (or `ap(k) + b(k)*ps`) |
| `atmosphere_hybrid_height_coordinate` | `z = a(k) + b(k)*orog(n,j,i)` (altitude) |
| `atmosphere_sleve_coordinate` | `z = a(k)*ztop + b1(k)*zsurf1 + b2(k)*zsurf2` (altitude) |
| `ocean_sigma_coordinate` | `z = eta(n,j,i) + sigma(k)*(depth(j,i)+eta(n,j,i))` (depth-family) |
| `ocean_s_coordinate` | `z = eta*(1+s(k)) + depth_c*s(k) + (depth-depth_c)*C(k)` |
| `ocean_s_coordinate_g1` | `z = S(k,j,i) + eta*(1 + S/depth)`, `S = depth_c*s + (depth-depth_c)*C` |
| `ocean_s_coordinate_g2` | `z = eta + (eta+depth)*S(k,j,i)`, `S = (depth_c*s + depth*C)/(depth_c+depth)` |
| `ocean_sigma_z_coordinate` | sigma above `depth_c`, fixed `zlev` below (nsigma levels hybrid) |
| `ocean_double_sigma_coordinate` | two stacked sigma systems split at `f(j,i)` |

(Also `atmosphere_hybrid_sigma_ln_pressure_coordinate` in current drafts.)
The ocean formulas all involve `eta(n,j,i)` — the time-varying free surface —
which is why derived ocean depth is 4D and why every implementation defers
evaluation (§3.3, §4, §2.5).

### 6.3 Chapter 7 — bounds for multidimensional coordinates [ch07.adoc]

A bounds variable adds one trailing (fastest in CDL) vertex dimension: size 2
for 1D coords; for 2D lat(n,m) with four-sided cells, `lat_bnds(n,m,4)` with
"the vertices … ordered such that, when visiting the vertices in order, the
four-sided perimeter of the cell is traversed anticlockwise on the lon-lat
surface as seen from above", starting rule 0=(j-1,i-1), 1=(j-1,i+1),
2=(j+1,i+1), 3=(j+1,i-1); shared vertices of contiguous cells must be
"represented identically in each instance where it occurs". Bounds variables
inherit ("BI") attributes from the parent coordinate and should not restate
them.

### 6.4 Chapter 9 — discrete sampling geometries [ch09.adoc]

DSG files declare a global `featureType` ∈ {point, timeSeries, trajectory,
profile, timeSeriesProfile, trajectoryProfile}. Structure: an **instance
dimension** ("the dimension with subscript i identifies a particular feature
within a collection of features" — e.g. `station`) and **element
dimension(s)** (subscripts o/p "distinguish the data elements that compose a
single feature" — e.g. `obs`). All space-time coordinates become **1D
auxiliary coordinate variables** over these dims — `lat(station)`,
`lon(station)`, `time(station, obs)` — i.e. several coordinates sharing one
axis, none of them a coordinate variable in the NUG sense. Two rules matter
for a reader:

- In DSG files the association is *not* optional: "The `coordinates`
  attribute must be attached to every data variable to indicate the
  spatiotemporal coordinate variables that are needed to geo-locate the
  data." (Contrast §6.1: elsewhere listing is merely how aux coords are
  identified.)
- Feature instances are labeled via `cf_role` ∈ {`timeseries_id`,
  `profile_id`, `trajectory_id`} on an instance variable ("strongly
  recommended", §9.5).
- DSG also inverts the missing-data posture: unused storage requires "the
  data variable and all its auxiliary coordinate variables (spatial and
  time)" to contain missing values (§9.6) — aux coords with missing values
  are structural there, not an anomaly.

Four representations: orthogonal multidimensional (all features share element
coordinates), incomplete multidimensional (padded with missing), contiguous
ragged (count variable), indexed ragged (index variable). The ragged ones are
the third exception to the dims-subset rule (§6.1). For ClimaAnalysis's
scope, the orthogonal/incomplete forms are just "one points axis + several 1D
coordinates over it" — the adopted model covers them with no extra machinery;
ragged forms would require gather/scatter first.

### 6.5 The construct taxonomy

See §5.1 — the auxiliary coordinate / domain ancillary / cell measure / field
ancillary discriminations are the CF data model's (Hassell et al. 2017), and
iris's class hierarchy (§3.1) implements them one-to-one.

---

## 7. What models actually write

### 7.1 WRF — the canonical offender

Raw wrfout files: coordinates are 2D-and-time-varying
`XLAT(Time, south_north, west_east)` / `XLONG(...)` with units
`"degree_north"/"degree_east"` (singular — still CF-legal), plus staggered
variants `XLAT_U/XLONG_U`, `XLAT_V/XLONG_V`
(WRF Registry.EM_COMMON, https://github.com/wrf-model/WRF/blob/master/Registry/Registry.EM_COMMON).
Correcting a common belief: **modern WRF does write a `coordinates`
attribute** — hardwired per stagger as `coordinates = "XLONG XLAT XTIME"`
(`.._U`/`.._V` variants) in
https://github.com/wrf-model/WRF/blob/master/share/wrf_ext_write_field.F,
with the in-code comment "It is a step towards CF". The files are still not
CF: no global `Conventions` attribute, time is a `Times(Time, DateStrLen)`
char array without units/calendar, no `standard_name`s anywhere, no vertical
coordinate in `coordinates`, and grid semantics ride in non-CF attributes
(FieldType, MemoryOrder, stagger); salem's summary: "the coordinate names
are exotic and do not correspond to the dimension names" and times need "a
special parser" (https://salem.readthedocs.io/en/stable/wrf.html).
`xwrf.postprocess()` repairs the rest on the xarray side: decodes times, adds
"CF and COMODO-compliant attributes" (standard_name, grid_mapping),
reconstructs the projection as a `wrf_projection` variable, **computes
regular 1D x/y projection coordinates from the projection metadata** (WRF
grids are regular in projected space — the 2D lat/lon are derived, not
fundamental), destaggers, and computes diagnostics
(https://xwrf.readthedocs.io/en/latest/tutorials/Overview.html,
https://xarray.dev/blog/introducing-xwrf). Lesson: when a curvilinear grid is
secretly "regular grid + known projection", recovering the 1D form beats
supporting 2D coords.

### 7.2 ROMS — parametric vertical + staggered curvilinear, done right

ROMS output is the CF showcase, verified in the model writer itself
(https://github.com/myroms/roms/blob/develop/ROMS/Utility/def_info.F,
def_var.F): `s_rho`/`s_w` get `standard_name = "ocean_s_coordinate_g1"` or
`_g2` keyed on `Vtransform`, `positive="up"`, and
`formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc"`;
`def_var.F` assembles per-Arakawa-C-point `coordinates` strings
(`"lon_rho lat_rho s_rho ocean_time"`, `lon_u/lat_u`, `lon_v/lat_v`,
`lon_psi/lat_psi`) with the 2D `lon_*/lat_*(eta_*, xi_*)` declared in
varinfo.yaml; it even writes an SGRID mesh-topology variable ("SGRID
conventions for staggered data on structured grids"). The xarray ROMS example
loads all of it with default decoding and plots via
`ds.salt.isel(...).plot(x="lon_rho", y="lat_rho")`
(https://docs.xarray.dev/en/stable/examples/ROMS_ocean_model.html). Because
`zeta(ocean_time, eta_rho, xi_rho)` is a formula term, derived depth is 4D
and time-varying — the example computes it manually
(`z_rho = zeta + (zeta + h) * Zo_rho`, dims `(s_rho, xi_rho, eta_rho,
ocean_time)`). This is the hardest generalization target for any
formula_terms implementation.

### 7.3 MPAS and UGRID — unstructured

MPAS writes its Voronoi mesh per its own registry-based spec, not CF:
`latCell/lonCell(nCells)` **in radians**, `latVertex/lonVertex(nVertices)`,
connectivity (`verticesOnCell`, `cellsOnCell`), global `on_a_sphere`
(https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf). The
atmosphere core's Registry.xml contains **zero** `coordinates` attributes
(https://github.com/MPAS-Dev/MPAS-Model/blob/master/src/core_atmosphere/Registry.xml),
so history files carry no CF linkage at all — CDO errors out
("cdfScanVarAttr Variable not found >lons<",
https://code.mpimet.mpg.de/boards/2/topics/13355) and NCAR ships
`convert_mpas` to interpolate to lat-lon first; consumers hard-code the
names. UGRID (now
incorporated by reference into CF ≥1.11) instead declares a mesh-topology
variable (`cf_role = "mesh_topology"`, `node_coordinates`,
`face_node_connectivity`) and data variables point at it via `mesh` +
`location` attributes. uxarray's whole architecture is normalizing
MPAS/SCRIP/Exodus/ESMF/ICON/HEALPix "in the UGRID conventions at the data
loading step" (§8). Data-side, an unstructured variable is just
`var(time, nCells)` — a points axis with several 1D coords over it, the DSG
shape again.

### 7.4 FV3 / cubed sphere (UFS incarnation)

The UFS write component ("quilting", the default) outputs history "not … on
the 6 tiles, but instead as a single global gaussian grid file"; only with
`quilting=.false.` does FV3 write "tiled output in the native projection"
(https://ufs-weather-model.readthedocs.io/en/latest/InputsOutputs.html). The
writer (https://github.com/NOAA-EMC/fv3atm/blob/develop/io/module_write_netcdf.F90)
uses dims `grid_xt`/`grid_yt` (+ `tile` for native cubed-sphere), 1D
dimension variables carrying FMS-style `cartesian_axis = "X"/"Y"` rather
than CF `axis`, separate 2D `lon/lat(grid_yt, grid_xt)` ("T-cell
longitude/latitude"; 3D with `tile` for native output), and — for the native
case — `coordinates = "lon lat"` plus `grid_mapping = "cubed_sphere"` on
every field. "Cubed-sphere dynamical cores produce native-grid output from
an irregular grid with two-dimensional lat-lon coordinates", viewable in
ncview/Panoply via that CF mechanism
(https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/discussions/208);
production pipelines regrid with GFDL's `fregrid` (FRE-NCtools). The
operational pattern — like ClimaAtmos — is to remap to rectilinear lat-lon
at diagnostics time rather than make consumers grid-aware.

### 7.5 CMIP6 — hybrid levels, exact CDL

Authoritative source: the CMOR tables — `CMIP6_coordinate.json` entry
`standard_hybrid_sigma` (out_name `lev`, formula `p = a*p0 + b*ps`,
`z_factors = "p0: p0 a: a b: b ps: ps"`, bounds mapping via a_bnds/b_bnds)
and the `alternate_hybrid_sigma` variant (`p = ap + b*ps`, ap in Pa), with
term variables defined in `CMIP6_formula_terms.json`
(https://github.com/PCMDI/cmip6-cmor-tables/blob/master/Tables/CMIP6_coordinate.json).
As emitted, from a published CESM2 CMIP6 header
(https://www.wdc-climate.de/ui/header?acronym=C6CMNRCES2amr10111AEzggn00226):

```
double lev(lev):  standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
                  units = "1", positive = "down", axis = "Z"
                  formula = "p = a*p0 + b*ps"
                  formula_terms = "p0: p0 a: a b: b ps: ps", bounds = "lev_bnds"
double lev_bnds(lev, nbnd): formula_terms = "p0: p0 a: a_bnds b: b_bnds ps: ps"
double a(lev), b(lev), a_bnds(lev,nbnd), b_bnds(lev,nbnd); float p0, ps(time,lat,lon)
zg:coordinates = "time lev lat lon"
```

Note (a) **bounds carry their own formula_terms** (a_bnds/b_bnds) — bounds of
a parametric coordinate are computed, not stored; (b) some models use the
`ap: ap b: b ps: ps` variant (ap = a*p0 pre-multiplied); (c) CMOR lists even
plain dimension coordinates in `coordinates` — the "permissible, but
optional" listing of §6.1 in the wild.

### 7.6 Consumers vs files that omit the `coordinates` attribute

- **xarray**: association is purely attribute-driven (§2.3), so unreferenced
  aux coords stay data variables and plotting/selection silently degrade to
  index space; users must `set_coords` by hand. The WRF/MPAS contrast makes
  the point: raw WRF `XLAT/XLONG` *do* become coordinates (WRF hardwires the
  attribute, §7.1) while raw MPAS `latCell` does not.
- **iris**: unassociated aux coords load as *separate cubes* [verified §7.7].
- **netCDF-Java**: the only mainstream reader that recovers gracefully — if
  no convention matches (`Conventions` attr, then registered builders'
  `isMine()`), it falls back to `DefaultConventions`, which types axes by
  COARDS units (`degrees_east` family → Lon, …) and variable names
  (`lat`/`latitude`, …) and assembles implicit coordinate systems from
  dimension matching
  (https://github.com/Unidata/netcdf-java/blob/develop/cdm-core/src/main/java/ucar/nc2/internal/dataset/conv/DefaultConventions.java).
- **CDO**: reconstructs a `griddes` from CF metadata; without it, variables
  fall back to "generic" grids — "Data variables without or wrong grid
  description are called generic grids" — and most spatial operators refuse
  (`sellonlatbox: Unsupported grid type: generic`). Documented remedies:
  `setgrid,<gridfile>` with a hand-written griddes, and (ICON-specific)
  auto-fetching the external grid file named by the global `grid_file_uri`
  attribute (CDO User Guide §1.6; https://code.mpimet.mpg.de/boards/1/topics/752).
  CDO also "does not handle time-dependent coordinate variables" — WRF's
  `XLONG(Time,…)` must be de-time-dimensioned first
  (https://code.mpimet.mpg.de/boards/1/topics/10771).
- **Panoply/ncview**: georeference via the CF 2D-coordinates mechanism when
  the attribute exists (per the GFDL discussion in §7.4). NASA GES DISC
  states the failure mode: "If a variable is not CF-compliant, Panoply
  displays the metadata, but an image may not be displayed"
  (https://disc.gsfc.nasa.gov/information/howto?title=Quick+View+Data+with+Panoply);
  Panoply's own help pages were offline during this survey, so its
  WRF-specific handling could not be verified.

The pattern: omitting `coordinates` downgrades every consumer except
netCDF-Java to index space or manual repair — per-model shim packages (xwrf,
xroms/salem, uxarray, xgcm) exist largely to re-add this metadata.

### 7.7 ClimaDiagnostics (the local problem statement) [verified in source]

`/home/kphan/.julia/packages/ClimaDiagnostics/yUbF5/src/netcdf_writer.jl:718-723` writes
data variables with only `short_name/long_name/units/comments/start_date` —
**no `coordinates` attribute**. With topography,
`netcdf_writer_coordinates.jl:705-738` emits `z_reference` (1D dim coord) plus
`z_physical(z_reference, lat, lon)` carrying only `units="m"` (no
standard_name, no axis, no positive), with the comment (line 736): "We do not
output this name because it is not an axis." Consequences, verified on
`ta_with_topography.nc`:

- xarray: `z_physical` and `date` load as *data variables* even with
  `decode_coords="all"` (only `time_bnds`/`date_bnds` promote, via `bounds`).
- iris: `z_physical` loads as a separate 3D cube; `ta` has no vertical aux.
- cf_xarray: axes X/Y/Z/T resolve (the 1D dims carry `axis` attrs) but
  `ds.cf["vertical"]` KeyErrors; nothing links `z_physical` to `ta`.

So for Clima output the discovery chain **must** include a name-convention
tier; no standards-based tool will ever associate `z_physical` as written.

---

## 8. Other consumers in brief

- **CDO** (User Guide 2.6.3, §1.6): grids are first-class values — "CDO
  supports structured grids like regular lon/lat or curvilinear grids and
  also unstructured grids". A grid description (`cdo griddes`) has `gridtype`
  ∈ {lonlat, gaussian, projection, curvilinear, unstructured} with
  full-length `xvals/yvals` for curvilinear/unstructured and
  `nvertex`+`xbounds/ybounds` for corners ("optional if area weights are not
  needed"). Remapping capability is grid-type-dependent: `remapbil`/`remapbic`
  "only works on quadrilateral curvilinear source grids"; `remapcon` (YAC)
  is "completely general … for any grid on a sphere"; `remapnn`/`remapdis`
  are spherical nearest-neighbor for anything. Selection = remap, never
  label-lookup (https://code.mpimet.mpg.de/projects/cdo/embedded/cdo.pdf).
- **NCO**: hyperslabbing is dimension-centric (`-d dim,min,max`). The
  `-X`/`--auxiliary lon_min,lon_max,lat_min,lat_max` switch adds value-based
  selection but **only for 1-D cell-based (unstructured) grids**, finding
  lat/lon via `standard_name` through the `coordinates` attribute; the manual
  is explicit: "This feature cannot be used to select regions of 2D grids
  (instead use the ncap2 where statement for such grids)". Regridding is
  external weight files (`ncremap`, SCRIP/ESMF/TempestRemap)
  (https://nco.sourceforge.net/nco.html#Auxiliary-Coordinates).
- **GDAL**: RFC 4 "geolocation arrays" — per-pixel (or subsampled) lon/lat
  arrays in the `GEOLOCATION` metadata domain (`X_DATASET`, `Y_DATASET`,
  `PIXEL/LINE_OFFSET/STEP`), consumed by `gdalwarp -geoloc`. The netCDF
  driver produces an affine geotransform only from equally-spaced 1D coord
  vars; otherwise it emits GEOLOCATION metadata. The inverse mapping
  (geoloc→pixel) is a backmap/quadtree search, not an index
  (https://gdal.org/en/stable/development/rfc/rfc4_geolocate.html,
  https://gdal.org/en/stable/drivers/raster/netcdf.html).
- **R stars**: dimension objects may hold "a matrix with longitudes or
  latitudes for all cells (in case of curvilinear grids)" — i.e. 2D coord
  matrices stored *as the dimension's values*, flagged `curvilinear`.
  Supported: plotting (point or cell polygons), `st_downsample`, and
  `st_warp` to a regular grid; selection/aggregation on curvilinear objects is
  not offered — warping is the sanctioned escape hatch. Discovery is delegated
  to GDAL's `GEOLOCATION` metadata domain (`X_DATASET`/`Y_DATASET` pointers)
  (https://r-spatial.github.io/stars/articles/stars4.html).
- **pyresample**: `SwathDefinition` (explicit 2D lon/lat arrays) vs
  `AreaDefinition` (CRS + extent + shape). The docs state the cost model
  plainly: with a swath, "any information we want about the swath will either
  require looking at every coordinate or it will make some assumption"
  (https://pyresample.readthedocs.io/en/latest/concepts/geometries.html).
  Resampling between geometries is KDTree nearest-neighbor / bilinear.
- **uxarray**: the exception to "no grid object": `UxDataset`/`UxDataArray`
  wrap xarray and carry a sidecar `Grid` ("Encapsulates a xarray.Dataset for
  storing the grid definition"), normalizing MPAS/CAM-SE/ICON/ESMF/GEOS/
  HEALPix/SCRIP/Exodus **into UGRID conventions at load time**, then provides
  grid-aware operators (remapping, subsetting, cross-sections, zonal
  averaging) that plain xarray cannot express
  (https://uxarray.readthedocs.io/en/latest/getting-started/overview.html).
  This is what "supporting unstructured for real" costs: a whole package.

---

## 9. Julia ecosystem: verified state and gaps

### 9.1 CommonDataModel.jl / NCDatasets.jl

`CommonDataModel.coord(v, standard_name)` (v0.4.4,
`/home/kphan/.julia/packages/CommonDataModel/WcahY/src/cfconventions.jl:130-179`) —
the only CF-ish discovery in Julia — **never reads the `coordinates`
attribute** (no reference anywhere in the source). Its chain:

1. scan *all* dataset variables for matching `standard_name` whose
   `dimnames ⊆ dimnames(v)` (lines 147-151);
2. fall back to units regexes — time `r".*since.*"`, the eight
   degrees_east/north spellings — again requiring the dims-subset rule, and
   preferring the candidate with the **most** dimensions: "prefer e.g. vectors
   over scalars; this is necessary for ROMS model output" (lines 154-176).

So it happily returns 2D `lon_rho(eta,xi)` — discovery handles aux coords,
by heuristic rather than by the CF attribute. Meanwhile
`@select`/`select` (src/select.jl:335-394) builds its coordinate list from
*any 1D variable whose dim is shared with the variable* ("Any 1D variable with
the same dimension name can be used in `@select`", docstring line 182 — this
covers DSG station files) but `coordinate_value` hard-asserts `ndims(ncv)==1`
(line 364): 2D coordinates cannot participate in selection. `bounds()`
(cfconventions.jl:192) follows the `bounds` attribute;
`ancillaryvariables()` (line 20) follows `ancillary_variables`. No
formula_terms support anywhere.

### 9.2 DimensionalData.jl (v0.29.27)

Strictly one `Lookup` per dimension. The multidimensional offerings solve
*different* problems:

- `MergedLookup <: MultiDimensionalLookup <: Lookup{T,1}`
  (src/Dimensions/merged.jl): a **1D** vector of coordinate *tuples* over one
  merged dimension (`mergedims(da, (X,Y) => :space)`) — the
  flatten-to-points model. Selection is a linear scan (`At`, `Between`,
  `Where`); `Near` **throws** `ArgumentError("Near is not implemented for
  coordinates")` (merged.jl:127); `reducelookup(::MergedLookup) =
  NoLookup(OneTo(1))`. The source carries "TODO: Do we still need `Coord` as
  a dimension?" (line 130) — experimental status in-code.
- `Transformed <: Unaligned` (src/Lookups/lookup_arrays.jl:539+): coordinates
  defined *functionally* by an affine transform of the index space
  (CoordinateTransformations.jl) — handles rotated/affine grids exactly, but
  cannot represent measured curvilinear arrays.
- Reduction semantics (contrast with both xarray and iris): Julia reductions
  keep singleton dims, and DD computes a representative length-1 lookup —
  `reducelookup(::AbstractSampled)` yields the span center for `Regular`,
  keeps outer bounds for `Irregular` (lookup_arrays.jl:778-800). Closer to
  iris than to xarray, but with no bounds provenance on aux coords (there are
  none).

Unstructured/curvilinear support was requested in
https://github.com/rafaqz/DimensionalData.jl/issues/1105; maintainer: "It
doesn't currently, but we have been discussing something like this - more
generally we need an abstraction without the assumption of a structured grid
or being part of julias `AbstractArray` interface."

### 9.3 Rasters.jl (v0.14.7)

Dimension construction is name-based only
(src/sources/commondatamodel.jl:218-233): a dim gets a lookup iff a
*same-named* variable exists, else `NoLookup(OneTo(len))`. The `coordinates`
attribute is never parsed; only `grid_mapping` is copied into metadata
(lines 190-199). A curvilinear file therefore loads with bare index dims and
`xc`/`yc` as unrelated layers. Related activity:

- https://github.com/rafaqz/Rasters.jl/pull/854 "WIP: geometry lookup" —
  xvec-style vector data cubes (a dimension whose lookup is a geometry
  vector), with CF DSG point/line/polygon I/O on the TODO list. Aimed at
  point/geometry axes, not curvilinear rasters.
- https://github.com/rafaqz/Rasters.jl/issues/342 — user asking for 2D
  coordinate arrays of a raster; answer: `DimPoints` (lazy tuple
  materialization of *rectilinear* coords).

### 9.4 YAXArrays.jl, GeoMakie/Makie

- YAXArrays inherits DD's model; uxarray-style support was requested in
  https://github.com/JuliaDataCubes/YAXArrays.jl/issues/551 and deferred:
  "this would happen on the DimensionalData level".
- Plotting is not the blocker: Makie's `surface`/`mesh` accept 2D coordinate
  matrices, and `heatmap`-with-curvilinear-x/y is a long-open request
  (https://github.com/MakieOrg/Makie.jl/issues/742, open since 2020, 44
  comments; recommended workaround `surface(X, Y, zeros(size(Z)); color=Z,
  shading=NoShading)`). CliMA-adjacent:
  https://github.com/CliMA/Oceananigans.jl/issues/5651 (curvilinear-mesh
  plotting for terrain-following/spherical grids). GeoMakie's source has no
  curvilinear-specific code (checked v across `/home/kphan/.julia/packages/GeoMakie/`).
- GitHub-wide searches for `formula_terms` and "auxiliary coordinate" across
  Julia repos return essentially nothing data-model-related — the concepts
  have no Julia implementation to date.

---

## 10. What nobody does

1. **Label/orthogonal indexing on N-D coordinates.** `sel(lat=slice(...))` on
   2D lat exists in no surveyed system. xarray: open feature request (#10572);
   iris: dedicated `CoordinateMultiDimError`; CDO/pyresample: nearest-neighbor
   remap instead; Julia: N/A. The value of an N-D coordinate at a point is
   well-defined; its *inverse* is not (non-monotonic, non-injective, no
   ordering) — every system that tried stopped at nearest-neighbor trees.
2. **Reducing an aux coordinate along one of several spanned dims** as a
   default behavior. iris collapses only whole coords (and then collapses all
   dims they span [verified]); xarray drops. Nobody computes "mean lat per
   remaining column" implicitly — it is ill-defined without weights.
3. **Automatic formula_terms evaluation on load.** Everyone makes derived
   vertical coordinates opt-in (iris factories are built at load but lazy;
   cf_xarray requires an explicit `decode_vertical_coords()`; netCDF-Java
   builds transforms on request; xarray core: nothing). Time-varying terms
   (ROMS `zeta`) make the derived coord 4D — materializing eagerly is
   prohibitive.
4. **Guessing associations beyond dims-subset + identity attributes.** No
   tool invents a coordinate link from shape alone without at least a units/
   standard_name/axis/name signal (CommonDataModel's scan is the most
   aggressive surveyed, and it still requires standard_name or CF units).
5. **A grid object inside the variable.** Only uxarray has a grid object at
   all, and it lives beside the data, not in it. The CF data model made the
   choice explicit: coordinate systems are emergent subsets of the attached
   coordinate constructs, and the Domain "is not a construct of the data
   model … but an abstract concept" (Hassell et al. 2017, §5.1 above) — the
   flat coordinate-list model won everywhere.

---

## 11. Implications for ClimaAnalysis

The adopted model — `coordinate = (name, spans ⊆ var dims, values,
attributes)`, dimension coordinate ⇔ `spans == (name,)` — is exactly the
convergent model (xarray coords, iris AuxCoord+coord_dims, CF ch5 subset
rule). Decisions the survey settles:

1. **Propagation semantics to copy** (per-operation):
   - *Slice along a spanned axis*: slice the coordinate along that axis
     (universal).
   - *Integer index / window that drops an axis*: keep the coordinate with
     that axis removed; if all spanned axes drop, keep a **scalar coordinate**
     (xarray + iris agree [verified]) — it preserves "where this slice was
     taken" for titles, provenance, and future concat/merge.
   - *Reduce over a spanned axis*: **do not silently drop.** xarray's drop is
     its most complained-about aux-coord behavior (#8317, #9168, #3510).
     Iris collapses 1D coords into range-bounded scalars; and notably even
     cf-python, the CF reference implementation, *removes* ≥2D aux coords on
     collapse — but deliberately, alongside an intact 1D-collapse story. So
     the defensible spec is: 1D coords over the reduced axis → scalar summary
     (or drop with provenance); multi-span coords touching the reduced axis →
     remove, but record it (a `reduced_over` provenance attribute or a
     warning). Never keep stale values.
   - *Reduce over an unspanned axis*: keep unchanged (universal).
   - *Concat/cat along a spanned axis*: concatenate coordinate values along
     it. *Along an unspanned axis*: require equality (iris `ConcatenateError`
     is the right strictness for a curated package; xarray's default
     "broadcast the conflict into a new dim" is a documented trap).
   - *Transpose/permute*: either physically permute values (xarray) or remap
     spans (iris). With a `spans` tuple stored per coordinate, the iris
     approach is free and copy-less.
   - *Arithmetic between vars*: require equal aux coords on shared spans;
     error (not silent drop) on conflict — silent drop is xarray's second
     most-cited wart.
2. **Discovery chain to implement** (each tier only fills gaps left by the
   previous; every tier enforces `spans ⊆ dims(var)` — CF's one hard rule):
   1. dimension coordinates by name==dim (CF's "always look first" rule);
   2. the variable's `coordinates` attribute (split on blanks; ignore missing
      names silently, as xarray does);
   3. dataset-level heuristics à la CommonDataModel/cf_xarray:
      `standard_name`, `axis`, `positive`, CF units regexes — preferring
      higher-rank matches (the ROMS lesson in CommonDataModel);
   4. name-convention table (`z_physical`↔`z_reference`, `xc`/`yc`,
      `date(time)`) — **mandatory**, because ClimaDiagnostics omits the
      attribute and no standard tool can associate `z_physical` [verified
      §7.7]. Also promote `bounds` targets (chapter 7) but as bounds, not
      coordinates.
   Record which tier fired (an `origin` field/attr) for debuggability.
3. **Selection: explicitly do not support** label selection on multi-span
   coordinates. Error with the iris message shape ("`z_physical` spans
   (z_reference, lat, lon), not just itself — interpolate/resample first"),
   and offer the two sanctioned escapes: (a) mask/`where`-style filtering,
   (b) nearest-neighbor as a *separate, explicit* API (KDTree over selected
   coords, the xoak/NDPointIndex/pyresample pattern) if demand appears. Note
   the xarray lesson [verified]: even NDPointIndex is nearest-only and does
   not survive slicing — a bolted-on index is not part of the data model, so
   ClimaAnalysis loses nothing by deferring it.
4. **Keep formula_terms as a recipe, not data** (phase 2, as scoped): the CF
   data model is emphatic that formula terms are *domain ancillaries*, not
   coordinates (§5.1), and all three engines that implement them defer
   evaluation (iris lazy factories, netCDF-Java `makeVerticalTransform()` on
   demand per time step, cf_xarray's explicit `decode_vertical_coords()`).
   Store the parametric standard_name + term→variable map; provide an
   explicit materialize verb. The ROMS time-varying case (4D derived depth
   via moving `zeta`) is the reason not to compute eagerly. Design the
   attachment so a lazy `values` can slot into the same coordinate struct
   later (iris proves factories can behave exactly like stored coords under
   slicing, §3.3).
5. **Writer-side fix**: ClimaAnalysis should *write* `coordinates =
   "z_physical"` (and standard_name/positive on `z_physical`) when saving,
   and ClimaDiagnostics should be nudged to do the same — one attribute makes
   Clima output legible to xarray/iris/CDO forever. xarray's encoding
   round-trip (§2.3) shows this costs nothing.
6. **Bounds**: adopt CF's trailing-vertex-dimension layout `(dims..., nv)`
   internally (n×2 for 1D as in API.md; (n,m,4) reserved for future 2D) so
   files round-trip without reshaping.
7. **Scalar coordinates** (0 remaining spans) earn their keep in labels/
   provenance; they are cheap (a 0-d value + attrs) and both reference
   implementations keep them. Support them from day one rather than
   special-casing "coordinate must have ≥1 span".

---

*Sources are cited inline. Local verification scripts: `xr_semantics.py`,
`xr_part2.py`, `iris_semantics.py`, `iris_part2.py`, `iris_part3.py`,
`topo_consumers.py` in the session scratchpad (venv with xarray 2026.7.0,
iris 3.16.0, cf_xarray). Local sources cited by absolute path + line number.*
