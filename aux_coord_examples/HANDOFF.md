# Handoff: Auxiliary coordinates in NetCDF output and how to model them in ClimaAnalysis

**Context:** Kevin (kphan2@caltech.edu) is designing support in **ClimaAnalysis.jl** for grids whose
coordinates are not simple 1D per-axis vectors — i.e. CF *auxiliary coordinates* (topography-warped
heights, curvilinear 2D lat/lon, parametric vertical coordinates). This document summarizes an
investigation into how CliMA packages and the wider ecosystem produce/consume such files, plus a
proposed data model. Example files live in this directory (`/home/kphan2/aux_coord_examples/`).

---

## 1. Files in this directory

| file | what it is |
|---|---|
| `ta_with_topography.nc` | Generated with ClimaDiagnostics.jl `main` (cubed sphere, `Hypsography.LinearAdaption`, `LevelsMethod()`). Contains `z_reference(z_reference)` (1D dimension coord, nominal levels 500–9500 m) **and** `z_physical(z_reference, lat, lon)` (3D auxiliary variable, actual altitudes; lowest level ranges 652–3342 m across the domain). The data variable `ta` was deliberately set equal to physical altitude, and was verified to match `z_physical` exactly. |
| `rasm.nc` | Downloaded from github.com/pydata/xarray-data. Canonical **curvilinear grid** example: `Tair(time, y, x)` with `Tair:coordinates = "yc xc"`, where `xc(y,x)`/`yc(y,x)` are 2D lon/lat auxiliary coordinates and `x`/`y` are **bare index dimensions with no coordinate variables at all**. |
| `coord_model_demo.jl` | Runnable ~90-line Julia sketch of the proposed data model (needs only NCDatasets). Loads both files above and demonstrates the generic coordinate-lookup contract. Output reproduced in §6. |

The generation script for `ta_with_topography.nc` built an `ExtrudedFiniteDifferenceSpace`
(EquiangularCubedSphere, GLL{4}, 10 vertical levels 0–10 km, `z_surface = 1000*(cosd(lat)+cosd(lon)+1)`),
then used `Writers.NetCDFWriter(space, dir; num_points=(36,18,10), z_sampling_method=Writers.LevelsMethod(),
start_date=...)` + `ScheduledDiagnostic` + `interpolate_field!`/`write_field!` — mirroring
ClimaDiagnostics `test/writers.jl` (the hypsography test only exists on the `pfull-coords` branch).

## 2. What CliMA's writers actually produce (verified in source)

Source examined: `/home/kphan2/worktree/ClimaDiagnostics.jl/{main,pfull-coords}` and
`/home/kphan2/worktree/ClimaAtmos.jl/main`.

Four vertical-coordinate flavors exist:

1. **Flat hypsography + `LevelsMethod`** → plain 1D `z` dimension coordinate (units m, axis Z).
   This is what the ClimaLand monthly files the user has (er/et/gpp/... on 1°×1°, 15 soil levels) look like.
2. **Topography + `LevelsMethod`** → 1D `z_reference` dimension coordinate **plus** a free-floating
   `z_physical(<hdims>..., z_reference)` variable, units "m", written **once** at file creation
   (time-independent), from the interpolated coordinate field.
   Implementation: `add_space_coordinates_maybe!` for `FiniteDifferenceSpace` with
   `interpolated_physical_z` in `src/netcdf_writer_coordinates.jl` (main), around lines 541–587;
   dispatch decided by `Spaces.grid(space).hypsography isa Grids.Flat` around line 501.
   **Crucial:** ClimaDiagnostics does **NOT** set a CF `coordinates = "z_physical"` attribute on the
   data variables — comment at `netcdf_writer_coordinates.jl:584`: "We do not output this name because
   it is not an axis." So `z_physical` is undiscoverable by CF-aware tooling; consumers must know the
   name. ClimaAnalysis currently just lists it in `ALTITUDE_NAMES = ["z", "z_reference", "z_physical",
   "height"]` (`src/outvar_dimensions.jl:6`).
3. **`FakePressureLevelsMethod`** → still a plain 1D `z` in meters, exponentially spaced; metadata
   indistinguishable from case 1.
4. **`RealPressureLevelsMethod`** (ClimaDiagnostics `pfull-coords` branch; wired into ClimaAtmos at
   `src/simulation/AtmosSimulations.jl:92` behind a `pressure_coordinates` diagnostic spec flag) →
   fields interpolated onto fixed pressure levels (default `era5_pressure_levels()`); writes a 1D
   `pressure_level` **dimension coordinate** with `units="Pa"`, `standard_name="air_pressure"`,
   `stored_direction="increasing"`. Not an auxiliary coordinate. Uses a `PfullInterpolator` built
   from a `compute_pfull!` callback.

Also: ClimaAtmos writes surface elevation as a **separate file** `orog_inst.nc`
(`orog(lat, lon, time)`, long_name "Surface Altitude") rather than as an aux variable inside 3D files.
Existing runs on disk (e.g. `ClimaAtmos.jl/main/output/amip_target_diagedmf/`) used
`FakePressureLevelsMethod`, so **no pre-existing local file had `z_physical`** — that's why
`ta_with_topography.nc` was generated.

Other quirks noted in ClimaLand output (from ncdump inspection): data variables stored as
`(lat, lon, time)` / `(z, lat, lon, time)` with **time as the fastest-varying (last-in-C-order)
dimension** despite time being UNLIMITED (unconventional; CF tools expect time-first); `time:units="s"`
with no epoch (calendar info only in a separate `date` variable with "seconds since <start>"); some
unit bugs (`sie` claims W m^-2 for an energy density) and empty units strings for dimensionless vars.

## 3. Additional online examples worth testing against

- `rasm.nc` (here) — curvilinear, proper CF `coordinates` attribute.
- `ROMS_example.nc` (same xarray-data repo) — CF **parametric vertical coordinate**: `s_rho` with
  `standard_name="ocean_s_coordinate_g2"` and `formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h
  depth_c: hc"`; physical depth is **4D and time-varying** (must be *computed*, since `zeta` is
  the moving sea surface). This is the hardest generalization target.
- Any CMIP6 CESM atmosphere file — hybrid sigma-pressure: `lev` with
  `formula_terms = "a: hyam b: hybm p0: P0 ps: PS"`, requiring 3D surface pressure `PS(time,lat,lon)`.
- CF conventions **Appendix D** catalogs all formula-terms vertical coordinate types.

## 4. How other ecosystems handle this (survey)

**Python**
- **xarray**: coords are DataArrays with dims ⊆ variable dims; "non-dimension coordinates" = aux
  coords, any rank. `decode_coords=True` promotes variables named in the CF `coordinates` attribute
  ( `"all"` also handles `bounds`, `grid_mapping`). Aux coords propagate through ops and drive 2D
  plotting, but **`.sel` only works on dimension coordinates**; 2D selection needs `.where` masks,
  flexible/custom indexes, or **xoak** (KD-tree nearest-neighbor over 2D lat/lon). No formula_terms
  evaluation in core.
- **cf_xarray**: name resolution by `standard_name`/`units`/`axis` (`ds.cf["latitude"]`) and
  `decode_vertical_coords()` which **computes** physical z from `formula_terms` — the direct analog
  of materializing `z_physical`.
- **iris** (Met Office): most CF-faithful model. `Cube` holds `DimCoord`/`AuxCoord` as first-class
  objects each knowing which data dims they span, plus **`AuxCoordFactory`** ("derived coordinates"):
  e.g. `HybridHeightFactory` lazily combines 1D level coord + 2D orography aux coord into 3D altitude.
- **cfdm / cf-python**: complete CF data-model implementation (aux coords, domain ancillaries, cell
  measures). **MetPy**: grid_mapping → CRS. **xgcm**: staggered grids. **UXarray**: UGRID unstructured.
  **xESMF/ESMF**: regridding consumes 2D coords + bounds.

**Java/C/C++/R**
- **netCDF-Java CDM** (THREDDS): explicit 3-layer architecture — raw data access → *coordinate
  systems layer* (`CoordinateAxis` may be multidimensional; `CoordinateSystem` assembled from CF
  attrs + heuristics) → feature types; vertical transforms compute 3D pressure/height from
  formula_terms on demand. Best architectural reference for this problem.
- **CDO**: curvilinear grids are a first-class grid type (2D lon/lat + bounds; `cdo griddes`);
  SCRIP remapping uses them; understands hybrid level coefficients. **NCO**: dimension-based,
  delegates regridding. **GDAL**: "geolocation arrays" (2D lon/lat driving warps). **VTK/ParaView**:
  CF reader converts 2D coords / some vertical transforms into curvilinear structured grids.
  **R stars**: supports curvilinear (2D coordinate matrices); terra does not.

**Julia (verified locally)**
- **NCDatasets.jl / CommonDataModel.jl**: `coord(ncv, "latitude")` finds coordinates by
  `standard_name`, falling back to units regexes, requiring candidate dims ⊆ variable dims
  (`~/.julia/packages/CommonDataModel/1dfiZ/src/cfconventions.jl:130`; even has a ROMS special case).
  Does **not** parse the CF `coordinates` attribute; no formula_terms support. `@select` is 1D-only.
- **DimensionalData.jl / Rasters.jl / YAXArrays.jl**: strictly one 1D lookup per axis; curvilinear
  geolocation unsupported (Rasters has CRS but not 2D coords; DD multidim lookups experimental at best).
- **Conclusion: no Julia package handles auxiliary/parametric coordinates well. ClimaAnalysis would
  be the first.**

## 5. Proposed data model (the convergent design)

> **A coordinate = (name, dims ⊆ var.dims, values, attribs).**
> A *dimension coordinate* is not a separate concept — it's the special case `dims == (name,)`.

This is literally the definition used by xarray (coord DataArrays), iris (`AuxCoord` + `coord_dims`),
netCDF-Java (`CoordinateAxis`), and CF itself (the `coordinates` attribute may only name variables
whose dims ⊆ the data variable's dims). It uniformly covers:

| example | name | dims |
|---|---|---|
| ClimaLand monthly files | `lat` | `(lat,)` |
| `ta_with_topography.nc` | `z_reference` | `(z_reference,)` |
| `ta_with_topography.nc` | `z_physical` | `(z_reference, lat, lon)` |
| `rasm.nc` | `yc` | `(y, x)` |
| ROMS physical depth | `z_rho` | `(time, s_rho, eta, xi)` |

Julia sketch (full runnable version: `coord_model_demo.jl` in this directory):

```julia
struct Coordinate{T, N, A <: AbstractArray{T, N}}
    name::String
    dims::NTuple{N, String}   # which axes of the PARENT VARIABLE it spans
    values::A                 # array of exactly that shape
    attribs::Dict{String, Any}
end

struct Var{T, N, A}
    name::String
    dims::NTuple{N, String}
    data::A
    coords::Dict{String, Coordinate}
    attribs::Dict{String, Any}
end

is_dimension_coord(c) = c.dims == (c.name,)

# Core lookup contract: full data index -> this coordinate's value there,
# regardless of how many axes the coordinate spans.
function coordvalue(var, cname, I::Tuple)
    c = var.coords[cname]
    axpos = map(d -> findfirst(==(d), var.dims), c.dims)
    c.values[map(p -> I[p], axpos)...]
end
```

Discovery in the loader: promote every other variable in the file whose dims ⊆ var dims
(`issubset(cdims, vdims)`). Guard label-selection: `nearest_index` errors informatively
("spans (...), not just itself — interpolate/regrid first") for non-dimension coords.

### Design principles distilled from the survey

1. **Separate "variables in the file" from "the interpreted grid"** (netCDF-Java, iris): a decode
   step produces coordinate objects; don't overload the axis/values model.
2. **Discovery is a fallback chain**: CF `coordinates` attribute → heuristics
   (`standard_name`/`axis`/units) → **name conventions** (`z_reference`+`z_physical` pairing).
   The last tier is mandatory because ClimaDiagnostics omits the `coordinates` attribute.
3. **Derived coordinates as lazy factories** (iris `AuxCoordFactory`, cf_xarray
   `decode_vertical_coords`): store the recipe (formula_terms, or "z_reference + orography"),
   compute the multidim array on demand. Bolts on by making `values` lazy; struct unchanged.
4. **Don't promise indexing on aux coords.** Every mature package punts: nearest-neighbor trees,
   masks, or "regrid to rectilinear first." Keep `slice`/`window` on dimension coordinates only,
   with explicit transforms (`to_pressure_coordinates`, `interpolate_to_levels`) for the rest.
5. Semantic queries become "give me the coordinate whose standard/short name is X" rather than
   "which axis is altitude" — the answer may span 1 axis or 3.

## 6. Demo output (proof the model works on the real files)

```
ta dims: ("time", "lon", "lat", "z_reference")  size: (1, 36, 18, 10)
date         spans ("time",)                     -> auxiliary
lat          spans ("lat",)                      -> dimension
lon          spans ("lon",)                      -> dimension
time         spans ("time",)                     -> dimension
z_physical   spans ("lon", "lat", "z_reference") -> auxiliary
z_reference  spans ("z_reference",)              -> dimension

point at index (1, 10, 5, 2):
  lon = -87.4   lat = -47.6   z_reference = 1500.0   z_physical = 2960.8

nearest z_reference to 3000 m -> index 3
z_physical selection -> ERROR: cannot select by 'z_physical': it spans
  ("lon", "lat", "z_reference"), not just itself. Interpolate/regrid first.

Tair dims: ("x", "y", "time")  size: (275, 205, 36)
time spans ("time",) -> dimension
xc   spans ("x","y") -> auxiliary
yc   spans ("x","y") -> auxiliary
point at index (100, 150, 3): lon = 118.33, lat = 57.19, Tair = -25.55
```

Notable emergent behavior: `time_bnds`/`date_bnds` were auto-excluded (they carry the `nv` dim,
absent from `ta`) — the subset rule does the filtering; and `date(time)` was correctly classified
as an *auxiliary* coordinate (a second 1D labeling of the time axis, named differently than its dim),
a case not explicitly designed for. The z_reference=1500 vs z_physical=2960.8 contrast at one column
is the terrain warping, answered through a single generic code path.

## 7. Relevant local paths

- ClimaDiagnostics worktrees: `/home/kphan2/worktree/ClimaDiagnostics.jl/{main,pfull-coords,itime,thread}`
  — topography writer: `main/src/netcdf_writer_coordinates.jl` (~lines 475–587);
  pressure coords: `pfull-coords/src/netcdf_writer_coordinates.jl` (RealPressureLevelsMethod,
  `pressure_level` dim written ~lines 347–380); hypsography test: `pfull-coords/test/writers.jl:411`,
  space setup `pfull-coords/test/TestTools.jl:89`.
- ClimaAtmos: `/home/kphan2/worktree/ClimaAtmos.jl/main/src/simulation/AtmosSimulations.jl:42-110`
  (writer wiring, incl. pressure writer).
- ClimaAnalysis worktrees: `/home/kphan2/worktree/ClimaAnalysis.jl/*` — current altitude handling
  in `src/outvar_dimensions.jl:6`.
- CommonDataModel heuristic: `~/.julia/packages/CommonDataModel/1dfiZ/src/cfconventions.jl:130`.
