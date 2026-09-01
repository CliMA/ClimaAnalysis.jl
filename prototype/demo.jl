# Demo of the axes-vs-coordinates prototype against two real files:
#   - aux_coord_examples/ta_with_topography.nc  (ClimaDiagnostics, LinearAdaption
#     hypsography: 1D z_reference dim + 3D z_physical aux coord, no `coordinates`
#     attribute, data deliberately equal to z_physical)
#   - aux_coord_examples/rasm.nc                (canonical curvilinear grid:
#     2D xc/yc aux coords over bare x/y index dims)
# The .nc files are not committed; (re)create them first with
#   include("aux_coord_examples/generate_examples.jl")
#
# Run from the repo root, in a fresh session:  include("prototype/demo.jl")

include(joinpath(@__DIR__, "aux_grid.jl"))
using .AuxGrid
using NCDatasets
using .AuxGrid: coord   # explicit: NCDatasets also exports a `coord`

const EXAMPLES = joinpath(@__DIR__, "..", "aux_coord_examples")

npass = Ref(0)
nfail = Ref(0)
function check(ok, msg)
    ok ? (npass[] += 1) : (nfail[] += 1)
    println(ok ? "  ✓ " : "  ✗ FAILED: ", msg)
end
banner(s) = (println(); println("─"^72); println("▶ ", s); println("─"^72))

# ---------------------------------------------------------------- loading

banner("1. Load ta_with_topography.nc — discovery chain in action")
ta = read_var(joinpath(EXAMPLES, "ta_with_topography.nc"), "ta")
show(stdout, MIME"text/plain"(), ta)
check(axis_names(ta.grid) == ("time", "lon", "lat", "z_reference"), "axes are the 4 data dims")
check("z_physical" in coord_names(ta.grid), "z_physical promoted (name convention — no CF link in file)")
check("date" in coord_names(ta.grid), "date promoted as 1D aux coord over time")
check(coord(ta, "z_physical").spans == ("lon", "lat", "z_reference"), "z_physical spans 3 axes")

banner("2. Load rasm.nc — curvilinear grid, bare index dims")
rasm = read_var(joinpath(EXAMPLES, "rasm.nc"), "Tair")
show(stdout, MIME"text/plain"(), rasm)
check(coord(rasm, "xc").spans == ("x", "y") && coord(rasm, "yc").spans == ("x", "y"),
    "xc/yc are 2D aux coords over (x, y)")
check(coord(rasm, "x").values == 1:275, "x is a bare index axis (no coordinate variable in file)")

# ------------------------------------------------------- generic lookup

banner("3. Generic coordinate lookup: one contract for 1D and N-D coords")
I = CartesianIndex(1, 10, 5, 2)
zr = coordvalue(ta, "z_reference", I)
zp = coordvalue(ta, "z_physical", I)
println("  at index $(Tuple(I)): lon = $(coordvalue(ta, "lon", I)), lat = $(coordvalue(ta, "lat", I)),")
println("  z_reference = $zr (nominal level) vs z_physical = $(round(zp, digits = 1)) (terrain-warped)")
check(zp > zr, "terrain lifts the physical altitude above the nominal level")

# ta was written with data == z_physical: coordinate_field makes the check a one-liner
err = maximum(abs.(ta.data .- coordinate_field(ta, "z_physical")))
check(err == 0, "ta.data == coordinate_field(ta, \"z_physical\") exactly (file ground truth)")

# --------------------------------------------------- subsetting semantics

banner("4. Subsetting: aux coords follow the axes they span")
w = subset(ta; z_reference = 1:3, lon = 5:14)
check(size(w.data) == (1, 10, 18, 3), "data windowed")
check(size(coord(w, "z_physical").values) == (10, 18, 3),
    "z_physical windowed along its shared axes (lon, z_reference)")
check(coord(w, "z_physical").values ==
      coord(ta, "z_physical").values[5:14, :, 1:3], "…with exactly the right slab")

s = subset(ta; lat = 5)                       # integer index drops the axis
check(axis_names(s.grid) == ("time", "lon", "z_reference"), "integer index dropped lat axis")
check(coord(s, "z_physical").spans == ("lon", "z_reference") &&
      ndims(coord(s, "z_physical").values) == 2, "z_physical lost the lat axis: 3D → 2D")
check(s.data == coordinate_field(s, "z_physical"), "data/coordinate correspondence survives")

p = subset(rasm; x = 100, y = 150)            # drop both spanned axes
xc0 = coord(p, "xc")
check(xc0.spans == () && ndims(xc0.values) == 0,
    "dropping every spanned axis leaves a scalar coordinate (xc = $(round(xc0.values[], digits = 2)))")

# ---------------------------------------------------- value-based slicing

banner("5. Value-based slice: works on dimension coords, guarded on aux coords")
sl = slice(ta; lon = -87.0, z_reference = 3100.0)
println("  slice(ta; lon = -87.0, z_reference = 3100.0) → lon = $(coord(sl, "lon").values[]), z_reference = $(coord(sl, "z_reference").values[])")
check(axis_names(sl.grid) == ("time", "lat"), "sliced axes dropped")
check(coord(sl, "lon").spans == () && coord(sl, "z_reference").spans == (),
    "sliced values survive as scalar coordinates (vs today's attributes[\"slice_lon\"] hack)")
check(coord(sl, "z_physical").spans == ("lat",), "z_physical reduced to the surviving lat axis")
msg = try
    slice(ta; z_physical = 3000.0)
    ""
catch e
    sprint(showerror, e)
end
println("  slice(ta; z_physical = 3000.0) →\n    ERROR: ", msg)
check(occursin("spans", msg) && occursin("interp_to_coord_levels", msg),
    "aux-coord slicing rejected with an actionable message")

# --------------------------------------------------------------- reductions

banner("6. Reductions: aux coords spanning a reduced axis are dropped")
a1 = average(ta; over = "time")
check("z_physical" in coord_names(a1.grid) && "date" ∉ coord_names(a1.grid),
    "averaging over time drops date (spans time) but keeps z_physical (doesn't)")
a2 = average(ta; over = ("lon", "lat"))
check("z_physical" ∉ coord_names(a2.grid),
    "averaging over lon/lat drops z_physical (spanned a reduced axis)")
check(size(a2.data) == (1, 10), "reduced data shape correct")

# ------------------------------------------------------------- permutedims

banner("7. permutedims: aux coordinate arrays are untouched")
tp = permutedims(ta, ("z_reference", "lat", "lon", "time"))
check(size(tp.data) == (10, 18, 36, 1), "data permuted")
check(coord(tp, "z_physical").values === coord(ta, "z_physical").values,
    "z_physical array is the identical object (spans carry its own axis order)")
check(coordvalue(tp, "z_physical", CartesianIndex(2, 5, 10, 1)) ==
      coordvalue(ta, "z_physical", CartesianIndex(1, 10, 5, 2)),
    "coordvalue agrees across the permutation")

# --------------------------------------------------------------------- cat

banner("8. cat along an axis an aux coord spans")
lo = subset(ta; z_reference = 1:5)
hi = subset(ta; z_reference = 6:10)
glued = cat(lo, hi; along = "z_reference")
check(glued.data == ta.data, "data reassembled")
check(coord(glued, "z_physical").values == coord(ta, "z_physical").values,
    "z_physical reassembled by concatenating along ITS OWN matching axis")
# two runs with different topography must refuse to concatenate along time
tweak = Var(
    copy(ta.attributes),
    Grid(ta.grid.dims, map(
        c -> c.name == "z_physical" ? Coord(c.name, c.spans, c.values .+ 1.0, c.attributes) : c,
        ta.grid.coords,
    )),
    ta.data,
)
bad = try
    cat(ta, tweak; along = "time")
    ""
catch e
    sprint(showerror, e)
end
check(occursin("'z_physical' differs", bad),
    "cat refuses two vars with different terrain (aux coord mismatch off the cat axis)")
t2 = cat(ta, ta; along = "time")
check(size(t2.data) == (2, 36, 18, 10) && size(coord(t2, "z_physical").values) == (36, 18, 10),
    "cat along time: z_physical NOT duplicated (doesn't span time), date concatenated")

# ---------------------------------------- the payoff: physical-level interp

banner("9. interp_to_coord_levels: terrain-following → fixed physical altitudes")
targets = [2000.0, 5000.0, 8000.0]
fixed = interp_to_coord_levels(ta, "z_physical", targets; along = "z_reference")
show(stdout, MIME"text/plain"(), fixed)
check(axis_names(fixed.grid) == ("time", "lon", "lat", "z_physical"),
    "z_reference axis replaced by a real z_physical dimension coordinate")
check("z_physical" ∉ coord_names(fixed.grid), "the 3D aux coord was consumed")
# ground truth: data == z_physical, so interpolated values must equal the targets
for (t, target) in enumerate(targets)
    lvl = view(fixed.data, 1, :, :, t)
    ok = view(lvl, .!isnan.(lvl))
    maxerr = isempty(ok) ? 0.0 : maximum(abs.(ok .- target))
    nnan = count(isnan, lvl)
    above = count(>(target), view(coord(ta, "z_physical").values, :, :, 1))
    println("  level $target m: max |error| = $(round(maxerr, sigdigits = 3)), ",
        "NaN columns = $nnan (columns whose terrain starts above target: $above)")
    check(maxerr < 1e-9, "interpolated values equal the target altitude ($target m)")
    check(nnan == above, "NaNs exactly where the column starts above $target m")
end

# ------------------------------------------------- curvilinear operations

banner("10. Curvilinear rasm: nearest-point lookup on 2D coords")
pt = nearest_point(rasm; xc = 118.33, yc = 57.19)
d2 = (coord(rasm, "xc").values .- 118.33) .^ 2 .+ (coord(rasm, "yc").values .- 57.19) .^ 2
Ibest = argmin(d2)
println("  nearest to (lon 118.33, lat 57.19): cell $(Tuple(Ibest)), ",
    "xc = $(round(coord(pt, "xc").values[], digits = 2)), ",
    "yc = $(round(coord(pt, "yc").values[], digits = 2)), ",
    "Tair[time=3] = $(round(pt.data[3], digits = 2))")
check(Tuple(Ibest) == (100, 150),
    "probe lands on cell (x=100, y=150) — real-file geometry, preserved by the generator")
check(coord(pt, "xc").values[] == coord(rasm, "xc").values[Ibest] &&
      coord(pt, "yc").values[] == coord(rasm, "yc").values[Ibest] &&
      isequal(pt.data[3], rasm.data[Ibest, 3]),
    "subset at the nearest cell carries exactly that cell's coords and data")

banner("11. Curvilinear rasm: 2D-coordinate-weighted reduction")
wts = cosd.(coordinate_field(rasm, "yc"))
gm = weighted_average(rasm, wts; over = ("x", "y"))
check(axis_names(gm.grid) == ("time",) && isempty(coord_names(gm.grid)),
    "area-weighted mean over the curvilinear plane (xc/yc consumed and dropped)")
println("  cos(lat)-weighted mean Tair, first 3 months: ",
    join(round.(gm.data[1:3], digits = 2), ", "), " °C")
check(all(-60 .< gm.data .< 30), "values are physically sensible")

wsub = subset(rasm; x = 50:150, y = 100:180)
check(size(coord(wsub, "xc").values) == (101, 81), "curvilinear window: 2D coords windowed with the data")

# ------------------------------------------------------ discovery guard

banner("12. Discovery guard: sibling data variable is NOT promoted")
tmp = joinpath(mktempdir(), "two_vars.nc")
cp(joinpath(EXAMPLES, "ta_with_topography.nc"), tmp; force = true)
NCDataset(tmp, "a") do ds
    defVar(ds, "tb", rand(36, 18, 10, 1), ("lon", "lat", "z_reference", "time"),
        attrib = Dict("units" => "K", "long_name" => "sibling data variable"))
end
ta2 = read_var(tmp, "ta")
check("tb" ∉ coord_names(ta2.grid) && "z_physical" in coord_names(ta2.grid),
    "tb (dims ⊆ ta dims!) rejected; z_physical still promoted — the chain, not the subset rule, decides")

# -------------------------------------------------------------------- done

println()
println("═"^72)
println("  $(npass[]) checks passed, $(nfail[]) failed")
println("═"^72)
