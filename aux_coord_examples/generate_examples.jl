# Regenerates the two example NetCDF files with the exact structure of the
# originals (dimensions, variable shapes and C-order, attributes, CF quirks),
# so the large binaries need not live in git history.
#
#   - ta_with_topography.nc reproduces the ClimaDiagnostics recipe exactly
#     (EquiangularCubedSphere sampled on 36×18×10, Hypsography.LinearAdaption
#     with z_surface = 1000*(cosd(lat)+cosd(lon)+1), LevelsMethod; the data
#     variable ta is deliberately set equal to z_physical) — values match the
#     original writer output, including the missing `coordinates` attribute.
#   - rasm.nc (originally 17 MB of VIC model output from
#     github.com/pydata/xarray-data) gets synthetic values on a synthetic
#     curvilinear grid; structure is identical (bare x/y index dims, 2D xc/yc
#     with `Tair:coordinates = "yc xc"`, noleap calendar, _FillValue holes).
#     The grid is anchored so the cell nearest (118.33°E, 57.19°N) is
#     (x=100, y=150), matching the real file's geometry at the demo's probe.
#
# Run from the repo root:  include("aux_coord_examples/generate_examples.jl")

using NCDatasets

let path = joinpath(@__DIR__, "ta_with_topography.nc")
    nlon, nlat, nz = 36, 18, 10
    lon = collect(range(-180.0, 180.0; length = nlon))
    lat = collect(range(-90.0, 90.0; length = nlat))
    z_ref = collect(500.0:1000.0:9500.0)
    z_top = 10_000.0
    z_surface = [1000.0 * (cosd(la) + cosd(lo) + 1.0) for lo in lon, la in lat]
    z_phys = [
        z_ref[k] + z_surface[i, j] * (1.0 - z_ref[k] / z_top)
        for i in 1:nlon, j in 1:nlat, k in 1:nz
    ]
    NCDataset(path, "c") do ds
        defDim(ds, "time", Inf)
        defDim(ds, "nv", 2)
        defVar(ds, "time", [0.0], ("time",), attrib = Dict(
            "units" => "s", "axis" => "T", "standard_name" => "time",
            "long_name" => "Time", "bounds" => "time_bnds"))
        defVar(ds, "lon", lon, ("lon",), attrib = Dict(
            "units" => "degrees_east", "axis" => "X",
            "standard_name" => "longitude", "long_name" => "Longitude"))
        defVar(ds, "lat", lat, ("lat",), attrib = Dict(
            "units" => "degrees_north", "axis" => "Y",
            "standard_name" => "latitude", "long_name" => "Latitude"))
        defVar(ds, "z_reference", z_ref, ("z_reference",), attrib = Dict(
            "units" => "m", "axis" => "Z"))
        # like ClimaDiagnostics: units only — no standard_name, no `coordinates`
        # link from ta (see OVERVIEW.md §7.7)
        defVar(ds, "z_physical", z_phys, ("lon", "lat", "z_reference"),
            attrib = Dict("units" => "m"))
        defVar(ds, "time_bnds", reshape([0.0, 0.0], 2, 1), ("nv", "time"),
            attrib = Dict("comments" => "time bounds for each time value", "units" => "s"))
        defVar(ds, "date", [0.0], ("time",), attrib = Dict(
            "units" => "seconds since 2010-01-01T00:00:00", "bounds" => "date_bnds"))
        defVar(ds, "date_bnds", reshape([0.0, 0.0], 2, 1), ("nv", "time"),
            attrib = Dict("comments" => "date bounds for each date value",
                "units" => "seconds since 2010-01-01T00:00:00"))
        defVar(ds, "ta", reshape(z_phys, 1, nlon, nlat, nz),
            ("time", "lon", "lat", "z_reference"), attrib = Dict(
                "short_name" => "ta",
                "long_name" => "Example variable on topography-warped grid",
                "units" => "K", "comments" => "",
                "start_date" => "2010-01-01T00:00:00"))
    end
    println("wrote $path ($(round(filesize(path) / 1024, digits = 1)) KiB)")
end

let path = joinpath(@__DIR__, "rasm.nc")
    nx, ny, nt = 275, 205, 36
    # sheared plane anchored at the probe cell: xc/yc genuinely vary along both axes
    xc = [118.33 + 0.4 * (i - 100) - 0.15 * (j - 150) for i in 1:nx, j in 1:ny]
    yc = [57.19 + 0.25 * (j - 150) + 0.05 * (i - 100) for i in 1:nx, j in 1:ny]
    # monthly midpoints Sep 1980 … Aug 1983, noleap calendar, days since 0001-01-01
    mlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    times = Float64[]
    let y = 1980, m = 9, start = 365.0 * (1980 - 1) + sum(mlen[1:8])
        for _ in 1:nt
            push!(times, start + mlen[m] / 2)
            start += mlen[m]
            m = m == 12 ? (y += 1; 1) : m + 1
        end
    end
    ocean = [sin(i / 23) + cos(j / 17) > 1.2 for i in 1:nx, j in 1:ny]
    @assert !ocean[100, 150] "probe cell must stay land"
    Tair = Array{Union{Missing, Float64}}(undef, nx, ny, nt)
    for t in 1:nt, j in 1:ny, i in 1:nx
        Tair[i, j, t] = ocean[i, j] ? missing :
            28.0 * cosd(yc[i, j]) - 18.0 + 12.0 * cospi(2 * (t - 8) / 12) +
            2.0 * sin(i / 9) * cos(j / 7)
    end
    NCDataset(path, "c") do ds
        defVar(ds, "Tair", Tair, ("x", "y", "time");
            fillvalue = 9.96920996838687e36,
            attrib = Dict("units" => "C", "long_name" => "Surface air temperature",
                "type_preferred" => "double", "time_rep" => "instantaneous",
                "coordinates" => "yc xc"))
        defVar(ds, "time", times, ("time",), attrib = Dict(
            "long_name" => "time", "units" => "days since 0001-01-01",
            "calendar" => "noleap"))
        # x and y deliberately get no coordinate variables (bare index dims)
        defVar(ds, "xc", xc, ("x", "y"), attrib = Dict(
            "long_name" => "longitude of grid cell center", "units" => "degrees_east"))
        defVar(ds, "yc", yc, ("x", "y"), attrib = Dict(
            "long_name" => "latitude of grid cell center", "units" => "degrees_north"))
        ds.attrib["title"] = "synthetic stand-in for pydata/xarray-data rasm.nc (structure-identical, synthetic values)"
        ds.attrib["convention"] = "CF-1.4"
    end
    println("wrote $path ($(round(filesize(path) / 1024^2, digits = 1)) MiB)")
end
