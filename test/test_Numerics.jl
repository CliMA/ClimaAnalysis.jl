using Test
import ClimaAnalysis
import ClimaAnalysis.Numerics.Interpolation as Intp
import OrderedCollections: OrderedDict

@testset "integration weights for lon and lat" begin
    # Integration weights for lon (not equispaced)
    lon = [-180.0, -45.0, 100.0, 180.0]
    lon_weights = [135.0, 145.0, 80.0, 0.0] .* (π / 180.0)
    @test all(
        isapprox.(
            lon_weights,
            ClimaAnalysis.Numerics._integration_weights_lon_left(lon),
        ),
    )

    # Integration weights for lat (not equispaced)
    lat = [-90.0, 20.0, 45.0, 90.0]
    lat_weights = [110.0, 25.0, 45.0, 0.0] .* (π / 180.0) .* cosd.(lat)
    @test all(
        isapprox.(
            lat_weights,
            ClimaAnalysis.Numerics._integration_weights_lat_left(lat),
        ),
    )

    # Integration weights for generic dim (not equispaced)
    z = [-180.0, -45.0, 100.0, 180.0]
    z_weights = [135.0, 145.0, 80.0, 0.0]
    @test all(
        isapprox.(
            z_weights,
            ClimaAnalysis.Numerics._integration_weights_generic_left(z),
        ),
    )

    # Integration weights for lon (not equispaced)
    lon = collect(range(-180.0, 180.0, 5))
    lon_weights = [90.0 for _ in lon] .* (π / 180.0)
    @test all(
        isapprox.(
            lon_weights,
            ClimaAnalysis.Numerics._integration_weights_lon_equispaced(lon),
        ),
    )

    # Integration weights for lat (not equispaced)
    lat = collect(range(-90.0, 90.0, 5))
    lat_weights = [45.0 for _ in lat] .* (π / 180.0) .* cosd.(lat)
    @test all(
        isapprox.(
            lat_weights,
            ClimaAnalysis.Numerics._integration_weights_lat_equispaced(lat),
        ),
    )

    # Integration weights for generic dim (not equispaced)
    z = collect(range(-180.0, 180.0, 5))
    z_weights = [90.0 for _ in z]
    @test all(
        isapprox.(
            z_weights,
            ClimaAnalysis.Numerics._integration_weights_generic_equispaced(z),
        ),
    )
end

@testset "Integrating on lon and lat" begin
    # Integrating only lon (non equispaced)
    lon = collect(range(-180.0, 179.0, 100))
    # Force array to be non equispaced for testing _integration_weights_lon
    push!(lon, 180.0)
    lon_data = ones(length(lon))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lon(lon_data, lon, dims = 1)[1],
        2.0π,
        atol = 0.01,
    )

    # Integrating only lat (non equispaced)
    lat = collect(range(-90.0, 89.0, 100))
    # Force array to be non equispaced for testing _integration_weights_lat
    push!(lat, 90.0)
    lat_data = ones(length(lat))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lat(lat_data, lat, dims = 1)[1],
        2.0,
        atol = 0.01,
    )

    # Integrating both lon and lat
    data = ones(length(lat), length(lon))
    integrated_lat = ClimaAnalysis.Numerics._integrate_lat(data, lat, dims = 1)
    integrated_latlon =
        ClimaAnalysis.Numerics._integrate_lon(integrated_lat, lon, dims = 1)

    integrated_lon = ClimaAnalysis.Numerics._integrate_lon(data, lon, dims = 2)
    integrated_lonlat =
        ClimaAnalysis.Numerics._integrate_lat(integrated_lon, lat, dims = 1)

    # Order of integration should not matter
    @test isapprox(integrated_latlon[1], integrated_lonlat[1])
    @test isapprox(integrated_latlon[1], 4π, atol = 0.01)

    # Error checking for length of lon and lat and values in lon and lat
    @test_throws "Cannot integrate when latitude is a single point" ClimaAnalysis.Numerics._integrate_lat(
        lat_data,
        [0.0],
        dims = 1,
    )

    @test_throws "Cannot integrate when latitude is a single point" ClimaAnalysis.Numerics._integrate_lat(
        lon_data,
        [0.0],
        dims = 1,
    )

    # Integrating only lon (equispaced)
    lon = collect(range(-179.5, 179.5, 360))
    lon_data = ones(length(lon))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lon(lon_data, lon, dims = 1)[1],
        2.0π,
        atol = 0.01,
    )

    # Integrating only lat (equispaced)
    lat = collect(range(-89.5, 89.5, 180))
    lat_data = ones(length(lat))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lat(lat_data, lat, dims = 1)[1],
        2.0,
        atol = 0.01,
    )

    # Integrating with NaNs
    lon = collect(range(-179.5, 179.5, 360))
    lon_data = ones(length(lon))
    lon_data[:] .= NaN
    @test ClimaAnalysis.Numerics._integrate_lon(lon_data, lon, dims = 1)[1] ==
          0.0

    lon = collect(range(-179.5, 179.5, 360))
    lon_data = ones(length(lon))
    lon_data[1:180] .= NaN
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lon(lon_data, lon, dims = 1)[1],
        1.0π,
        atol = 0.01,
    )
end

@testset "Integrating on generic dimension" begin
    # Integrating (equispaced)
    z = [0.0, 1.0, 2.0]
    z_data = ones(length(z))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_dim(z_data, z, dims = 1)[1],
        3.0,
    )

    # Integrating (not equispaced)
    z = [0.0, 1.0, 2.0, 2.5]
    z_data = ones(length(z))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_dim(z_data, z, dims = 1)[1],
        2.5,
    )
end

@testset "Regridder construction" begin
    data = reshape(1.0:9.0, (3, 3))
    grid = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    extp = (Intp.Throw(), Intp.Flat())
    stag = (Intp.Node(), Intp.Center())

    regridder = Intp.Regridder(data, grid, extp, stag)
    @test regridder.src_grid == grid
    @test regridder.src_data == data
    @test regridder.extrapolation == extp
    @test regridder.staggering == stag
    @test regridder.method == Intp.LinearInterpolation()

    # Extrapolation and staggering can be any iterable
    regridder = Intp.Regridder(
        data,
        grid,
        [Intp.Throw(), Intp.Flat()],
        [Intp.Node(), Intp.Center()],
    )
    @test regridder.extrapolation == extp
    @test regridder.staggering == stag

    # Number of dimensions of everything must match
    @test_throws "Number of extrapolation conditions" Intp.Regridder(
        data,
        grid,
        (Intp.Throw(),),
        stag,
    )
    @test_throws "Number of staggering" Intp.Regridder(
        data,
        grid,
        extp,
        (Intp.Node(),),
    )
    @test_throws "Number of dimensions in src_grid" Intp.Regridder(
        data,
        (grid[1],),
        extp,
        stag,
    )
    @test_throws "do not match size(src_data)" Intp.Regridder(
        data,
        ([1.0, 2.0], [4.0, 5.0, 6.0]),
        extp,
        stag,
    )

    # Coordinates must be strictly increasing with at least two coordinates
    @test_throws "less than two coordinates" Intp.Regridder(
        ones(3, 1),
        ([1.0, 2.0, 3.0], [4.0]),
        extp,
        stag,
    )
    @test_throws "not strictly increasing" Intp.Regridder(
        data,
        ([1.0, 1.0, 2.0], [4.0, 5.0, 6.0]),
        extp,
        stag,
    )
    @test_throws "not strictly increasing" Intp.Regridder(
        data,
        ([3.0, 2.0, 1.0], [4.0, 5.0, 6.0]),
        extp,
        stag,
    )

    # NaN and missing values in coordinates are not supported
    @test_throws "Missing or NaN values" Intp.Regridder(
        data,
        ([1.0, NaN, 3.0], [4.0, 5.0, 6.0]),
        extp,
        stag,
    )
    @test_throws "Missing or NaN values" Intp.Regridder(
        data,
        ([1.0, 2.0, 3.0], Union{Missing, Float64}[4.0, missing, 6.0]),
        extp,
        stag,
    )

    # Period must be positive and not smaller than the span of the coordinates
    @test_throws "is not positive" Intp.Regridder(
        data,
        grid,
        (Intp.Periodic(-360.0), Intp.Flat()),
        stag,
    )
    @test_throws "larger than the period" Intp.Regridder(
        data,
        grid,
        (Intp.Periodic(1.5), Intp.Flat()),
        stag,
    )

    # Longitudes whose span slightly exceeds 360 degrees should still get a
    # periodic boundary condition
    noisy_lon = [0.0, 90.0, 180.0, 270.0, 360.0 + 1.0e-6]
    noisy = Intp.Regridder(
        [1.0, 2.0, 3.0, 4.0, 1.0],
        (noisy_lon,),
        (Intp.Periodic(360.0),),
        (Intp.Node(),),
    )
    # mod(405, 360) is 45
    @test Intp.interpolate(noisy, (405.0,)) ≈ 1.5
end

@testset "Interpolating a point or a vector of points" begin
    # 1D
    coords = ([1.0, 2.0, 3.0],)
    data = [3.0, 1.0, 0.0]

    throw_node = Intp.Regridder(data, coords, (Intp.Throw(),), (Intp.Node(),))
    @test Intp.interpolate(throw_node, (1.0,)) == 3.0
    @test Intp.interpolate(throw_node, (3.0,)) == 0.0
    @test Intp.interpolate(throw_node, (1.5,)) == 2.0
    @test Intp.interpolate(throw_node, [(1.0,), (1.5,), (2.5,)]) ==
          [3.0, 2.0, 0.5]
    @test_throws DomainError Intp.interpolate(throw_node, (0.0,))
    @test_throws DomainError Intp.interpolate(throw_node, (4.0,))

    # Values outside of the domain evaluate to the value at the boundary
    flat_node = Intp.Regridder(data, coords, (Intp.Flat(),), (Intp.Node(),))
    @test Intp.interpolate(flat_node, (0.0,)) == 3.0
    @test Intp.interpolate(flat_node, (4.0,)) == 0.0
    @test Intp.interpolate(flat_node, (1.5,)) == 2.0
    @test isnan(Intp.interpolate(flat_node, (NaN,)))

    # With Center staggering, the domain extends half a cell beyond the first
    # and last coordinates and values in the margins are linearly extended
    throw_center =
        Intp.Regridder(data, coords, (Intp.Throw(),), (Intp.Center(),))
    @test Intp.interpolate(throw_center, (0.75,)) == 3.5
    @test Intp.interpolate(throw_center, (3.25,)) == -0.25
    @test_throws DomainError Intp.interpolate(throw_center, (0.25,))
    @test_throws DomainError Intp.interpolate(throw_center, (3.75,))

    # 2D
    coords = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    data = reshape(1.0:9.0, (3, 3))
    throw2d = Intp.Regridder(
        data,
        coords,
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(throw2d, (1.0, 4.0)) == 1.0
    @test Intp.interpolate(throw2d, (3.0, 6.0)) == 9.0
    @test Intp.interpolate(throw2d, (2.0, 5.0)) == 5.0
    @test Intp.interpolate(throw2d, (1.5, 4.5)) == 3.0
    @test Intp.interpolate(throw2d, (1.5, 5.5)) == 6.0

    # Extrapolation conditions apply per dimension
    flat_throw = Intp.Regridder(
        data,
        coords,
        (Intp.Flat(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test_throws DomainError Intp.interpolate(flat_throw, (0.0, 8.0))
    @test Intp.interpolate(flat_throw, (0.0, 5.5)) == 5.5

    flat_flat = Intp.Regridder(
        data,
        coords,
        (Intp.Flat(), Intp.Flat()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(flat_flat, (0.0, 8.0)) == 7.0

    # 3D
    coords = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0])
    data = reshape(1.0:27.0, (3, 3, 3))
    throw3d = Intp.Regridder(
        data,
        coords,
        (Intp.Throw(), Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(throw3d, (1.0, 5.0, 7.0)) == 4.0
    @test Intp.interpolate(throw3d, (1.5, 5.2, 7.5)) ≈ 9.6

    # Non equispaced coordinates
    coords = ([1.0, 3.0, 7.0], [4.0, 5.0, 7.0])
    data = reshape(1.0:9.0, (3, 3))
    noneq = Intp.Regridder(
        data,
        coords,
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(noneq, (2.0, 4.5)) == 3.0
    @test Intp.interpolate(noneq, (5.0, 6.0)) == 7.0
end

@testset "Periodic boundary condition" begin
    # Note that the data at the endpoints do not agree
    lon = [0.0, 90.0, 180.0, 270.0, 360.0]
    data = [1.0, 2.0, 3.0, 4.0, 10.0]
    node = Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Node(),))
    @test Intp.interpolate(node, (0.0,)) == 1.0
    @test Intp.interpolate(node, (360.0,)) == 10.0
    @test Intp.interpolate(node, (720.0,)) == 1.0
    @test Intp.interpolate(node, (45.0,)) == 1.5
    @test Intp.interpolate(node, (405.0,)) == 1.5
    @test Intp.interpolate(node, (-45.0,)) == 7.0
    @test Intp.interpolate(node, (-90.0,)) == 4.0

    # Center staggering without a duplicated endpoint; points between the
    # last and first coordinates interpolate across the wrap
    lon = [45.0, 135.0, 225.0, 315.0]
    data = [1.0, 2.0, 3.0, 4.0]
    center =
        Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Center(),))
    @test Intp.interpolate(center, (45.0,)) == 1.0
    # The points 0.0 and 360.0 represents the same location
    @test Intp.interpolate(center, (0.0,)) == 2.5
    @test Intp.interpolate(center, (360.0,)) == 2.5
    @test Intp.interpolate(center, (337.5,)) == 3.25
    @test Intp.interpolate(center, (405.0,)) == 1.0
end

@testset "Regridding" begin
    src_grid = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    data = reshape(1.0:9.0, (3, 3))
    regridder = Intp.Regridder(
        data,
        src_grid,
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )

    # Regridding onto the same grid reproduces the data
    @test Intp.regrid(regridder, src_grid) == data

    # Regridding agrees with interpolating each point
    dest_grid = ([1.0, 1.5, 2.0, 2.5, 3.0], [4.0, 4.5, 5.0])
    expected = [
        Intp.interpolate(regridder, (x, y)) for
        x in dest_grid[1], y in dest_grid[2]
    ]
    @test Intp.regrid(regridder, dest_grid) == expected

    # The same regridder can be reused with a different destination grid
    another_dest_grid = ([1.5, 2.5], [4.5, 5.5])
    another_expected = [
        Intp.interpolate(regridder, (x, y)) for
        x in another_dest_grid[1], y in another_dest_grid[2]
    ]
    @test Intp.regrid(regridder, another_dest_grid) == another_expected

    # In-place regridding
    dest_data = zeros(5, 3)
    Intp.regrid!(dest_data, regridder, dest_grid)
    @test dest_data == expected

    # OrderedDict adapters
    dest_dict = OrderedDict("lon" => dest_grid[1], "lat" => dest_grid[2])
    @test Intp.regrid(regridder, dest_dict) == expected
    fill!(dest_data, 0.0)
    Intp.regrid!(dest_data, regridder, dest_dict)
    @test dest_data == expected

    # Check size and dimensions
    @test_throws DimensionMismatch Intp.regrid!(
        zeros(2, 2),
        regridder,
        dest_grid,
    )
    @test_throws ArgumentError Intp.regrid!(
        zeros(5),
        regridder,
        (dest_grid[1],),
    )

    # Test with points in no particular order
    @test Intp.regrid(regridder, ([3.0, 1.0, 2.0], [4.5, 4.5])) ==
          [4.5 4.5; 2.5 2.5; 3.5 3.5]

    # Empty destination grids
    @test size(Intp.regrid(regridder, (Float64[], [4.0]))) == (0, 1)
    @test size(Intp.regrid(regridder, (Float64[], Float64[]))) == (0, 0)

    # A point has to have as many coordinates as the regridder has dimensions
    @test_throws MethodError Intp.interpolate(regridder, (1.0,))
    @test_throws MethodError Intp.interpolate(regridder, (1.0, 4.0, 7.0))
end

@testset "NaN-aware regridding" begin
    @test Intp.NaNLinearInterpolation().threshold == 0.5
    @test Intp.NaNLinearInterpolation(0).threshold === 0.0
    @test Intp.NaNLinearInterpolation(1.0).threshold === 1.0
    @test Intp.NaNLinearInterpolation(1 // 2).threshold === 0.5
    # Threshold must be between 0 and 1
    @test_throws "must be between 0 and 1" Intp.NaNLinearInterpolation(-0.5)
    @test_throws "must be between 0 and 1" Intp.NaNLinearInterpolation(1.5)
    @test_throws "must be between 0 and 1" Intp.NaNLinearInterpolation(NaN)

    # extp and stag apply to every dimension unless given as tuples
    function make_regridder(
        data,
        coords;
        extp = Intp.Throw(),
        stag = Intp.Node(),
        threshold = 0.5,
    )
        n = length(coords)
        return Intp.Regridder(
            data,
            coords,
            extp isa Tuple ? extp : ntuple(_ -> extp, n),
            stag isa Tuple ? stag : ntuple(_ -> stag, n),
            method = Intp.NaNLinearInterpolation(threshold),
        )
    end

    # 1D
    coords = ([1.0, 2.0, 3.0],)
    data = [1.0, NaN, 3.0]
    nan_regridder(threshold) = make_regridder(data, coords; threshold)

    # NaN with zero weight does not affect the result for any threshold
    for threshold in (0.0, 0.5, 1.0)
        @test Intp.interpolate(nan_regridder(threshold), (1.0,)) == 1.0
    end

    # NaN weight is 0.5: below, at, and above the threshold
    @test isnan(Intp.interpolate(nan_regridder(0.25), (1.5,)))
    # Exceeding the threshold means strictly greater
    @test Intp.interpolate(nan_regridder(0.5), (1.5,)) == 1.0
    @test Intp.interpolate(nan_regridder(0.75), (1.5,)) == 1.0
    # Weights are renormalized by the sum of the non-NaN weights
    @test Intp.interpolate(nan_regridder(0.75), (1.25,)) == 1.0
    @test Intp.interpolate(nan_regridder(0.75), (2.75,)) == 3.0

    # NaN weight is 1, so the result is NaN even when the threshold is 1
    @test isnan(Intp.interpolate(nan_regridder(1.0), (2.0,)))

    # Regridding onto the source grid reproduces the data, including NaNs
    @test Intp.regrid(nan_regridder(0.5), coords) ≈ data nans = true

    # Interpolating a vector of points
    @test isequal(
        Intp.interpolate(nan_regridder(0.5), [(1.0,), (1.5,), (2.0,)]),
        [1.0, 1.0, NaN],
    )

    # Renormalization returns the only non-NaN value when the threshold is not
    # exceeded
    sparse = make_regridder([NaN, 2.0, NaN], coords)
    @test Intp.interpolate(sparse, (1.5,)) == 2.0
    @test Intp.interpolate(sparse, (2.5,)) == 2.0
    @test isnan(Intp.interpolate(sparse, (1.25,)))

    # regrid! into a Float32 destination array converts values and NaNs
    dest32 = zeros(Float32, 2)
    Intp.regrid!(dest32, nan_regridder(0.25), ([1.0, 1.5],))
    @test dest32[1] === 1.0f0
    @test isnan(dest32[2])

    # All NaN source data produces NaN for any threshold
    all_nan = make_regridder([NaN, NaN, NaN], coords, threshold = 1.0)
    @test all(isnan, Intp.regrid(all_nan, ([1.0, 1.5, 2.0, 3.0],)))

    # Without NaNs, NaNLinearInterpolation agrees with LinearInterpolation
    src_grid = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    data2d = collect(reshape(1.0:9.0, (3, 3)))
    dest_grid = ([1.0, 1.5, 2.5, 3.0], [4.0, 4.5, 5.5, 6.0])
    linear = Intp.Regridder(
        data2d,
        src_grid,
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    nan_linear = make_regridder(data2d, src_grid)
    @test Intp.regrid(nan_linear, dest_grid) == Intp.regrid(linear, dest_grid)

    # TODO: Review this
    # 2D with a single NaN: the query point (1.5, 4.5) weights all four
    # corners by 0.25
    data_nan = collect(reshape(1.0:9.0, (3, 3)))
    data_nan[1, 1] = NaN
    nan2d(threshold) = make_regridder(data_nan, src_grid; threshold)
    # (2 + 4 + 5) * 0.25 / 0.75
    @test Intp.interpolate(nan2d(0.25), (1.5, 4.5)) ≈ 11.0 / 3.0
    @test isnan(Intp.interpolate(nan2d(0.2), (1.5, 4.5)))
    # Points on the grid away from the NaN are unaffected
    @test Intp.interpolate(nan2d(0.0), (2.0, 5.0)) == 5.0
    @test isnan(Intp.interpolate(nan2d(0.0), (1.0, 4.0)))

    # Asymmetric weights: the query point (1.5, 4.25) weights the NaN corner
    # by 0.375, so the result is (2 * 0.375 + 4 * 0.125 + 5 * 0.125) / 0.625
    @test Intp.interpolate(nan2d(0.375), (1.5, 4.25)) == 3.0
    @test isnan(Intp.interpolate(nan2d(0.25), (1.5, 4.25)))

    # 2D with two NaNs: NaN weight at (1.5, 4.5) is 0.5
    data_nan2 = collect(reshape(1.0:9.0, (3, 3)))
    data_nan2[1, 1] = NaN
    data_nan2[2, 2] = NaN
    nan2d_two(threshold) = make_regridder(data_nan2, src_grid; threshold)
    # (2 + 4) * 0.25 / 0.5
    @test Intp.interpolate(nan2d_two(0.5), (1.5, 4.5)) == 3.0
    @test isnan(Intp.interpolate(nan2d_two(0.25), (1.5, 4.5)))

    # 2D with three NaN corners: only one valid corner remains
    data_nan3 = collect(reshape(1.0:9.0, (3, 3)))
    data_nan3[1, 1] = NaN
    data_nan3[2, 1] = NaN
    data_nan3[1, 2] = NaN
    nan2d_three(threshold) = make_regridder(data_nan3, src_grid; threshold)
    @test Intp.interpolate(nan2d_three(0.75), (1.5, 4.5)) == 5.0
    @test isnan(Intp.interpolate(nan2d_three(0.5), (1.5, 4.5)))

    # 2D with all four NaN corners: NaN even when the threshold is 1
    data_nan4 = copy(data_nan3)
    data_nan4[2, 2] = NaN
    nan2d_four = make_regridder(data_nan4, src_grid, threshold = 1.0)
    @test isnan(Intp.interpolate(nan2d_four, (1.5, 4.5)))

    # 3D with one NaN corner: the center of the cell weights all eight corners
    # by 0.125
    data3d = collect(reshape(1.0:8.0, (2, 2, 2)))
    data3d[1, 1, 1] = NaN
    nan3d(threshold) =
        make_regridder(data3d, ([1.0, 2.0], [1.0, 2.0], [1.0, 2.0]); threshold)
    # (2 + 3 + ... + 8) * 0.125 / 0.875
    @test Intp.interpolate(nan3d(0.125), (1.5, 1.5, 1.5)) == 5.0
    @test isnan(Intp.interpolate(nan3d(0.1), (1.5, 1.5, 1.5)))

    # regrid, regrid!, and interpolate agree
    expected = [
        Intp.interpolate(nan2d(0.25), (x, y)) for
        x in dest_grid[1], y in dest_grid[2]
    ]
    regridded = Intp.regrid(nan2d(0.25), dest_grid)
    @test regridded ≈ expected nans = true
    dest_data = zeros(4, 4)
    Intp.regrid!(dest_data, nan2d(0.25), dest_grid)
    @test dest_data ≈ expected nans = true

    # Extrapolation conditions apply per dimension with NaN handling
    lon_grid = [0.0, 90.0, 180.0, 270.0]
    data_mixed = collect(reshape(1.0:12.0, (4, 3)))
    data_mixed[1, 1] = NaN
    mixed_extp(threshold) = make_regridder(
        data_mixed,
        (lon_grid, [4.0, 5.0, 6.0]);
        extp = (Intp.Periodic(360.0), Intp.Flat()),
        threshold,
    )
    # The wrap gap interpolates between the last (4.0) and first (NaN) points
    @test Intp.interpolate(mixed_extp(0.5), (315.0, 4.0)) == 4.0
    @test isnan(Intp.interpolate(mixed_extp(0.25), (315.0, 4.0)))
    # Outside the flat dimension, away from the NaN
    @test Intp.interpolate(mixed_extp(0.5), (315.0, 8.0)) == 10.5
    # Outside the flat dimension, next to the NaN
    @test Intp.interpolate(mixed_extp(0.5), (45.0, 3.0)) == 2.0

    # NaN-aware interpolation with periodic and flat boundary conditions
    lon = [0.0, 90.0, 180.0, 270.0]
    data_periodic = [NaN, 2.0, 3.0, 4.0]
    periodic =
        make_regridder(data_periodic, (lon,), extp = Intp.Periodic(360.0))
    # Interpolating across the wrap gap (270, 360) between 4.0 and NaN
    @test Intp.interpolate(periodic, (315.0,)) == 4.0
    @test isnan(Intp.interpolate(periodic, (350.0,)))

    # NaN at the duplicated periodic endpoint
    lon_dup = [0.0, 90.0, 180.0, 270.0, 360.0]
    periodic_dup = make_regridder(
        [1.0, 2.0, 3.0, 4.0, NaN],
        (lon_dup,),
        extp = Intp.Periodic(360.0),
    )
    @test Intp.interpolate(periodic_dup, (0.0,)) == 1.0
    @test isnan(Intp.interpolate(periodic_dup, (360.0,)))
    # Wrapped queries interpolate against the NaN endpoint
    @test Intp.interpolate(periodic_dup, (-45.0,)) == 4.0
    # Queries wrapped by whole periods
    @test Intp.interpolate(periodic_dup, (720.0,)) == 1.0
    @test Intp.interpolate(periodic_dup, (-405.0,)) == 4.0

    # Interpolating across the wrap gap with Center staggering
    lon_center = [45.0, 135.0, 225.0, 315.0]
    periodic_center(threshold) = make_regridder(
        [1.0, 2.0, 3.0, NaN],
        (lon_center,);
        extp = Intp.Periodic(360.0),
        stag = Intp.Center(),
        threshold,
    )
    # The point 0.0 is midway between the last (NaN) and first coordinates
    @test Intp.interpolate(periodic_center(0.5), (0.0,)) == 1.0
    @test isnan(Intp.interpolate(periodic_center(0.25), (0.0,)))

    flat = make_regridder([NaN, 2.0, 3.0], coords, extp = Intp.Flat())
    @test isnan(Intp.interpolate(flat, (0.0,)))
    @test Intp.interpolate(flat, (4.0,)) == 3.0
    # A NaN query coordinate produces NaN
    @test isnan(Intp.interpolate(flat, (NaN,)))

    # Flat with Center staggering behaves the same as with Node
    flat_center = make_regridder(
        [NaN, 2.0, 3.0],
        coords,
        extp = Intp.Flat(),
        stag = Intp.Center(),
    )
    @test isnan(Intp.interpolate(flat_center, (0.5,)))
    @test Intp.interpolate(flat_center, (3.75,)) == 3.0

    # With Center staggering and Throw, the half-cell margins extrapolate with
    # weights outside [0, 1]
    throw_center(data, threshold) =
        make_regridder(data, coords; stag = Intp.Center(), threshold)
    # A NaN at the near boundary has weight greater than one, so the result is
    # NaN for any threshold
    @test isnan(Intp.interpolate(throw_center([NaN, 2.0, 3.0], 1.0), (0.75,)))
    @test isnan(Intp.interpolate(throw_center([1.0, 2.0, NaN], 1.0), (3.25,)))
    # A NaN at the far point has negative weight, so it is excluded and the
    # result degrades to the boundary value
    @test Intp.interpolate(throw_center([1.0, NaN, 3.0], 0.0), (0.75,)) == 1.0
    @test Intp.interpolate(throw_center([1.0, NaN, 3.0], 0.0), (3.25,)) == 3.0
    # A missing value propagates even from an extrapolated corner
    @test ismissing(
        Intp.interpolate(
            throw_center(Union{Missing, Float64}[1.0, missing, 3.0], 0.0),
            (0.75,),
        ),
    )
    # Out-of-domain queries still throw regardless of NaNs
    @test_throws DomainError Intp.interpolate(
        throw_center([NaN, 2.0, 3.0], 1.0),
        (0.25,),
    )

    # Missing values propagate through the weighted sum, like with
    # LinearInterpolation
    grid4 = ([1.0, 2.0, 3.0, 4.0],)
    data_missing = Union{Missing, Float64}[1.0, 2.0, missing, 4.0]
    missing_regridder(threshold) =
        make_regridder(data_missing, grid4; threshold)
    # Stencils away from the missing value are unaffected
    @test Intp.interpolate(missing_regridder(0.5), (1.5,)) == 1.5
    # Missing propagates for any threshold and any weight, including zero
    @test ismissing(Intp.interpolate(missing_regridder(1.0), (2.5,)))
    @test ismissing(Intp.interpolate(missing_regridder(1.0), (3.0,)))
    @test ismissing(Intp.interpolate(missing_regridder(1.0), (4.0,)))
    # The eltype of the destination data includes Missing
    out = Intp.regrid(missing_regridder(0.5), ([1.5, 2.5],))
    @test eltype(out) == Union{Missing, Float64}
    @test out[1] == 1.5
    @test ismissing(out[2])

    # NaN handling applies alongside missing propagation: the NaN is excluded,
    # the missing propagates
    nan_missing =
        make_regridder(Union{Missing, Float64}[1.0, NaN, missing, 4.0], grid4)
    @test Intp.interpolate(nan_missing, (1.5,)) == 1.0
    @test ismissing(Intp.interpolate(nan_missing, (2.5,)))
    @test ismissing(Intp.interpolate(nan_missing, (3.5,)))

    # All missing source data produces missing
    all_missing = make_regridder(
        Vector{Union{Missing, Float64}}(missing, 3),
        coords,
        threshold = 1.0,
    )
    @test all(ismissing, Intp.regrid(all_missing, ([1.0, 1.5, 2.0],)))

    # Integer source data with missing keeps Missing in the eltype
    int_missing = make_regridder(
        Union{Missing, Int}[1, missing, 3],
        coords,
        threshold = 0.75,
    )
    @test ismissing(Intp.interpolate(int_missing, (1.25,)))
    @test eltype(Intp.regrid(int_missing, ([1.25],))) == Union{Missing, Float64}

    # Float32 source data keeps its element type
    nan32 = make_regridder(Float32[1.0, NaN, 3.0], ([1.0f0, 2.0f0, 3.0f0],))
    @test Intp.interpolate(nan32, (1.5f0,)) === 1.0f0
    @test eltype(Intp.regrid(nan32, (Float32[1.5],))) == Float32

    # Integer source data cannot contain NaN and produces floating point data,
    # matching LinearInterpolation
    int_nan = make_regridder(reshape(1:9, (3, 3)), src_grid)
    @test Intp.interpolate(int_nan, (1.5, 4.5)) === 3.0
    @test eltype(Intp.regrid(int_nan, src_grid)) == Float64
    @test Intp.regrid(int_nan, src_grid) == reshape(1:9, (3, 3))
end

@testset "Interpolating and regridding agree" begin
    # Interpolation happens in the promotion of the float types of the source
    # data and source grid, so regrid! gives the same values as regrid and
    # interpolate assuming that the destination array is the same type
    src_data = Float32[1.0, 1.0000001, 2.0]
    regridder = Intp.Regridder(
        src_data,
        ([1.0, 2.0, 3.0],),
        (Intp.Flat(),),
        (Intp.Node(),),
    )

    intp_val = Intp.interpolate(regridder, (1.5,))
    regridded_val = Intp.regrid(regridder, ([1.5],))[1]
    @test intp_val isa Float64
    @test regridded_val isa Float64
    @test intp_val == regridded_val

    dest64 = zeros(1)
    Intp.regrid!(dest64, regridder, ([1.5],))
    @test dest64[1] == intp_val

    dest32 = zeros(Float32, 1)
    Intp.regrid!(dest32, regridder, ([1.5],))
    @test dest32[1] == Float32(intp_val)
end

@testset "Regridding with different types" begin
    # Coordinates and points with different float types
    data = [1.0 3.0; 2.0 4.0]
    regridder = Intp.Regridder(
        data,
        ([1.0f0, 2.0f0], Float16[3.0, 4.0]),
        (Intp.Flat(), Intp.Flat()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(regridder, (1.5f0, 3.5)) == 2.5
    @test Intp.interpolate(regridder, (1.5, 4.5f0)) == 3.5

    # Mixed float types for coordinates and source data lead to the right float
    # type for destination data
    regridder_mixed = Intp.Regridder(
        Float32[1.0, 2.0, 4.0],
        ([1.0f0, 2.0f0, 3.0f0],),
        (Intp.Flat(),),
        (Intp.Node(),),
    )
    @test Intp.interpolate(regridder_mixed, (2.5,)) === 3.0f0
    @test eltype(Intp.regrid(regridder_mixed, ([2.5],))) == Float32

    # The eltype of the destination data is the promotion of the eltypes of the
    # source grid and the source data
    regridder_promote = Intp.Regridder(
        Float32[1.0, 2.0, 4.0],
        ([1.0, 2.0, 3.0],),
        (Intp.Flat(),),
        (Intp.Node(),),
    )
    @test Intp.interpolate(regridder_promote, (2.5,)) === 3.0
    @test eltype(Intp.regrid(regridder_promote, ([2.5],))) == Float64

    # Integer source data produces floating point destination data
    regridder_int = Intp.Regridder(
        reshape(1:9, (3, 3)),
        ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(regridder_int, (2.0, 5.0)) === 5.0
    out = Intp.regrid(regridder_int, ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]))
    @test eltype(out) == Float64
    @test out == reshape(1:9, (3, 3))

    # Element type of destination data is Union{Missing, T} if element type of
    # source data is Union{Missing, T} too
    lon_m = convert(Vector{Union{Missing, Float64}}, [1.0, 2.0, 3.0])
    data_m = convert(
        Matrix{Union{Missing, Float64}},
        collect(reshape(1.0:9.0, (3, 3))),
    )
    regridder_m = Intp.Regridder(
        data_m,
        (lon_m, [4.0, 5.0, 6.0]),
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test Intp.interpolate(regridder_m, (1.5, 4.5)) === 3.0
    out = Intp.regrid(regridder_m, (lon_m, [4.0, 4.5, 5.0]))
    @test eltype(out) == Union{Missing, Float64}
    @test out == [1.0 2.5 4.0; 2.0 3.5 5.0; 3.0 4.5 6.0]

    # Actual missing values propagate through interpolation and regridding
    data_missing =
        convert(Vector{Union{Missing, Float64}}, [10.0, 20.0, 30.0, 40.0])
    data_missing[3] = missing
    regridder_missing = Intp.Regridder(
        data_missing,
        ([1.0, 2.0, 3.0, 4.0],),
        (Intp.Throw(),),
        (Intp.Node(),),
    )
    @test Intp.interpolate(regridder_missing, (1.5,)) == 15.0
    @test ismissing(Intp.interpolate(regridder_missing, (2.5,)))
    out = Intp.regrid(regridder_missing, ([1.5, 2.5, 3.5],))
    @test eltype(out) == Union{Missing, Float64}
    @test out[1] == 15.0
    @test ismissing(out[2])
    @test ismissing(out[3])
end
