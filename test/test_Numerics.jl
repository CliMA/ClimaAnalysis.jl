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
