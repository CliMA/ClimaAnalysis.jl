using Test
import ClimaAnalysis
import ClimaAnalysis.Interpolation as Intp
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
    @test regridder.stencil_cache isa
          NTuple{2, Vector{Tuple{Int, Int, Float64}}}

    # Extrapolation and staggering can be any iterable
    regridder = Intp.Regridder(
        data,
        grid,
        [Intp.Throw(), Intp.Flat()],
        [Intp.Node(), Intp.Center()],
    )
    @test regridder.extrapolation == extp
    @test regridder.staggering == stag

    # Weight type is promoted from the eltypes of the coordinates
    regridder = Intp.Regridder(
        [1.0 3.0; 2.0 4.0],
        ([1.0f0, 2.0f0], Float16[3.0, 4.0]),
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
    )
    @test regridder.stencil_cache isa
          NTuple{2, Vector{Tuple{Int, Int, Float32}}}

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

    # Spans within floating-point tolerance of the period are accepted, since
    # coordinates spanning the full period often carry accumulated round-off
    noisy_lon = [0.0, 90.0, 180.0, 270.0, 360.0 + 1.0e-6]
    noisy = Intp.Regridder(
        [1.0, 2.0, 3.0, 4.0, 1.0],
        (noisy_lon,),
        (Intp.Periodic(360.0),),
        (Intp.Node(),),
    )
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

    # A NaN query coordinate produces NaN instead of clamping to a boundary
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
    # Node staggering with a duplicated endpoint (0 and 360 are the same
    # physical point)
    lon = [0.0, 90.0, 180.0, 270.0, 360.0]
    data = [1.0, 2.0, 3.0, 4.0, 1.0]
    node = Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Node(),))
    @test Intp.interpolate(node, (45.0,)) == 1.5
    @test Intp.interpolate(node, (405.0,)) == 1.5
    @test Intp.interpolate(node, (-45.0,)) == 2.5
    @test Intp.interpolate(node, (-90.0,)) == 4.0
    @test Intp.interpolate(node, (720.0,)) == 1.0

    # Center staggering without a duplicated endpoint; points between the
    # last and first coordinates interpolate across the wrap
    lon = [45.0, 135.0, 225.0, 315.0]
    data = [1.0, 2.0, 3.0, 4.0]
    center =
        Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Center(),))
    @test Intp.interpolate(center, (45.0,)) == 1.0
    @test Intp.interpolate(center, (0.0,)) == 2.5
    @test Intp.interpolate(center, (360.0,)) == 2.5
    @test Intp.interpolate(center, (337.5,)) == 3.25
    @test Intp.interpolate(center, (405.0,)) == 1.0
end

@testset "Agreement with Interpolations.jl" begin
    # The expected values below were produced with Interpolations.jl, which is
    # what ClimaAnalysis used before this module existed. The grids are the ones
    # that Var._find_extp_and_staggering produces for global longitudes and
    # latitudes, since those are the cases where the two implementations differ
    # the most in how they get to the answer (Interpolations.jl needed an extra
    # longitude point and a copy of the data to handle cell centers)
    lon = collect(0.5:1.0:359.5)
    data = collect(1.0:360.0)
    center =
        Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Center(),))
    for (x, expected) in [
        (0.0, 180.5),
        (0.5, 1.0),
        (180.0, 180.5),
        (359.5, 360.0),
        (360.0, 180.5),
        (360.5, 1.0),
        (-0.5, 360.0),
        (-10.0, 350.5),
    ]
        @test Intp.interpolate(center, (x,)) == expected
    end

    # Longitudes that are equispaced nodes spanning all 360 degrees are also
    # treated as cell centers
    lon = collect(0.0:1.0:359.0)
    node_as_center =
        Intp.Regridder(data, (lon,), (Intp.Periodic(360.0),), (Intp.Center(),))
    for (x, expected) in [
        (0.0, 1.0),
        (359.0, 360.0),
        (359.5, 180.5),
        (360.0, 1.0),
        (361.0, 2.0),
        (-1.0, 360.0),
        (-0.5, 180.5),
    ]
        @test Intp.interpolate(node_as_center, (x,)) == expected
    end

    lat = collect(-89.5:1.0:89.5)
    data = collect(1.0:180.0)
    flat = Intp.Regridder(data, (lat,), (Intp.Flat(),), (Intp.Center(),))
    for (x, expected) in [
        (-90.0, 1.0),
        (-89.5, 1.0),
        (-89.0, 1.5),
        (89.5, 180.0),
        (90.0, 180.0),
        (100.0, 180.0),
    ]
        @test Intp.interpolate(flat, (x,)) == expected
    end
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
    ret = Intp.regrid!(dest_data, regridder, dest_grid)
    @test ret === dest_data
    @test dest_data == expected

    # OrderedDict adapters
    dest_dict = OrderedDict("lon" => dest_grid[1], "lat" => dest_grid[2])
    @test Intp.regrid(regridder, dest_dict) == expected
    fill!(dest_data, 0.0)
    Intp.regrid!(dest_data, regridder, dest_dict)
    @test dest_data == expected

    # Size and dimension checking
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

    # The destination coordinates are arbitrary: unlike the source coordinates,
    # they do not have to be sorted or distinct
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
    # Interpolation happens in the float type of the source data, so regrid!
    # gives the same values as regrid and interpolate no matter what the eltype
    # of the destination array is
    src = Float32[1.0, 1.0000001, 2.0]
    regridder =
        Intp.Regridder(src, ([1.0, 2.0, 3.0],), (Intp.Flat(),), (Intp.Node(),))

    expected = Intp.interpolate(regridder, (1.5,))
    @test expected === Intp.regrid(regridder, ([1.5],))[1]

    dest64 = zeros(1)
    Intp.regrid!(dest64, regridder, ([1.5],))
    @test dest64[1] == expected

    dest32 = zeros(Float32, 1)
    Intp.regrid!(dest32, regridder, ([1.5],))
    @test dest32[1] === expected
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

    # The output eltype comes from the eltype of the source data
    data32 = Float32[3.0, 1.0, 0.0]
    regridder32 = Intp.Regridder(
        data32,
        ([1.0, 2.0, 3.0],),
        (Intp.Throw(),),
        (Intp.Node(),),
    )
    @test Intp.interpolate(regridder32, (1.5,)) === 2.0f0
    @test eltype(Intp.regrid(regridder32, ([1.0, 1.5, 2.0],))) == Float32

    # The weights are computed in the float type of the coordinates and the
    # result in the float type of the data, so a query that is wider than the
    # coordinates does not widen either one
    regridder_mixed = Intp.Regridder(
        Float32[1.0, 2.0, 4.0],
        ([1.0f0, 2.0f0, 3.0f0],),
        (Intp.Flat(),),
        (Intp.Node(),),
    )
    @test regridder_mixed.stencil_cache isa
          NTuple{1, Vector{Tuple{Int, Int, Float32}}}
    @test Intp.interpolate(regridder_mixed, (2.5,)) === 3.0f0
    @test eltype(Intp.regrid(regridder_mixed, ([2.5],))) == Float32

    # Integer data produces floats
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

    # Union{Missing, T} eltypes work; the eltype of the result follows the
    # eltype of the source data and not whether a missing value is present
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
    @test regridder_m.stencil_cache isa
          NTuple{2, Vector{Tuple{Int, Int, Float64}}}
    @test Intp.interpolate(regridder_m, (1.5, 4.5)) === 3.0
    out = Intp.regrid(regridder_m, (lon_m, [4.0, 4.5, 5.0]))
    @test eltype(out) == Union{Missing, Float64}
    @test out == [1.0 2.5 4.0; 2.0 3.5 5.0; 3.0 4.5 6.0]

    # NaNLinearInterpolation never returns missing, so Missing is dropped
    regridder_m_nan = Intp.Regridder(
        data_m,
        (lon_m, [4.0, 5.0, 6.0]),
        (Intp.Throw(), Intp.Throw()),
        (Intp.Node(), Intp.Node()),
        Intp.NaNLinearInterpolation(0.0),
    )
    out = Intp.regrid(regridder_m_nan, (lon_m, [4.0, 4.5, 5.0]))
    @test eltype(out) == Float64
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

@testset "NaN-aware interpolation" begin
    @test_throws "is not between 0 and 1" Intp.Regridder(
        [1.0, 2.0],
        ([1.0, 2.0],),
        (Intp.Throw(),),
        (Intp.Node(),),
        Intp.NaNLinearInterpolation(1.5),
    )

    data = [10.0 NaN; 30.0 40.0]
    grid = ([1.0, 2.0], [1.0, 2.0])
    extp = (Intp.Throw(), Intp.Throw())
    stag = (Intp.Node(), Intp.Node())

    # With max_nan_fraction = 1, the result is the weighted average of the
    # non-NaN points
    avg =
        Intp.Regridder(data, grid, extp, stag, Intp.NaNLinearInterpolation(1.0))
    @test Intp.interpolate(avg, (1.5, 1.5)) ≈ 80.0 / 3.0
    @test Intp.interpolate(avg, (1.0, 1.0)) == 10.0

    # With max_nan_fraction = 0, any NaN with a nonzero weight produces NaN
    strict =
        Intp.Regridder(data, grid, extp, stag, Intp.NaNLinearInterpolation(0.0))
    @test isnan(Intp.interpolate(strict, (1.5, 1.5)))
    @test Intp.interpolate(strict, (1.0, 1.0)) == 10.0

    # All-NaN stencils produce NaN for every query point, including points
    # whose stencil weights do not sum exactly to one in floating point
    xs = collect(range(0.0, 1.0, length = 11))
    all_nan = Intp.Regridder(
        fill(NaN, 11, 11),
        (xs, xs),
        extp,
        stag,
        Intp.NaNLinearInterpolation(1.0),
    )
    pts = [(i / 300.0, j / 300.0) for i in 0:300 for j in 0:300]
    @test all(isnan, Intp.interpolate(all_nan, pts))

    # Missing values are treated like NaN
    data_m = convert(Matrix{Union{Missing, Float64}}, [10.0 20.0; 30.0 40.0])
    data_m[1, 2] = missing
    missing_as_nan = Intp.Regridder(
        data_m,
        grid,
        extp,
        stag,
        Intp.NaNLinearInterpolation(1.0),
    )
    @test Intp.interpolate(missing_as_nan, (1.5, 1.5)) ≈ 80.0 / 3.0
end
