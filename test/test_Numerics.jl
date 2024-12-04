using Test
import ClimaAnalysis

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

@testset "Get indices and sign" begin
    indices = ((1, 2),)
    @test ClimaAnalysis.Numerics.get_indices(indices, 0) == (2,)
    @test ClimaAnalysis.Numerics.get_indices(indices, 1) == (1,)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 0) == (1,)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 1) == (2,)
    @test ClimaAnalysis.Numerics.get_sign(indices, 0) == -1
    @test ClimaAnalysis.Numerics.get_sign(indices, 1) == 1

    indices = ((1, 2), (3, 4))
    @test ClimaAnalysis.Numerics.get_indices(indices, 0) == (2, 4)
    @test ClimaAnalysis.Numerics.get_indices(indices, 1) == (1, 4)
    @test ClimaAnalysis.Numerics.get_indices(indices, 2) == (2, 3)
    @test ClimaAnalysis.Numerics.get_indices(indices, 3) == (1, 3)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 0) == (1, 3)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 1) == (2, 3)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 2) == (1, 4)
    @test ClimaAnalysis.Numerics.get_complement_indices(indices, 3) == (2, 4)
    @test ClimaAnalysis.Numerics.get_sign(indices, 0) == 1
    @test ClimaAnalysis.Numerics.get_sign(indices, 1) == -1
    @test ClimaAnalysis.Numerics.get_sign(indices, 2) == -1
    @test ClimaAnalysis.Numerics.get_sign(indices, 3) == 1
end

@testset "Find indices for cell" begin
    val1 = 5
    ax1 = [0, 10]
    @test ClimaAnalysis.Numerics.find_cell_indices_for_ax(val1, ax1) == (1, 2)

    val2 = 6
    ax2 = [0, 4, 10]
    @test ClimaAnalysis.Numerics.find_cell_indices_for_ax(val2, ax2) == (2, 3)

    @test ClimaAnalysis.Numerics.find_cell_indices_for_axes(
        (val1, val2),
        (ax1, ax2),
    ) == ((1, 2), (2, 3))
end

@testset "Extrapolation conditions" begin
    throw = ClimaAnalysis.Numerics.extp_cond_throw()
    flat = ClimaAnalysis.Numerics.extp_cond_flat()
    periodic = ClimaAnalysis.Numerics.extp_cond_periodic()

    ax = [0, 1, 2, 3]
    @test_throws ErrorException throw.get_val_for_point(10, ax)
    @test_throws ErrorException throw.get_val_for_point(-1, ax)
    @test throw.get_val_for_point(1.5, ax) == 1.5

    @test flat.get_val_for_point(10, ax) == 3
    @test flat.get_val_for_point(-1, ax) == 0
    @test flat.get_val_for_point(1.5, ax) == 1.5

    @test periodic.get_val_for_point(10, ax) == 1
    @test periodic.get_val_for_point(-1, ax) == 2
    @test periodic.get_val_for_point(1.5, ax) == 1.5
    @test periodic.get_val_for_point(3, ax) == 3
end

@testset "Extrapolate to new point" begin
    throw = ClimaAnalysis.Numerics.extp_cond_throw()
    flat = ClimaAnalysis.Numerics.extp_cond_flat()
    periodic = ClimaAnalysis.Numerics.extp_cond_periodic()

    ax1 = [0, 1, 2, 3]
    @test ClimaAnalysis.Numerics.extp_to_point((1,), (ax1,), (throw,)) == (1,)

    ax2 = [4, 5, 6, 7]
    @test ClimaAnalysis.Numerics.extp_to_point(
        (-1, 8),
        (ax1, ax2),
        (flat, periodic),
    ) == (0, 5)
end

@testset "Interpolation" begin
    throw = ClimaAnalysis.Numerics.extp_cond_throw()
    flat = ClimaAnalysis.Numerics.extp_cond_flat()
    periodic = ClimaAnalysis.Numerics.extp_cond_periodic()

    # 1D case
    axes = ([1.0, 2.0, 3.0],)
    data = [3.0, 1.0, 0.0]

    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.0,),
        axes,
        data,
        (throw,),
    ) == 3.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (3.0,),
        axes,
        data,
        (throw,),
    ) == 0.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5,),
        axes,
        data,
        (throw,),
    ) == 2.0

    # 1D case with extrapolation conditions
    @test_throws ErrorException ClimaAnalysis.Numerics.linear_interpolate(
        (0.0,),
        axes,
        data,
        (throw,),
    )
    @test_throws ErrorException ClimaAnalysis.Numerics.linear_interpolate(
        (4.0,),
        axes,
        data,
        (throw,),
    )
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (0.0,),
        axes,
        data,
        (flat,),
    ) == 3.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (4.0,),
        axes,
        data,
        (flat,),
    ) == 0.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (0.0,),
        axes,
        data,
        (periodic,),
    ) == 1.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (4.0,),
        axes,
        data,
        (periodic,),
    ) == 1.0

    # 2D case
    axes = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    data = reshape(1:9, (3, 3))

    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.0, 4.0),
        axes,
        data,
        (throw, throw),
    ) == 1.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (3.0, 6.0),
        axes,
        data,
        (throw, throw),
    ) == 9.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (2.0, 5.0),
        axes,
        data,
        (throw, throw),
    ) == 5.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5, 4.5),
        axes,
        data,
        (throw, throw),
    ) == 3.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5, 5.5),
        axes,
        data,
        (throw, throw),
    ) == 6.0

    # 2D cases with extrapolation conditions
    @test_throws ErrorException ClimaAnalysis.Numerics.linear_interpolate(
        (4.0, 5.0),
        axes,
        data,
        (throw, flat),
    )
    @test_throws ErrorException ClimaAnalysis.Numerics.linear_interpolate(
        (2.0, 7.0),
        axes,
        data,
        (flat, throw),
    )
    @test_throws ErrorException ClimaAnalysis.Numerics.linear_interpolate(
        (0.0, 8.0),
        axes,
        data,
        (throw, throw),
    )
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (0.0, 8.0),
        axes,
        data,
        (flat, flat),
    ) == 7.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (4.0, 7.0),
        axes,
        data,
        (periodic, periodic),
    ) == 5.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (3.0, 6.0),
        axes,
        data,
        (periodic, periodic),
    ) == 9.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (4.0, 7.0),
        axes,
        data,
        (flat, periodic),
    ) == 6.0

    # 3D cases with extrapolation conditions
    axes = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0])
    data = reshape(1:27, (3, 3, 3))
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.0, 5.0, 7.0),
        axes,
        data,
        (throw, throw, throw),
    ) == 4.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5, 5.2, 7.5),
        axes,
        data,
        (throw, throw, throw),
    ) ≈ 9.6

    # Non equispaced
    axes = ([1.0, 3.0, 7.0], [4.0, 5.0, 7.0])
    data = reshape(1:9, (3, 3))
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (2.0, 4.5),
        axes,
        data,
        (throw, throw, throw),
    ) == 3.0
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (5.0, 6.0),
        axes,
        data,
        (throw, throw, throw),
    ) == 7.0

    # Different types
    # Axes have different types and inputs have different types
    axes = ([1.0f0, 2.0f0], [Float16(3.0), Float16(4.0)])
    data = [[1.0, 2.0] [3.0, 4.0]]
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5f0, 3.5),
        axes,
        data,
        (throw, throw),
    ) == 2.5
    @test ClimaAnalysis.Numerics.linear_interpolate(
        (1.5, 4.5f0),
        axes,
        data,
        (flat, flat),
    ) == 3.5
    @test ClimaAnalysis.Numerics.linear_interpolate(
        [1.5, 4.5f0],
        axes,
        data,
        (flat, flat),
    ) == 3.5

    # Single number
    axes = ([1.0, 2.0, 3.0],)
    data = [3.0, 1.0, 0.0]

    @test ClimaAnalysis.Numerics.linear_interpolate(
        1.0,
        axes,
        data,
        (throw,),
    ) == 3.0
end
