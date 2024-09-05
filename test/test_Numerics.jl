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

    # Integration weights for lon (not equispaced)
    lon = collect(range(-180.0, 180.0, 5))
    lon_weights = [90.0 for _ in lon] .* (π / 180.0)
    @test all(
        isapprox.(
            lon_weights,
            ClimaAnalysis.Numerics._integration_weights_lon_equispaced(lon),
        ),
    )

    lat = collect(range(-90.0, 90.0, 5))
    lat_weights = [45.0 for _ in lat] .* (π / 180.0) .* cosd.(lat)
    @test all(
        isapprox.(
            lat_weights,
            ClimaAnalysis.Numerics._integration_weights_lat_equispaced(lat),
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

    # Integrating only lon (non equispaced)
    lon = collect(range(-179.5, 179.5, 360))
    lon_data = ones(length(lon))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lon(lon_data, lon, dims = 1)[1],
        2.0π,
        atol = 0.01,
    )

    # Integrating only lat (non equispaced)
    lat = collect(range(-89.5, 89.5, 180))
    lat_data = ones(length(lat))
    @test isapprox(
        ClimaAnalysis.Numerics._integrate_lat(lat_data, lat, dims = 1)[1],
        2.0,
        atol = 0.01,
    )

end
