using Test
import ClimaAnalysis

import Dates
import Statistics: mean, std

import ClimaAnalysis.Template:
    TemplateVar,
    make_template_var,
    add_attribs,
    add_dim,
    add_time_dim,
    add_lon_dim,
    add_lat_dim,
    add_data,
    ones_data,
    zeros_data,
    one_to_n_data,
    initialize

@testset "GroupAll and Reduce" begin
    time = [0.0, 1.0, 2.0, 3.0]
    lat = [-90.0, 0.0, 90.0]
    lon = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_attribs(start_date = Dates.DateTime(2010), short_name = "lwu") |>
        initialize

    dimnames = ClimaAnalysis.dim_names(var)
    reductions = [sum, prod, minimum, maximum, mean, std]
    for reduction in reductions
        for (i, dimname) in enumerate(dimnames)
            reduced_var =
                var |>
                ClimaAnalysis.GroupAll(dimname) |>
                ClimaAnalysis.Reduce(reduction) |>
                ClimaAnalysis.combine
            @test reduced_var.data == reduction(var.data, dims = i)
            @test length(reduced_var.dims[dimname]) == 1
            @test reduced_var.dims[dimname][1] == first(var.dims[dimname])
            @test reduced_var.attributes == var.attributes
            @test reduced_var.dim_attributes == var.dim_attributes
        end
    end

    @test_throws ErrorException var |>
                                ClimaAnalysis.GroupAll("cool_dim") |>
                                ClimaAnalysis.Reduce(mean) |>
                                ClimaAnalysis.combine
end

@testset "SplitBySeason and Reduce" begin
    time = [0.0, 1.0, 2.0, 3.0]
    lat = [-90.0, 0.0, 90.0]
    lon = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]
    single_season_var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_attribs(start_date = Dates.DateTime(2010)) |>
        initialize

    # OutputVar with a single season
    reduced_var =
        single_season_var |>
        ClimaAnalysis.SplitSeason() |>
        ClimaAnalysis.Reduce(mean) |>
        ClimaAnalysis.combine
    @test reduced_var.data == mean(single_season_var.data, dims = 1)

    # OutputVar with multiple seasons
    # DJF - 12, 1, 2
    # MAM - 3, 4, 5
    # JJA - 6, 7, 8
    # SON - 9, 10, 11
    ref_date = Dates.DateTime(2010)
    time =
        ClimaAnalysis.Utils.date_to_time.(
            Ref(ref_date),
            [
                Dates.DateTime(2010, 1),
                Dates.DateTime(2010, 2),
                Dates.DateTime(2010, 4),
                Dates.DateTime(2010, 6),
                Dates.DateTime(2010, 7),
                Dates.DateTime(2010, 8),
                Dates.DateTime(2011, 10),
            ],
        )
    lat = [-90.0, 0.0, 90.0]
    lon = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]
    multiple_season_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("time", time, units = "s") |>
        add_dim("lon", lon, units = "degrees") |>
        add_attribs(start_date = ref_date) |>
        initialize

    reduced_var =
        multiple_season_var |>
        ClimaAnalysis.SplitSeason() |>
        ClimaAnalysis.Reduce(mean) |>
        ClimaAnalysis.combine
    @test reduced_var.data == cat(
        mean(multiple_season_var.data[:, 1:2, :], dims = 2), # DJF
        mean(multiple_season_var.data[:, 3:3, :], dims = 2), # JJA
        mean(multiple_season_var.data[:, 4:6, :], dims = 2), # MAM
        mean(multiple_season_var.data[:, 7:7, :], dims = 2), # SON
        dims = 2,
    )
    @test length(reduced_var.dims["time"]) == 4
    @test reduced_var.dims["time"] ==
          multiple_season_var.dims["time"][[1, 3, 4, 7]]
    @test reduced_var.attributes == multiple_season_var.attributes
    @test reduced_var.dim_attributes == multiple_season_var.dim_attributes

    # OutputVar with dates out of order
    diff_order_var = ClimaAnalysis.select(
        multiple_season_var,
        by = ClimaAnalysis.Index(),
        time = [3, 4, 1, 6, 2, 5, 7],
    )
    reduced_var2 =
        diff_order_var |>
        ClimaAnalysis.SplitSeason() |>
        ClimaAnalysis.Reduce(mean) |>
        ClimaAnalysis.combine
    @test reduced_var.attributes == reduced_var2.attributes
    @test reduced_var.dims == reduced_var2.dims
    @test reduced_var.dim_attributes == reduced_var2.dim_attributes
    @test reduced_var.data == reduced_var2.data

    # OutputVar with repeated dates
    repeated_dates_var =
        TemplateVar() |>
        add_dim("time", [0.0, 0.0, 1.0], units = "s") |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_attribs(start_date = Dates.DateTime(2010)) |>
        initialize

    @test_throws ErrorException repeated_dates_var |>
                                ClimaAnalysis.SplitSeason() |>
                                ClimaAnalysis.Reduce(mean) |>
                                ClimaAnalysis.combine

    # Time dimension does not exist
    no_time_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        initialize

    @test_throws ErrorException no_time_var |>
                                ClimaAnalysis.SplitSeason() |>
                                ClimaAnalysis.Reduce(mean) |>
                                ClimaAnalysis.combine
end

@testset "Invalid split and apply operations chaining" begin
    time = [0.0, 1.0]
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_attribs(start_date = Dates.DateTime(2010), short_name = "lwu") |>
        initialize

    @test_throws r"A split operation is already set" var |>
                                                     ClimaAnalysis.GroupAll(
                                                         "time",
                                                     ) |>
                                                     ClimaAnalysis.GroupAll(
                                                         "time",
                                                     )

    @test_throws r"An apply operation is already set" var |>
                                                      ClimaAnalysis.GroupAll(
                                                          "time",
                                                      ) |>
                                                      ClimaAnalysis.Reduce(
                                                          mean,
                                                      ) |>
                                                      ClimaAnalysis.Reduce(mean)
end
