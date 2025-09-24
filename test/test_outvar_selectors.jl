using Test
import ClimaAnalysis

import OrderedCollections: OrderedDict
import Dates

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

@testset "Selectors" begin
    lon = [-20.0, -10.0, 0.0, 10.0]
    lat = [-30.0, 1.0]
    var =
        TemplateVar() |>
        add_dim("lon", lon, units = "deg") |>
        add_dim("lat", lat, units = "deg") |>
        add_attribs(long_name = "hi") |>
        initialize

    var_2d_dim =
        TemplateVar() |>
        add_dim("2ddim", [[1.0, 2.0] [3.0, 4.0]], units = "idk") |>
        add_attribs(two_dims = "yes") |>
        initialize

    @testset "NearestValue" begin
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            -12.0,
            ClimaAnalysis.NearestValue(),
        ) == 2
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            -10.0,
            ClimaAnalysis.NearestValue(),
        ) == 2
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            -15.0,
            ClimaAnalysis.NearestValue(),
        ) == 1
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            100.0,
            ClimaAnalysis.NearestValue(),
        ) == 4
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            -100.0,
            ClimaAnalysis.NearestValue(),
        ) == 1
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var_2d_dim,
            "2ddim",
            1,
            ClimaAnalysis.NearestValue(),
        )
    end

    @testset "MatchValue" begin
        @test ClimaAnalysis.Var.get_index(
            var,
            "lat",
            -30.0,
            ClimaAnalysis.MatchValue(),
        ) == 1
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var,
            "lat",
            100.0,
            ClimaAnalysis.MatchValue(),
        )
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var,
            "lat",
            -100.0,
            ClimaAnalysis.MatchValue(),
        )
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var_2d_dim,
            "2ddim",
            1.0,
            ClimaAnalysis.MatchValue(),
        )
    end

    @testset "Index" begin
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            1,
            ClimaAnalysis.Index(),
        ) == 1
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            2,
            ClimaAnalysis.Index(),
        ) == 2
        @test ClimaAnalysis.Var.get_index(
            var,
            "lon",
            4,
            ClimaAnalysis.Index(),
        ) == 4
        @test ClimaAnalysis.Var.get_index(
            var_2d_dim,
            "2ddim",
            3,
            ClimaAnalysis.Index(),
        ) == 3
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var,
            "lon",
            5,
            ClimaAnalysis.Index(),
        )
        @test_throws ErrorException ClimaAnalysis.Var.get_index(
            var,
            "lon",
            0,
            ClimaAnalysis.Index(),
        ) == 1
    end
end

@testset "Slicing" begin
    z = 0.0:20.0 |> collect
    time = 100.0:110.0 |> collect

    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("z", z, b = 2) |>
        add_attribs(long_name = "hi") |>
        initialize
    data = var.data

    z_sliced = ClimaAnalysis.slice(var, z = 1.0)
    # 1.0 is the second index
    z_expected_data = data[:, 2]
    @test z_sliced.dims == OrderedDict(["time" => time])
    @test z_sliced.dim_attributes ==
          OrderedDict(["time" => Dict("units" => "s")])
    @test z_sliced.data == z_expected_data

    t_sliced = ClimaAnalysis.slice(var, time = 200.0)
    # 200 is the last index
    t_expected_data = data[end, :]
    @test t_sliced.dims == OrderedDict(["z" => z])
    @test t_sliced.dim_attributes == OrderedDict(["z" => Dict("b" => 2)])
    @test t_sliced.data == t_expected_data

    @test t_sliced.attributes["long_name"] == "hi time = 1m 50.0s"

    # Test with the general slice

    t_sliced = ClimaAnalysis.slice(var, time = 200.0)
    # 200 is the last index
    t_expected_data = data[end, :]
    @test t_sliced.dims == OrderedDict(["z" => z])
    @test t_sliced.dim_attributes == OrderedDict(["z" => Dict("b" => 2)])
    @test t_sliced.data == t_expected_data

    @test t_sliced.attributes["long_name"] == "hi time = 1m 50.0s"

    @test t_sliced.attributes["slice_time"] == "110.0"
    @test t_sliced.attributes["slice_time_units"] == "s"

    # Test with different name for time
    t_sliced_other_name = ClimaAnalysis.slice(var, t = 200.0)
    @test t_sliced_other_name.attributes == t_sliced.attributes
    @test t_sliced_other_name.data == t_sliced.data

    # Test with Dates.DateTime

    var_with_start_date = copy(var)
    push!(var_with_start_date.attributes, "start_date" => "2001-1-1")
    t_sliced_dates = ClimaAnalysis.slice(
        var_with_start_date,
        t = Dates.DateTime(2001) + Dates.Second(200),
    )
    @test t_sliced.data == t_sliced_dates.data
    # Dimension is not time
    @test_throws ErrorException ClimaAnalysis.slice(
        var_with_start_date,
        z = Dates.DateTime(2001),
    )

    # Test slicing with MatchValue and Index
    time = 100.0:110.0 |> collect
    z = 0.0:20.0 |> collect
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("z", z, b = 2) |>
        add_attribs(long_name = "hi", start_date = "2010-01-01") |>
        initialize
    data = var.data

    # Slicing by approximately matching with the value
    sliced_var =
        ClimaAnalysis.slice(var, t = 105.0, by = ClimaAnalysis.MatchValue())
    @test sliced_var.data == data[6, :]
    @test_throws ErrorException ClimaAnalysis.slice(
        var,
        t = 105.5,
        by = ClimaAnalysis.MatchValue(),
    )

    # Also test with dates
    sliced_var = ClimaAnalysis.slice(
        var,
        time = Dates.DateTime(2010, 1, 1, 0, 1, 42),
        by = ClimaAnalysis.MatchValue(),
    )
    @test sliced_var.data == data[3, :]

    # Slicing by index
    sliced_var = ClimaAnalysis.slice(var, t = 3, by = ClimaAnalysis.Index())
    @test sliced_var.data == data[3, :]
    @test_throws ErrorException ClimaAnalysis.slice(
        var,
        t = 101,
        by = ClimaAnalysis.Index(),
    )

    # Test index with OutputVar with a dimension that is 2d
    var_2d_dim =
        TemplateVar() |>
        add_dim("dim2d", [[1.0, 2.0] [3.0, 4.0]], units = "idk") |>
        add_attribs(two_dims = "yes") |>
        initialize
    sliced_var =
        ClimaAnalysis.slice(var_2d_dim, dim2d = 3, by = ClimaAnalysis.Index())
    @test sliced_var.data == fill(3)
end

@testset "Windowing" begin
    z = 0.0:20.0 |> collect
    time = 0.0:10.0 |> collect
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("z", z, b = 2) |>
        add_attribs(long_name = "hi") |>
        initialize
    data = var.data

    # Dimension not existing
    @test_throws ErrorException ClimaAnalysis.window(var, "lat")

    # Left right not ordered
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "time",
        left = 10,
        right = 1,
    )

    var_windowed = ClimaAnalysis.window(var, "time", left = 2.5, right = 5.1)
    expected_data = data[3:6, :]

    @test var_windowed.data == expected_data

    @test var_windowed.dims["time"] == time[3:6]

    # Test with different name for time
    var_windowed_t = ClimaAnalysis.window(var, "t", left = 2.5, right = 5.1)
    @test var_windowed_t.attributes == var_windowed.attributes
    @test var_windowed_t.data == var_windowed.data

    # Windowing with Dates.DateTime

    # First, `start_date` not available
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "time",
        left = Dates.DateTime(2001),
    )

    # Now with a valid start date
    var_with_start_date = copy(var)
    push!(var_with_start_date.attributes, "start_date" => "2001-1-1")
    var_windowed_dates = ClimaAnalysis.window(
        var_with_start_date,
        "time",
        left = Dates.DateTime(2001, 1, 1, 0, 0, 2, 500),
        right = Dates.DateTime(2001, 1, 1, 0, 0, 5, 100),
    )
    @test var_windowed.data == var_windowed_dates.data

    # Dates with a dimension that is not time
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "z",
        left = Dates.DateTime(2001),
    )

    # Test with MatchValue and Index
    time = 100.0:110.0 |> collect
    z = 0.0:20.0 |> collect
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("z", z, b = 2) |>
        add_attribs(long_name = "hi") |>
        initialize
    data = var.data

    # Windowing by approximately matching with the values
    windowed_var = ClimaAnalysis.window(
        var,
        "time",
        left = 105.0,
        right = 108.0,
        by = ClimaAnalysis.MatchValue(),
    )
    @test windowed_var.data == data[6:9, :]
    windowed_var = ClimaAnalysis.window(
        var,
        "time",
        left = 105.0,
        by = ClimaAnalysis.MatchValue(),
    )
    @test windowed_var.data == data[6:end, :]
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "time",
        left = 102.0,
        right = 103.5,
        by = ClimaAnalysis.MatchValue(),
    )

    # Windowing by index
    windowed_var = ClimaAnalysis.window(
        var,
        "time",
        left = 2,
        right = 5,
        by = ClimaAnalysis.Index(),
    )
    @test windowed_var.data == data[2:5, :]
    windowed_var =
        ClimaAnalysis.window(var, "time", right = 7, by = ClimaAnalysis.Index())
    @test windowed_var.data == data[begin:7, :]
    @test_throws ErrorException ClimaAnalysis.window(
        var,
        "time",
        left = 7.5,
        right = 8,
        by = ClimaAnalysis.Index(),
    )
end

@testset "Selecting" begin
    time = 0.0:7.0 |> collect
    lon = 0.0:4.0 |> collect
    lat = 40.0:55.0 |> collect
    var =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("lat", lat, units = "degrees") |>
        add_attribs(long_name = "hi", start_date = "2010-1-1") |>
        initialize

    # Select 1 dimension
    # With values
    nearest_var1 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.NearestValue(),
        lat = [20, 60],
    )
    @test nearest_var1.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[[begin, end]])
    @test isequal(nearest_var1.data, var.data[:, :, [begin, end]])
    @test ndims(nearest_var1.data) == 3

    nearest_var2 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.NearestValue(),
        latitude = 60,
    )
    @test nearest_var2.dims == OrderedDict("time" => time, "lon" => lon)
    @test isequal(nearest_var2.data, var.data[:, :, end])
    @test ndims(nearest_var2.data) == 2

    match_var1 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.MatchValue(),
        latitude = [55.0],
    )
    @test match_var1.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[[end]])
    @test isequal(match_var1.data, var.data[:, :, [end]])
    @test ndims(match_var1.data) == 3

    match_var2 =
        ClimaAnalysis.select(var, by = ClimaAnalysis.MatchValue(), lat = 53:55)
    @test match_var2.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[[14, 15, 16]])
    @test isequal(match_var2.data, var.data[:, :, [14, 15, 16]])
    @test ndims(match_var2.data) == 3

    # With indices
    index_var1 = ClimaAnalysis.select(var, by = ClimaAnalysis.Index(), lat = 1)
    @test isequal(index_var1.data, var.data[:, :, 1])
    @test index_var1.dims == OrderedDict("time" => time, "lon" => lon)
    @test ndims(index_var1.data) == 2

    index_var2 =
        ClimaAnalysis.select(var, by = ClimaAnalysis.Index(), lat = [1])
    @test index_var2.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[[1]])
    @test isequal(index_var2.data, var.data[:, :, [1]])
    @test ndims(index_var2.data) == 3

    index_var3 =
        ClimaAnalysis.select(var, by = ClimaAnalysis.Index(), lat = [1, 5])
    @test index_var3.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[[1, 5]])
    @test isequal(index_var3.data, var.data[:, :, [1, 5]])
    @test ndims(index_var3.data) == 3

    index_var4 =
        ClimaAnalysis.select(var, by = ClimaAnalysis.Index(), lat = 1:3)
    @test index_var4.dims ==
          OrderedDict("time" => time, "lon" => lon, "lat" => lat[1:3])
    @test isequal(index_var4.data, var.data[:, :, 1:3])
    @test ndims(index_var4.data) == 3

    # With dates
    dates_var1 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.MatchValue(),
        time = Dates.DateTime(2010),
    )
    @test dates_var1.dims == OrderedDict("lon" => lon, "lat" => lat)
    @test isequal(dates_var1.data, var.data[1, :, :])
    @test ndims(dates_var1.data) == 2

    dates_var2 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.MatchValue(),
        time = [Dates.DateTime(2010)],
    )
    @test dates_var2.dims ==
          OrderedDict("time" => time[[1]], "lon" => lon, "lat" => lat)
    @test isequal(dates_var2.data, var.data[[1], :, :])
    @test ndims(dates_var2.data) == 3

    dates_var3 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.MatchValue(),
        time = [Dates.DateTime(2010), Dates.DateTime(2010, 1, 1, 0, 0, 2)],
    )
    @test dates_var3.dims ==
          OrderedDict("time" => time[[1, 3]], "lon" => lon, "lat" => lat)
    @test isequal(dates_var3.data, var.data[[1, 3], :, :])
    @test ndims(dates_var3.data) == 3

    dates_var4 = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.MatchValue(),
        time = Dates.DateTime(2010, 1, 1, 0, 0, 0):Dates.Second(2):Dates.DateTime(
            2010,
            1,
            1,
            0,
            0,
            7,
        ),
    )
    @test dates_var4.dims ==
          OrderedDict("time" => time[[1, 3, 5, 7]], "lon" => lon, "lat" => lat)
    @test isequal(dates_var4.data, var.data[[1, 3, 5, 7], :, :])
    @test ndims(dates_var4.data) == 3

    # Duplicate indices / values
    duplicate_var = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        lat = [1, 2, 2, 3, 3, 5],
    )
    @test isequal(duplicate_var.data, var.data[:, :, [1, 2, 2, 3, 3, 5]])

    # Select multiple dimensions
    var_select_3d = ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        lat = [3, 4, 5],
        time = [1, 2],
        lon = 2,
    )
    @test var_select_3d.dims ==
          OrderedDict("time" => time[[1, 2]], "lat" => lat[[3, 4, 5]])
    @test isequal(var_select_3d.data, var.data[[1, 2], 2, [3, 4, 5]])
    @test ndims(var_select_3d.data) == 2

    # 2D dim for a OutputVar
    var_2d_dim =
        TemplateVar() |>
        add_dim("dim2d", [[1.0, 2.0, 3.0] [4.0, 5.0, 6.0]], units = "idk") |>
        add_attribs(two_dims = "yes") |>
        add_data(data = [7.0, 8.0, 9.0]) |>
        initialize
    var_2d_dim_selected = ClimaAnalysis.select(
        var_2d_dim,
        by = ClimaAnalysis.Index(),
        dim2d = [1],
    )
    @test isequal(var_2d_dim.data[[1]], var_2d_dim_selected.data)
    @test isequal(
        var_2d_dim_selected.dims["dim2d"],
        var_2d_dim.dims["dim2d"][:, [1]],
    )
    @test ndims(var_2d_dim_selected.data) == 1

    # Scalar and vector indexing with one dimension OutputVar
    var_1d =
        TemplateVar() |>
        add_dim("time", time, units = "s") |>
        add_attribs(long_name = "hi", start_date = "2010-1-1") |>
        initialize

    var_vector_index =
        ClimaAnalysis.select(var_1d, by = ClimaAnalysis.Index(), time = [1])
    @test var_vector_index.dims == OrderedDict("time" => time[[1]])
    @test isequal(var_vector_index.data, [1.0])
    @test ndims(var_vector_index.data) == 1

    var_scalar_index =
        ClimaAnalysis.select(var_1d, by = ClimaAnalysis.Index(), time = 1)
    @test isempty(var_scalar_index.dims)
    @test isequal(var_scalar_index.data, fill(1.0))
    @test ndims(var_scalar_index.data) == 0

    # Test view_select gets a view and select does not get a view
    var_copy = ClimaAnalysis.remake(var)
    view_var = ClimaAnalysis.view_select(
        var_copy,
        by = ClimaAnalysis.Index(),
        lat = [3, 4, 5],
    )
    no_view_var = ClimaAnalysis.select(
        var_copy,
        by = ClimaAnalysis.Index(),
        lat = [3, 4, 5],
    )
    @test isequal(view_var.data, no_view_var.data)
    view_var.data[1, 1, 1] = 42.0
    @test var_copy.data[1, 1, 3] == 42.0
    @test !isequal(view_var.data, no_view_var.data)

    # Error handling
    # Invalid indices
    @test_throws ErrorException ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        cool = 42,
    )
    @test_throws ErrorException ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        cool = 1:42,
    )
    @test_throws ErrorException ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        cool = [1, 42],
    )

    # Missing dimensions
    @test_throws ErrorException ClimaAnalysis.select(
        var,
        by = ClimaAnalysis.Index(),
        cool = 42,
    )
end
