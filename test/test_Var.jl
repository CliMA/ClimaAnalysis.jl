using Test
import ClimaAnalysis

import Interpolations as Intp
import Statistics: mean
import OrderedCollections: OrderedDict

@testset "General" begin
    # Add test for short constructor
    long = -180.0:180.0 |> collect
    data = copy(long)

    longvar = ClimaAnalysis.OutputVar(Dict("long" => long), data)

    @test longvar.dims["long"] == long

    # center_longitude!
    #
    # Check var without long
    dims = Dict("z" => long)
    var_error = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        dims,
        Dict{String, Any}(),
        data,
    )
    @test_throws ErrorException ClimaAnalysis.center_longitude!(
        var_error,
        180.0,
    )

    time = 0:10.0 |> collect
    dims = OrderedDict("lon" => long, "time" => time)
    data = collect(reshape(1:(361 * 11), (361, 11)))
    var_good = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        dims,
        Dict{String, Any}(),
        data,
    )
    ClimaAnalysis.center_longitude!(var_good, 90.0)
    # We are shifting by 91
    @test var_good.dims["lon"][180] == 90
    @test var_good.data[3, :] == data[3, :]
    @test var_good.data[180, 1] == 271

    @test_throws ErrorException ClimaAnalysis.OutputVar(
        Dict("time" => time),
        [1],
    )
end

@testset "Arithmetic operations" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data1 = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    var1 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data1)

    dim_attributes2 = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("lol" => 2),
        "lat" => Dict("a" => 1),
    ])

    var2 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes2, data1)

    data3 = 5.0 .+ collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))
    attribs3 = Dict("long_name" => "bob", "short_name" => "bula")
    var3 = ClimaAnalysis.OutputVar(attribs3, dims, dim_attributes, data3)

    # Check arecompatible
    @test !ClimaAnalysis.arecompatible(var1, var2)
    @test ClimaAnalysis.arecompatible(var1, var3)

    var1plus10 = var1 + 10

    @test var1plus10.data == data1 .+ 10
    @test ClimaAnalysis.short_name(var1plus10) == "bob + 10"
    @test ClimaAnalysis.long_name(var1plus10) == "hi + 10"

    tenplusvar1 = 10 + var1

    @test tenplusvar1.data == data1 .+ 10
    @test ClimaAnalysis.short_name(tenplusvar1) == "10 + bob"
    @test ClimaAnalysis.long_name(tenplusvar1) == "10 + hi"

    var1plusvar3 = var1 + var3

    @test var1plusvar3.data == data1 .+ data3
    @test ClimaAnalysis.short_name(var1plusvar3) == "bob + bula"
    @test ClimaAnalysis.long_name(var1plusvar3) == "hi + bob"
end

@testset "Reductions (sphere dims)" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    # Test copy
    var_copied = copy(var)
    fields = fieldnames(ClimaAnalysis.OutputVar)
    for field in fields
        @test getfield(var, field) == getfield(var_copied, field)
        @test getfield(var, field) !== getfield(var_copied, field)
    end

    # Test reduction
    lat_avg = ClimaAnalysis.average_lat(var)
    @test lat_avg.dims == OrderedDict(["lon" => long, "time" => time])
    @test lat_avg.dim_attributes ==
          OrderedDict(["lon" => Dict("b" => 2), "time" => Dict()])
    @test lat_avg.data == dropdims(mean(data, dims = 3), dims = 3)

    wei_lat_avg = ClimaAnalysis.weighted_average_lat(var)
    @test wei_lat_avg.dims == OrderedDict(["lon" => long, "time" => time])
    @test wei_lat_avg.dim_attributes ==
          OrderedDict(["lon" => Dict("b" => 2), "time" => Dict()])
    weights = ones(size(data))
    for i in eachindex(time)
        for j in eachindex(long)
            for k in eachindex(lat)
                weights[i, j, k] = cosd(lat[k])
            end
        end
    end
    weights ./= mean(cosd.(lat))
    expected_avg = dropdims(mean(data .* weights, dims = 3), dims = 3)
    @test wei_lat_avg.data == expected_avg

    wrong_dims = OrderedDict(["lat" => [0.0, 0.1]])
    wrong_dim_attributes = OrderedDict(["lat" => Dict("a" => 1)])
    wrong_var = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        wrong_dims,
        wrong_dim_attributes,
        [0.0, 0.1],
    )
    @test_logs (
        :warn,
        "Detected latitudes are small. If units are radians, results will be wrong",
    )

    lat_lon_avg = ClimaAnalysis.average_lon(lat_avg)
    @test lat_lon_avg.dims == OrderedDict(["time" => time])
    @test lat_lon_avg.dim_attributes == OrderedDict(["time" => Dict()])

    @test lat_lon_avg.data == dropdims(mean(lat_avg.data, dims = 2), dims = 2)

    lat_lon_time_avg = ClimaAnalysis.average_time(lat_lon_avg)
    @test lat_lon_time_avg.dims == OrderedDict()
    @test lat_lon_time_avg.dim_attributes == OrderedDict()

    @test lat_lon_time_avg.data[] == mean(data)

    @test lat_lon_time_avg.attributes["long_name"] ==
          "hi averaged over latitudes averaged over longitudes averaged over time"
end

@testset "Reductions (box dims)" begin
    x = 0.0:180.0 |> collect
    y = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data = reshape(1.0:(91 * 181 * 11), (11, 181, 91))

    # Identical test pattern to sphere setup, with `dims` modified.
    dims = OrderedDict(["time" => time, "x" => x, "y" => y])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "x" => Dict("b" => 2),
        "y" => Dict("a" => 1),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    y_avg = ClimaAnalysis.average_y(var)
    @test y_avg.dims == OrderedDict(["x" => x, "time" => time])
    @test y_avg.dim_attributes ==
          OrderedDict(["x" => Dict("b" => 2), "time" => Dict()])
    @test y_avg.data == dropdims(mean(data, dims = 3), dims = 3)

    y_x_avg = ClimaAnalysis.average_x(y_avg)
    xy_avg = ClimaAnalysis.average_xy(var)
    @test y_x_avg.data == xy_avg.data
    @test y_x_avg.dims == OrderedDict(["time" => time])
    @test y_x_avg.dim_attributes == OrderedDict(["time" => Dict()])

    @test y_x_avg.data == dropdims(mean(y_avg.data, dims = 2), dims = 2)

    y_x_time_avg = ClimaAnalysis.average_time(y_x_avg)
    xy_time_avg = ClimaAnalysis.average_time(xy_avg)
    @test y_x_time_avg.dims == OrderedDict()
    @test y_x_time_avg.dim_attributes == OrderedDict()

    @test y_x_time_avg.data[] == mean(data)

    @test y_x_time_avg.attributes["long_name"] ==
          "hi averaged over y averaged over x averaged over time"

    @test xy_time_avg.attributes["long_name"] ==
          "hi averaged horizontally averaged over time"
end

@testset "Slicing" begin
    z = 0.0:20.0 |> collect
    time = 100.0:110.0 |> collect

    data = reshape(1.0:(11 * 21), (11, 21))

    dims = OrderedDict(["time" => time, "z" => z])
    dim_attributes =
        OrderedDict(["time" => Dict("units" => "s"), "z" => Dict("b" => 2)])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

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
end

@testset "Windowing" begin
    z = 0.0:20.0 |> collect
    time = 0.0:10.0 |> collect

    data = reshape(1.0:(11 * 21), (11, 21))

    dims = OrderedDict(["time" => time, "z" => z])
    dim_attributes =
        OrderedDict(["time" => Dict("units" => "s"), "z" => Dict("b" => 2)])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

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
end

@testset "Extracting dimension" begin
    @test ClimaAnalysis.Var.find_dim_name(["a", "b"], ["c", "a"]) == "a"
    @test_throws ErrorException ClimaAnalysis.Var.find_dim_name(
        ["a", "b"],
        ["c", "d"],
    )

    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data = reshape(1.0:(91 * 181 * 11), (11, 181, 91))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    @test ClimaAnalysis.time_name(var) == "time"
    @test ClimaAnalysis.longitude_name(var) == "lon"
    @test ClimaAnalysis.latitude_name(var) == "lat"
    @test ClimaAnalysis.times(var) == time
    @test ClimaAnalysis.latitudes(var) == lat
    @test ClimaAnalysis.longitudes(var) == long
end

@testset "Interpolation" begin
    # 1D interpolation with linear data, should yield correct results
    long = -180.0:180.0 |> collect
    data = copy(long)

    longvar = ClimaAnalysis.OutputVar(Dict("long" => long), data)

    @test longvar.([10.5, 20.5]) == [10.5, 20.5]

    # Test error for data outside of range
    @test_throws BoundsError longvar(200.0)

    # 2D interpolation with linear data, should yield correct results
    time = 100.0:110.0 |> collect
    z = 0.0:20.0 |> collect

    data = reshape(1.0:(11 * 21), (11, 21))
    var2d = ClimaAnalysis.OutputVar(Dict("time" => time, "z" => z), data)
    @test var2d.([[105.0, 10.0], [105.5, 10.5]]) == [116.0, 122]
end
