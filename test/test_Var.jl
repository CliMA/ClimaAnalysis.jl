using Test
import ClimaAnalysis

import Interpolations as Intp
import NaNStatistics: nanmean
import NCDatasets: NCDataset
import OrderedCollections: OrderedDict
import Unitful: @u_str
import Dates

@testset "General" begin
    # Add test for short constructor
    long = -180.0:180.0 |> collect
    data = copy(long)

    longvar = ClimaAnalysis.OutputVar(Dict("long" => long), data)

    @test longvar.dims["long"] == long

    # Unitful
    attribs = Dict("long_name" => "hi", "units" => "m/s")
    dim_attributes = OrderedDict(["long" => Dict("units" => "m")])

    var_with_unitful = ClimaAnalysis.OutputVar(
        attribs,
        Dict("long" => long),
        dim_attributes,
        data,
    )

    @test ClimaAnalysis.units(var_with_unitful) == "m s^-1"
    @test var_with_unitful.attributes["units"] == u"m" / u"s"

    # Unparsable unit
    attribs = Dict("long_name" => "hi", "units" => "bob")
    var_without_unitful = ClimaAnalysis.OutputVar(
        attribs,
        Dict("long" => long),
        dim_attributes,
        data,
    )

    @test ClimaAnalysis.units(var_without_unitful) == "bob"
    @test var_without_unitful.attributes["units"] == "bob"

    # Reading directly from file
    ncpath = joinpath(@__DIR__, "topo_drag.res.nc")
    file_var = ClimaAnalysis.OutputVar(ncpath, "t11")
    NCDataset(ncpath) do nc
        @test nc["t11"][:, :, :] == file_var.data
    end

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

@testset "Interpolant boundary conditions" begin
    # Check boundary condtions for lon (equispaced and span), lat (equispaced and span), and
    # time
    lon = 0.5:1.0:359.5 |> collect
    lat = -89.5:1.0:89.5 |> collect
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
        "time" => Dict("units" => "test_units3"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test var.interpolant.et == (Intp.Periodic(), Intp.Flat(), Intp.Throw())

    # Not equispaced for lon and lat
    lon = 0.5:1.0:359.5 |> collect |> x -> push!(x, 42.0) |> sort
    lat = -89.5:1.0:89.5 |> collect |> x -> push!(x, 42.0) |> sort
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test var.interpolant.et == (Intp.Throw(), Intp.Throw(), Intp.Throw())

    # Does not span entire range for and lat
    lon = 0.5:1.0:350.5 |> collect
    lat = -89.5:1.0:80.5 |> collect
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test var.interpolant.et == (Intp.Throw(), Intp.Throw(), Intp.Throw())

    # Dates for the time dimension
    lon = 0.5:1.0:359.5 |> collect
    lat = -89.5:1.0:89.5 |> collect
    time = [
        Dates.DateTime(2020, 3, 1, 1, 1),
        Dates.DateTime(2020, 3, 1, 1, 2),
        Dates.DateTime(2020, 3, 1, 1, 3),
    ]
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
        "time" => Dict("units" => "test_units3"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test isnothing(var.interpolant)
end

@testset "empty" begin
    dims = OrderedDict{String, Vector{Float64}}()
    data = Float64[]
    empty_var = ClimaAnalysis.OutputVar(dims, data)
    @test ClimaAnalysis.isempty(empty_var)

    dims = OrderedDict{String, Vector{Float64}}()
    data = fill(1.0)
    empty_var = ClimaAnalysis.OutputVar(dims, data)
    @test !ClimaAnalysis.isempty(empty_var)

    long = 0.0:180.0 |> collect
    dims = OrderedDict(["long" => long])
    data = ones(size(long))
    dim_attributes = OrderedDict(["lon" => Dict("b" => 2)])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    not_empty_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)
    @test !ClimaAnalysis.isempty(not_empty_var)
end

@testset "Arithmetic operations" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data1 = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict("units" => "s"),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    var1 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data1)

    dim_attributes2 = OrderedDict([
        "time" => Dict("units" => "m"),
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
    @test var1plus10((0.0, 0.0, 0.0)) == var1((0.0, 0.0, 0.0)) + 10

    tenplusvar1 = 10 + var1

    @test tenplusvar1.data == data1 .+ 10
    @test ClimaAnalysis.short_name(tenplusvar1) == "10 + bob"
    @test ClimaAnalysis.long_name(tenplusvar1) == "10 + hi"
    @test tenplusvar1((0.0, 0.0, 0.0)) == 10 + var1((0.0, 0.0, 0.0))

    var1plusvar3 = var1 + var3

    @test var1plusvar3.data == data1 .+ data3
    @test ClimaAnalysis.short_name(var1plusvar3) == "bob + bula"
    @test ClimaAnalysis.long_name(var1plusvar3) == "hi + bob"

    # Test for element wise multiplication and division between OutputVars
    var_times = var1 * var3
    @test var_times.data == var1.data .* var3.data

    var_divide = var1 / var3
    @test var_divide.data == var1.data ./ var3.data
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
    @test lat_avg.data == dropdims(nanmean(data, dims = 3), dims = 3)

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
    weights ./= nanmean(cosd.(lat))
    expected_avg = dropdims(nanmean(data .* weights, dims = 3), dims = 3)
    @test wei_lat_avg.data ≈ expected_avg

    # Test reduction with NaN
    latnan = [1, 2, 3]
    datanan = [10.0, 20.0, NaN]

    dimsnan = OrderedDict(["lat" => latnan])
    dim_attributesnan = OrderedDict(["lat" => Dict("b" => 2)])
    attribsnan = Dict("lat_name" => "hi")
    varnan =
        ClimaAnalysis.OutputVar(attribsnan, dimsnan, dim_attributesnan, datanan)
    @test isnan(ClimaAnalysis.average_lat(varnan; ignore_nan = false).data[])
    @test ClimaAnalysis.average_lat(varnan; weighted = true).data[] ≈
          (datanan[1] * cosd(latnan[1]) + datanan[2] * cosd(latnan[2])) /
          (cosd(latnan[1]) + cosd(latnan[2]))

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

    @test lat_lon_avg.data ==
          dropdims(nanmean(lat_avg.data, dims = 2), dims = 2)

    lat_lon_time_avg = ClimaAnalysis.average_time(lat_lon_avg)
    @test lat_lon_time_avg.dims == OrderedDict()
    @test lat_lon_time_avg.dim_attributes == OrderedDict()

    @test lat_lon_time_avg.data[] == nanmean(data)

    @test lat_lon_time_avg.attributes["long_name"] ==
          "hi averaged over lat (0.0 to 90.0) averaged over lon (0.0 to 180.0) averaged over time (0.0 to 10.0)"
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
    @test y_avg.data == dropdims(nanmean(data, dims = 3), dims = 3)

    y_x_avg = ClimaAnalysis.average_x(y_avg)
    xy_avg = ClimaAnalysis.average_xy(var)
    @test y_x_avg.data == xy_avg.data
    @test y_x_avg.dims == OrderedDict(["time" => time])
    @test y_x_avg.dim_attributes == OrderedDict(["time" => Dict()])

    @test y_x_avg.data == dropdims(nanmean(y_avg.data, dims = 2), dims = 2)

    y_x_time_avg = ClimaAnalysis.average_time(y_x_avg)
    xy_time_avg = ClimaAnalysis.average_time(xy_avg)
    @test y_x_time_avg.dims == OrderedDict()
    @test y_x_time_avg.dim_attributes == OrderedDict()

    @test y_x_time_avg.data[] == nanmean(data)

    @test y_x_time_avg.attributes["long_name"] ==
          "hi averaged over y (0.0 to 90.0) averaged over x (0.0 to 180.0) averaged over time (0.0 to 10.0)"

    @test xy_time_avg.attributes["long_name"] ==
          "hi averaged horizontally over x (0.0 to 180.0) and y (0.0 to 90.0) averaged over time (0.0 to 10.0)"
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
    alt = 0.0:2.0 |> collect

    data = reshape(1.0:(3 * 91 * 181 * 11), (11, 181, 91, 3))

    dims =
        OrderedDict(["time" => time, "lon" => long, "lat" => lat, "z" => alt])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
        "z" => Dict(),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    @test ClimaAnalysis.time_name(var) == "time"
    @test ClimaAnalysis.longitude_name(var) == "lon"
    @test ClimaAnalysis.latitude_name(var) == "lat"
    @test ClimaAnalysis.altitude_name(var) == "z"
    @test ClimaAnalysis.times(var) == time
    @test ClimaAnalysis.latitudes(var) == lat
    @test ClimaAnalysis.longitudes(var) == long
    @test ClimaAnalysis.altitudes(var) == alt
    @test ClimaAnalysis.has_time(var)
    @test ClimaAnalysis.has_longitude(var)
    @test ClimaAnalysis.has_latitude(var)
    @test ClimaAnalysis.has_altitude(var)
    @test ClimaAnalysis.conventional_dim_name("long") == "longitude"
    @test ClimaAnalysis.conventional_dim_name("latitude") == "latitude"
    @test ClimaAnalysis.conventional_dim_name("t") == "time"
    @test ClimaAnalysis.conventional_dim_name("date") == "date"
    @test ClimaAnalysis.conventional_dim_name("z") == "altitude"
    @test ClimaAnalysis.conventional_dim_name("hi") == "hi"
end

@testset "Interpolation" begin
    # 1D interpolation with linear data, should yield correct results
    long = -175.0:175.0 |> collect
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

@testset "Dim of units and range" begin
    x = 0.0:180.0 |> collect
    y = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect
    data = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "x" => x, "y" => y])
    dim_attributes = OrderedDict([
        "time" => Dict("units" => "seconds"),
        "x" => Dict("units" => "km"),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    @test ClimaAnalysis.dim_units(var, "y") == ""
    @test ClimaAnalysis.dim_units(var, "x") == "km"
    @test ClimaAnalysis.range_dim(var, "x") == (0.0, 180.0)
    @test_throws ErrorException(
        "Var does not have dimension z, found [\"time\", \"x\", \"y\"]",
    ) ClimaAnalysis.dim_units(var, "z")
    @test_throws ErrorException(
        "Var does not have dimension z, found [\"time\", \"x\", \"y\"]",
    ) ClimaAnalysis.range_dim(var, "z")
end

@testset "Long name updates" begin
    # Setup to test x_avg, y_avg, xy_avg  
    x = 0.0:180.0 |> collect
    y = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect
    data = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "x" => x, "y" => y])
    dim_attributes = OrderedDict([
        "time" => Dict("units" => "seconds"),
        "x" => Dict("units" => "km"),
        "y" => Dict("units" => "km"),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    y_avg = ClimaAnalysis.average_y(var)
    @test y_avg.attributes["long_name"] == "hi averaged over y (0.0 to 90.0km)"

    x_avg = ClimaAnalysis.average_x(var)
    @test x_avg.attributes["long_name"] == "hi averaged over x (0.0 to 180.0km)"

    xy_avg = ClimaAnalysis.average_xy(var)
    @test xy_avg.attributes["long_name"] ==
          "hi averaged horizontally over x (0.0 to 180.0km) and y (0.0 to 90.0km)"

    # Setup to test average_lat and average_lon 
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data1 = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict("units" => "seconds"),
        "lon" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    attribs = Dict("long_name" => "hi")
    var1 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data1)

    lat_avg = ClimaAnalysis.average_lat(var1)
    lon_avg = ClimaAnalysis.average_lon(var1)
    lat_weighted_avg = ClimaAnalysis.weighted_average_lat(var1)

    @test lon_avg.attributes["long_name"] ==
          "hi averaged over lon (0.0 to 180.0test_units1)"
    @test lat_avg.attributes["long_name"] ==
          "hi averaged over lat (0.0 to 90.0test_units2)"
    @test lat_weighted_avg.attributes["long_name"] ==
          "hi weighted averaged over lat (0.0 to 90.0test_units2)"
end

@testset "Consistent units checking" begin
    x_long = 0.0:180.0 |> collect
    x_lat = 0.0:90.0 |> collect
    x_data = reshape(1.0:(181 * 91), (181, 91))
    x_dims = OrderedDict(["long" => x_long, "lat" => x_lat])
    x_attribs = Dict("long_name" => "hi")
    x_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)

    y_lon = 0.0:90.0 |> collect
    y_lat = 0.0:45.0 |> collect
    y_data = reshape(1.0:(91 * 46), (91, 46))
    y_dims = OrderedDict(["lon" => y_lon, "lat" => y_lat])
    y_attribs = Dict("long_name" => "hi")
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    y_var = ClimaAnalysis.OutputVar(y_attribs, y_dims, y_dim_attribs, y_data)
    @test_nowarn ClimaAnalysis.Var._check_dims_consistent(x_var, y_var)

    # Test if units are consistent between dimensions
    x_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units2"),
        "lat" => Dict("units" => "test_units1"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)
    @test_throws "Units for dimensions [\"long\", \"lat\"] in x is not consistent with units for dimensions [\"lon\", \"lat\"] in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
    )

    # Test if units are missing from any of the dimensions
    x_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units2"),
        "lat" => Dict("units" => ""),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)
    @test_throws "Units for dimensions [\"lat\"] are missing in x and units for dimensions [\"lat\"] are missing in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        x_var,
    )
    @test_throws "Units for dimensions [\"lat\"] are missing in x" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
    )
    @test_throws "Units for dimensions [\"lat\"] are missing in y" ClimaAnalysis.Var._check_dims_consistent(
        y_var,
        x_var,
    )

    # Test if type of dimensions agree
    x_data = reshape(1.0:(91 * 181), (91, 181))
    x_dims = OrderedDict(["lat" => x_lat, "long" => x_long])
    x_dim_attribs = OrderedDict([
        "lat" => Dict("units" => "test_units1"),
        "long" => Dict("units" => "test_units2"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)
    @test_throws "Dimensions do not agree between x ([\"latitude\", \"longitude\"]) and y ([\"longitude\", \"latitude\"])" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
    )

    # Test number of dimensions are the same
    x_data = reshape(1.0:(181), (181))
    x_dims = OrderedDict(["long" => x_long])
    x_attribs = Dict("long_name" => "hi")
    x_dim_attribs = OrderedDict(["long" => Dict("units" => "test_units1")])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)
    @test_throws "Number of dimensions do not match between x (1) and y (2)" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
    )

end

@testset "Reordering" begin
    # Reordering the dimensions of a var to match itself
    src_long = 0.0:180.0 |> collect
    src_lat = 0.0:90.0 |> collect
    src_data = ones(length(src_long), length(src_lat))
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_attribs = Dict("long_name" => "hi")
    src_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )
    reordered_var = ClimaAnalysis.reordered_as(src_var, src_var)
    @test reordered_var.attributes == src_var.attributes
    @test reordered_var.dims == src_var.dims
    @test reordered_var.dim_attributes == src_var.dim_attributes
    @test reordered_var.data == src_var.data

    # Reordering the dimensions of src_var to match a different order of dimensions in
    # dest_var
    dest_long = 20.0:180.0 |> collect
    dest_lat = 30.0:90.0 |> collect
    dest_data = zeros(length(dest_lat), length(dest_long))
    dest_dims = OrderedDict(["lat" => dest_lat, "long" => dest_long])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict([
        "lat" => Dict("units" => "test_units4"),
        "long" => Dict("units" => "test_units3"),
    ])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )
    reordered_var = ClimaAnalysis.reordered_as(src_var, dest_var)
    @test reordered_var.attributes == src_var.attributes
    @test reordered_var.dims ==
          OrderedDict(["lat" => src_lat, "long" => src_long])
    @test reordered_var.dim_attributes == OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => "test_units1"),
    ])
    @test reordered_var.data == ones(length(src_lat), length(src_long))

    # Reordering but dim_attributes is not available for every dimension
    src_dim_attribs_one = OrderedDict(["lat" => Dict("units" => "test_units2")])
    src_dim_attribs_empty = empty(src_dim_attribs_one)
    src_dim_attribs_extra = OrderedDict([
        "extra_info" => "hi",
        "lat" => Dict("units" => "test_units2"),
    ])
    src_var_one = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs_one,
        src_data,
    )
    src_var_empty = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs_empty,
        src_data,
    )
    src_var_extra = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs_extra,
        src_data,
    )
    reordered_var = ClimaAnalysis.reordered_as(src_var_one, dest_var)
    @test reordered_var.dim_attributes == src_dim_attribs_one
    reordered_var = ClimaAnalysis.reordered_as(src_var_empty, dest_var)
    @test reordered_var.dim_attributes == src_dim_attribs_empty
    reordered_var = ClimaAnalysis.reordered_as(src_var_extra, dest_var)
    @test reordered_var.dim_attributes == OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "extra_info" => "hi",
    ])

    # Error checking for dimensions not being the same in both
    src_long = 20.0:180.0 |> collect
    src_lat = 30.0:90.0 |> collect
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_data = ones(length(src_long), length(src_lat))
    dest_lat = 30.0:90.0 |> collect
    dest_dims = OrderedDict(["lat" => dest_lat])
    dest_data = ones(length(src_lat))
    src_var = ClimaAnalysis.OutputVar(src_dims, src_data)
    dest_var = ClimaAnalysis.OutputVar(dest_dims, dest_data)
    @test_throws ErrorException ClimaAnalysis.reordered_as(src_var, dest_var)
end

@testset "Resampling" begin
    src_long = 0.0:180.0 |> collect
    src_lat = 0.0:90.0 |> collect
    src_data = reshape(1.0:(181 * 91), (181, 91))
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_attribs = Dict("long_name" => "hi")
    src_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )

    dest_long = 0.0:90.0 |> collect
    dest_lat = 0.0:45.0 |> collect
    dest_data = reshape(1.0:(91 * 46), (91, 46))
    dest_dims = OrderedDict(["long" => dest_long, "lat" => dest_lat])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )

    @test src_var.data == ClimaAnalysis.resampled_as(src_var, src_var).data
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data == reshape(1.0:(181 * 91), (181, 91))[1:91, 1:46]
    @test_throws BoundsError ClimaAnalysis.resampled_as(dest_var, src_var)

    # BoundsError check
    src_long = 90.0:120.0 |> collect
    src_lat = 45.0:90.0 |> collect
    src_data = zeros(length(src_long), length(src_lat))
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )

    dest_long = 85.0:115.0 |> collect
    dest_lat = 50.0:85.0 |> collect
    dest_data = zeros(length(dest_long), length(dest_lat))
    dest_dims = OrderedDict(["long" => dest_long, "lat" => dest_lat])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )

    @test_throws BoundsError ClimaAnalysis.resampled_as(src_var, dest_var)
end

@testset "Units" begin
    long = -180.0:180.0 |> collect
    data = copy(long)

    # Unitful
    attribs = Dict("long_name" => "hi", "units" => "m/s")
    dim_attributes = OrderedDict(["long" => Dict("units" => "m")])

    var_with_unitful = ClimaAnalysis.OutputVar(
        attribs,
        Dict("long" => long),
        dim_attributes,
        data,
    )
    var_without_unitful = ClimaAnalysis.OutputVar(
        Dict{String, Any}(),
        Dict("long" => long),
        dim_attributes,
        data,
    )

    var_empty_unit = ClimaAnalysis.OutputVar(
        Dict{String, Any}("units" => ""),
        Dict("long" => long),
        dim_attributes,
        data,
    )

    @test ClimaAnalysis.has_units(var_with_unitful)

    # Convert to cm/s
    var_unitful_in_cms = ClimaAnalysis.convert_units(var_with_unitful, "cm/s")

    @test var_unitful_in_cms.data == 100 .* var_with_unitful.data

    # Unparsable because of new units
    @test_throws ErrorException ClimaAnalysis.convert_units(
        var_with_unitful,
        "bob",
    )

    # New units, using conversion function
    var_notunitful = ClimaAnalysis.convert_units(
        var_with_unitful,
        "bob",
        conversion_function = (data) -> 2 * data,
    )

    @test var_notunitful.data == 2 .* var_with_unitful.data

    # New units parsaeble, but with conversion function
    @test_logs (:warn, "Ignoring conversion_function, units are parseable.") ClimaAnalysis.convert_units(
        var_with_unitful,
        "cm/s",
        conversion_function = (data) -> 2 * data,
    )

end

@testset "Integrating on lat and lon" begin
    # Tests for integrate_lon
    lon = collect(range(-179.5, 179.5, 360))
    lon_data = ones(length(lon))
    lon_dims = OrderedDict(["lon" => lon])
    lon_attribs = Dict("long_name" => "hi")
    lon_dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(
        lon_attribs,
        lon_dims,
        lon_dim_attribs,
        lon_data,
    )
    var_integrated_lon = ClimaAnalysis.Var.integrate_lon(var)

    @test isapprox(var_integrated_lon.data[1], 2.0 * π, atol = 0.01)
    @test var_integrated_lon.dims == OrderedDict()
    @test var_integrated_lon.dim_attributes == OrderedDict()
    @test "hi integrated over lon (-179.5 to 179.5deg)" ==
          var_integrated_lon.attributes["long_name"]
    @test_throws "var does not has latitude as a dimension" ClimaAnalysis.Var.integrate_lat(
        var,
    )

    # Tests for integrate_lat
    lat = collect(range(-89.5, 89.5, 180))
    lat_data = ones(length(lat))
    lat_dims = OrderedDict(["lat" => lat])
    lat_attribs = Dict("long_name" => "hi")
    lat_dim_attribs = OrderedDict(["lat" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(
        lat_attribs,
        lat_dims,
        lat_dim_attribs,
        lat_data,
    )
    var_integrated_lat = ClimaAnalysis.Var.integrate_lat(var)

    @test isapprox(var_integrated_lat.data[1], 2.0, atol = 0.01)
    @test var_integrated_lat.dims == OrderedDict()
    @test var_integrated_lat.dim_attributes == OrderedDict()
    @test "hi integrated over lat (-89.5 to 89.5deg)" ==
          var_integrated_lat.attributes["long_name"]
    @test_throws "var does not has longitude as a dimension" ClimaAnalysis.Var.integrate_lon(
        var,
    )

    # Unit checking
    dim_attribs_no_units = OrderedDict([
        "lon" => Dict("units" => ""),
        "lat" => Dict("units" => ""),
    ])
    var_lon_no_units = ClimaAnalysis.OutputVar(
        lon_attribs,
        lon_dims,
        dim_attribs_no_units,
        lon_data,
    )
    @test_throws ErrorException ClimaAnalysis.Var.integrate_lon(
        var_lon_no_units,
    )

    var_lat_no_units = ClimaAnalysis.OutputVar(
        lat_attribs,
        lat_dims,
        dim_attribs_no_units,
        lat_data,
    )
    @test_throws ErrorException ClimaAnalysis.Var.integrate_lat(
        var_lat_no_units,
    )
end

@testset "Integrating on sphere" begin
    # Integrate out all dimensions (lat and lon) from OutputVar
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    integrated_var = ClimaAnalysis.Var.integrate_lonlat(var)
    @test isapprox(integrated_var.data[1], 4 * π, atol = 0.1)
    @test integrated_var.dims == OrderedDict()
    @test integrated_var.dim_attributes == OrderedDict()

    # Integrating out lon and lat to get time series data
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    time = collect(range(0.0, 10.0, 10))
    data = ones(length(lat), length(time), length(lon))
    dims = OrderedDict(["lat" => lat, "time" => time, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "days"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    integrated_var = ClimaAnalysis.Var.integrate_lonlat(var)

    @test all(
        isapprox.(integrated_var.data, [4 * π for _ in 1:10], atol = 0.01),
    )
    @test integrated_var.dims == OrderedDict(["time" => time])
    @test integrated_var.dim_attributes ==
          OrderedDict(["time" => Dict("units" => "days")])
    @test "hi integrated over lon (-179.5 to 179.5deg) and integrated over lat (-89.5 to 89.5deg)" ==
          integrated_var.attributes["long_name"]

    # Unit checking
    dim_attribs_no_lon = OrderedDict([
        "time" => Dict("units" => "days"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs_no_lon, data)
    @test_throws "The unit for lon is missing or is not degree" ClimaAnalysis.Var.integrate_lonlat(
        var,
    )

    dim_attribs_no_lat = OrderedDict([
        "time" => Dict("units" => "days"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs_no_lat, data)
    @test_throws "The unit for lat is missing or is not degree" ClimaAnalysis.Var.integrate_lonlat(
        var,
    )
end

@testset "split_by_season" begin
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    time = [0.0]
    push!(time, 5_184_000.0) # correspond to 2024-3-1
    push!(time, 5_184_001.0)
    push!(time, 13_132_800.0) # correspond to 2024-6-1
    push!(time, 13_132_802.0)
    push!(time, 13_132_803.0)
    data = ones(length(lat), length(time), length(lon))
    dims = OrderedDict(["lat" => lat, "time" => time, "lon" => lon])
    attribs = Dict("long_name" => "hi", "start_date" => "2024-1-1")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    MAM, JJA, SON, DJF = ClimaAnalysis.split_by_season(var)

    # Check size of data
    @test size(MAM.data) == (length(lat), 2, length(lon))
    @test size(JJA.data) == (length(lat), 3, length(lon))
    @test size(SON.data) == (0,)
    @test size(DJF.data) == (length(lat), 1, length(lon))

    # Check times are correct in OutputVars
    @test MAM.dims["time"] == [5_184_000.0, 5_184_001.0]
    @test JJA.dims["time"] == [13_132_800.0, 13_132_802.0, 13_132_803.0]
    @test DJF.dims["time"] == [0.0]

    # Check start date
    MAM.attributes["start_date"] == "2024-1-1"
    JJA.attributes["start_date"] == "2024-1-1"
    DJF.attributes["start_date"] == "2024-1-1"

    # Check empty OutputVar
    @test isempty(SON)

    # Check error handling
    attribs_no_start_date = Dict("long_name" => "hi")
    var =
        ClimaAnalysis.OutputVar(attribs_no_start_date, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.split_by_season(var)

    dim_attribs_no_sec = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "min"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs_no_sec, data)
    @test_throws ErrorException ClimaAnalysis.split_by_season(var)

    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lon))
    dims = OrderedDict(["lon" => lon])
    attribs = Dict("long_name" => "hi", "start_date" => "2024-1-1")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.split_by_season(var)
end

@testset "Compute bias" begin
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_ones = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_ones)

    data_twos = ones(length(lon), length(lat)) .* 2.0
    var_twos = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_twos)

    bias_var = ClimaAnalysis.bias(var_ones, var_twos)
    global_bias = ClimaAnalysis.global_bias(var_ones, var_twos)

    # Check global bias computation, short_name, long_name, and units
    @test bias_var.attributes["global_bias"] == -1.0
    @test global_bias == -1.0
    @test bias_var.data == ones(length(lon), length(lat)) * -1.0
    @test ClimaAnalysis.short_name(bias_var) == "sim-obs_short"
    @test ClimaAnalysis.long_name(bias_var) == "SIM - OBS short"
    @test ClimaAnalysis.units(bias_var) == "kg"

    # Flip order in bias and check computations
    bias_var = ClimaAnalysis.bias(var_twos, var_ones)
    global_bias = ClimaAnalysis.global_bias(var_twos, var_ones)
    @test bias_var.attributes["global_bias"] == 1.0
    @test global_bias == 1.0
    @test bias_var.data == ones(length(lon), length(lat)) * 1.0
end

@testset "Compute mse" begin
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_ones = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_ones)

    data_threes = ones(length(lon), length(lat)) .* 3.0
    var_threes =
        ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_threes)

    squared_error_var = ClimaAnalysis.squared_error(var_ones, var_threes)
    global_mse = ClimaAnalysis.global_mse(var_ones, var_threes)
    global_rmse = ClimaAnalysis.global_rmse(var_ones, var_threes)

    # Check global bias computation, short_name, long_name, and units
    @test squared_error_var.attributes["global_mse"] == (1.0 - 3.0)^2
    @test squared_error_var.attributes["global_rmse"] == 2.0
    @test global_mse == (1.0 - 3.0)^2
    @test global_rmse == 2.0
    @test squared_error_var.data == (data_ones - data_threes) .^ 2
    @test ClimaAnalysis.short_name(squared_error_var) == "(sim-obs)^2_short"
    @test ClimaAnalysis.long_name(squared_error_var) == "(SIM - OBS)^2 short"
    @test ClimaAnalysis.units(squared_error_var) == "kg^2"

    # Flip order in squared_error and check computations
    squared_error_var = ClimaAnalysis.squared_error(var_threes, var_ones)
    global_mse = ClimaAnalysis.global_mse(var_threes, var_ones)
    global_rmse = ClimaAnalysis.global_rmse(var_threes, var_ones)
    @test squared_error_var.attributes["global_mse"] == (3.0 - 1.0)^2
    @test squared_error_var.attributes["global_rmse"] == 2.0
    @test global_mse == (3.0 - 1.0)^2
    @test global_rmse == 2.0
    @test squared_error_var.data == (data_threes - data_ones) .^ 2

    # Check unit handling
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs_unitful =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg^2/m")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_unitful =
        ClimaAnalysis.OutputVar(attribs_unitful, dims, dim_attribs, data_ones)

    attribs_not_unitful = Dict(
        "long_name" => "idk",
        "short_name" => "short",
        "units" => "wacky/weird^2",
    )
    var_not_unitful = ClimaAnalysis.OutputVar(
        attribs_not_unitful,
        dims,
        dim_attribs,
        data_ones,
    )

    var_unitful = ClimaAnalysis.squared_error(var_unitful, var_unitful)
    var_not_unitful =
        ClimaAnalysis.squared_error(var_not_unitful, var_not_unitful)
    @test ClimaAnalysis.units(var_unitful) == "(kg^2 m^-1)^2"
    @test ClimaAnalysis.units(var_not_unitful) == "(wacky/weird^2)^2"
end

@testset "Units and dims check for error functions" begin
    # Missing units for data
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs_missing_data_units =
        Dict("long_name" => "idk", "short_name" => "short")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "rad"),
        "lat" => Dict("units" => "rad"),
    ])
    var_missing_data_units = ClimaAnalysis.OutputVar(
        attribs_missing_data_units,
        dims,
        dim_attribs,
        data_ones,
    )
    @test_throws ErrorException ClimaAnalysis.bias(
        var_missing_data_units,
        var_missing_data_units,
    )
    @test_throws ErrorException ClimaAnalysis.squared_error(
        var_missing_data_units,
        var_missing_data_units,
    )

    # Mismatch units for data
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs_kg =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    attribs_g =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "g")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_kg = ClimaAnalysis.OutputVar(attribs_kg, dims, dim_attribs, data_ones)
    var_g = ClimaAnalysis.OutputVar(attribs_g, dims, dim_attribs, data_ones)
    @test_throws ErrorException ClimaAnalysis.bias(var_kg, var_g)
    @test_throws ErrorException ClimaAnalysis.squared_error(var_kg, var_g)

    # Mismatch units for dims
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs_rad = OrderedDict([
        "lon" => Dict("units" => "rad"),
        "lat" => Dict("units" => "rad"),
    ])
    dim_attribs_deg = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_rad = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs_rad, data_ones)
    var_deg = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs_deg, data_ones)
    @test_throws ErrorException ClimaAnalysis.bias(var_rad, var_deg)
    @test_throws ErrorException ClimaAnalysis.squared_error(var_rad, var_deg)

    # Missing units for dims
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data_ones = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var_missing_dim_units =
        ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_ones)
    @test_throws ErrorException ClimaAnalysis.bias(
        var_missing_dim_units,
        var_missing_dim_units,
    )
    @test_throws ErrorException ClimaAnalysis.squared_error(
        var_missing_dim_units,
        var_missing_dim_units,
    )

    # Missing dims
    lon = collect(range(-179.5, 179.5, 360))
    data_missing_dim = ones(length(lon))
    dims_missing_dim = OrderedDict(["lon" => lon])
    dim_attribs_missing_dim = OrderedDict(["lon" => Dict("units" => "deg")])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var_missing = ClimaAnalysis.OutputVar(
        attribs,
        dims_missing_dim,
        dim_attribs_missing_dim,
        data_missing_dim,
    )
    @test_throws ErrorException ClimaAnalysis.bias(var_missing, var_missing)
    @test_throws ErrorException ClimaAnalysis.squared_error(
        var_missing,
        var_missing,
    )

    # Dimensions should be lon and lat
    lon = collect(range(-179.5, 179.5, 360))
    tal = collect(range(-89.5, 89.5, 180))
    data_tal = ones(length(lon), length(tal))
    dims_tal = OrderedDict(["lon" => lon, "tal" => tal])
    dim_attribs_tal = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "tal" => Dict("units" => "deg"),
    ])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var_tal =
        ClimaAnalysis.OutputVar(attribs, dims_tal, dim_attribs_tal, data_tal)
    @test_throws ErrorException ClimaAnalysis.bias(var_tal, var_tal)
    @test_throws ErrorException ClimaAnalysis.squared_error(var_tal, var_tal)
end

@testset "Setting units" begin
    # Unit exists (unitful)
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_units = ClimaAnalysis.set_units(var, "idk")
    @test ClimaAnalysis.units(var_units) == "idk"

    # Unit exists (not unitful)
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "wacky")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_units = ClimaAnalysis.set_units(var, "idk")
    @test ClimaAnalysis.units(var_units) == "idk"

    # Unit does not exist
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    attribs = Dict("long_name" => "idk", "short_name" => "short")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_units = ClimaAnalysis.set_units(var, "idk")
    @test ClimaAnalysis.units(var_units) == "idk"
end

@testset "Dates to seconds for vars" begin
    # Test for no start date
    time_arr = [
        Dates.DateTime(2020, 3, 1, 1, 1),
        Dates.DateTime(2020, 3, 1, 1, 2),
        Dates.DateTime(2020, 3, 1, 1, 3),
    ]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_s = ClimaAnalysis.Var._dates_to_seconds(var)
    @test ClimaAnalysis.times(var_s) == [0.0, 60.0, 120.0]
    @test var_s.attributes["start_date"] == "2020-03-01T01:01:00"

    # Test for a new start date
    var_s = ClimaAnalysis.Var._dates_to_seconds(
        var;
        new_start_date = "2020-03-01T01:03:00",
    )
    @test ClimaAnalysis.times(var_s) == [-120.0, -60.0, 0.0]
    @test var_s.attributes["start_date"] == "2020-03-01T01:03:00"

    # Test for a new start date as a DateTime object
    var_s = ClimaAnalysis.Var._dates_to_seconds(
        var;
        new_start_date = Dates.DateTime("2020-03-01T01:03:00"),
    )
    @test ClimaAnalysis.times(var_s) == [-120.0, -60.0, 0.0]
    @test var_s.attributes["start_date"] == "2020-03-01T01:03:00"

    # Test for shifting dates
    var_s = ClimaAnalysis.Var._dates_to_seconds(
        var,
        shift_by = t -> t - Dates.Day(15),
    )
    @test ClimaAnalysis.times(var_s) == [0.0, 60.0, 120.0]
    @test var_s.attributes["start_date"] == "2020-02-15T01:01:00"

    # Test for shifting dates and new date together
    var_s = ClimaAnalysis.Var._dates_to_seconds(
        var;
        new_start_date = "2020-03-01T01:00:00",
        shift_by = t -> t + Dates.Minute(4),
    )
    @test ClimaAnalysis.times(var_s) == [300.0, 360.0, 420.0]
    @test var_s.attributes["start_date"] == "2020-03-01T01:00:00"

    # Test constructor for OutputVar that uses _dates_to_seconds
    ncpath = joinpath(@__DIR__, "sample_nc/test_pr.nc")
    file_var = ClimaAnalysis.OutputVar(
        ncpath;
        new_start_date = nothing,
        shift_by = identity,
    )
    @test ClimaAnalysis.times(file_var) == [0.0, 1398902400.0]
    @test file_var.attributes["start_date"] == "1979-01-01T00:00:00"

    # Test constructor for DateTimeNoLeap for the type of the time dimension
    ncpath = joinpath(@__DIR__, "sample_nc/test_gpp.nc")
    file_var = ClimaAnalysis.OutputVar(
        ncpath,
        "gpp";
        new_start_date = nothing,
        shift_by = identity,
    )
    @test ClimaAnalysis.times(file_var) == [0.0]
    @test file_var.attributes["start_date"] == "1980-01-15T00:00:00"

    # Test for error handling
    # Use date dimension instead of time dimension
    date_arr = [
        Dates.DateTime(2020, 3, 1, 1, 1),
        Dates.DateTime(2020, 3, 1, 1, 2),
        Dates.DateTime(2020, 3, 1, 1, 3),
    ]
    data = ones(length(date_arr))
    dims = OrderedDict("date" => date_arr)
    dim_attribs = OrderedDict("date" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.Var._dates_to_seconds(var)

    # Cannot convert if the element type of time array is float
    time_arr = [0.0, 60.0, 120.0]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.Var._dates_to_seconds(var)
end

@testset "End of previous month" begin
    # Shift to beginning of month and shift back one month
    time_arr = [
        Dates.DateTime("2010-02-01T00:00:00"),
        Dates.DateTime("2010-03-01T00:02:00"),
        Dates.DateTime("2010-04-01T00:02:00"),
    ]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_s = ClimaAnalysis.Var._dates_to_seconds(var)
    var_times = ClimaAnalysis.shift_to_start_of_previous_month(var_s)

    @test ClimaAnalysis.times(var_times) == [0.0, 2678400.0, 5097600.0]
    @test var_times.attributes["start_date"] == "2010-01-01T00:00:00"

    # Error checking
    # Dates in time array
    time_arr = [
        Dates.DateTime(2020, 3, 1, 1, 1),
        Dates.DateTime(2020, 3, 1, 1, 2),
        Dates.DateTime(2020, 3, 1, 1, 3),
    ]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.shift_to_start_of_previous_month(
        var,
    )

    # Time is not a dimension
    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.shift_to_start_of_previous_month(
        var,
    )

    # Units is wrong
    time_arr = [
        Dates.DateTime("2010-02-01T00:00:00"),
        Dates.DateTime("2010-03-01T00:02:00"),
        Dates.DateTime("2010-04-01T00:02:00"),
    ]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("units" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_min = ClimaAnalysis.Var._dates_to_seconds(var)
    var_min.dim_attributes["time"]["units"] = "min"
    @test_throws ErrorException ClimaAnalysis.shift_to_start_of_previous_month(
        var_min,
    )

    # Duplicate dates after applying transformations
    time_arr = [
        Dates.DateTime("2010-01-02T00:00:00"),
        Dates.DateTime("2010-01-21T00:00:00"),
        Dates.DateTime("2010-01-31T00:00:00"),
    ]
    data = ones(length(time_arr))
    dims = OrderedDict("time" => time_arr)
    dim_attribs = OrderedDict("time" => Dict("blah" => "blah"))
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_s = ClimaAnalysis.Var._dates_to_seconds(var)
    @test_throws ErrorException ClimaAnalysis.shift_to_start_of_previous_month(
        var_s,
    )
end

@testset "Land and ocean masks" begin
    # Order of dimensions should not matter
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var_latlon = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var_lonlat = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    land_var_lonlat = ClimaAnalysis.apply_landmask(var_lonlat)
    ocean_var_lonlat = ClimaAnalysis.apply_oceanmask(var_lonlat)
    land_var_latlon = ClimaAnalysis.apply_landmask(var_latlon)
    ocean_var_latlon = ClimaAnalysis.apply_oceanmask(var_latlon)
    @test land_var_lonlat.data == land_var_latlon.data'
    @test ocean_var_lonlat.data == ocean_var_latlon.data'

    # Testing with another dimension
    lat = collect(range(-89.5, 89.5, 180))
    times = collect(range(0.0, 100, 2 * 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(times), length(lon))
    dims = OrderedDict(["lat" => lat, "time" => times, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    land_var = ClimaAnalysis.apply_landmask(var) |> ClimaAnalysis.average_time
    ocean_var = ClimaAnalysis.apply_oceanmask(var) |> ClimaAnalysis.average_time
    @test land_var.data |> transpose == land_var_lonlat.data
    @test ocean_var.data |> transpose == ocean_var_lonlat.data

    # Test error handling
    times = collect(range(0.0, 100, 2 * 180))
    data = ones(length(times))
    dims = OrderedDict(["time" => times])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["time" => Dict("units" => "s")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.apply_landmask(var)
    @test_throws ErrorException ClimaAnalysis.apply_oceanmask(var)
end

@testset "Bias and RMSE with masks" begin
    # Test bias and global_bias
    land_var = ClimaAnalysis.OutputVar(ClimaAnalysis.Var.LAND_MASK)
    ocean_var = ClimaAnalysis.OutputVar(ClimaAnalysis.Var.OCEAN_MASK)
    data_zero = zeros(land_var.data |> size)
    zero_var = ClimaAnalysis.OutputVar(
        land_var.attributes,
        land_var.dims,
        land_var.dim_attributes,
        data_zero,
    )

    # Trim data because periodic boundary condition on the edges
    @test ClimaAnalysis.bias(
        land_var,
        zero_var,
        mask = ClimaAnalysis.apply_oceanmask,
    ).data[
        begin:(end - 1),
        :,
    ] == data_zero[begin:(end - 1), :]
    @test ClimaAnalysis.bias(
        ocean_var,
        zero_var,
        mask = ClimaAnalysis.apply_landmask,
    ).data[
        begin:(end - 1),
        :,
    ] == data_zero[begin:(end - 1), :]

    # Not exactly zero because of the periodic boundary condition on the edges
    # which results in some ones in the data
    @test isapprox(
        ClimaAnalysis.global_bias(
            land_var,
            zero_var,
            mask = ClimaAnalysis.apply_oceanmask,
        ),
        0.0,
        atol = 1e-5,
    )
    @test isapprox(
        ClimaAnalysis.global_bias(
            ocean_var,
            zero_var,
            mask = ClimaAnalysis.apply_landmask,
        ),
        0.0,
        atol = 1e-5,
    )

    # Test squared error, global_mse, and global_rmse
    # Trim data because periodic boundary condition on the edges
    @test ClimaAnalysis.squared_error(
        land_var,
        zero_var,
        mask = ClimaAnalysis.apply_oceanmask,
    ).data[
        begin:(end - 1),
        :,
    ] == data_zero[begin:(end - 1), :]
    @test ClimaAnalysis.squared_error(
        ocean_var,
        zero_var,
        mask = ClimaAnalysis.apply_landmask,
    ).data[
        begin:(end - 1),
        :,
    ] == data_zero[begin:(end - 1), :]

    # Not exactly zero because of the periodic boundary condition on the edges
    # which results in some ones in the data
    @test isapprox(
        ClimaAnalysis.global_mse(
            land_var,
            zero_var,
            mask = ClimaAnalysis.apply_oceanmask,
        ),
        0.0,
        atol = 1e-5,
    )
    @test isapprox(
        ClimaAnalysis.global_mse(
            ocean_var,
            zero_var,
            mask = ClimaAnalysis.apply_landmask,
        ),
        0.0,
        atol = 1e-5,
    )

    @test isapprox(
        ClimaAnalysis.global_rmse(
            land_var,
            zero_var,
            mask = ClimaAnalysis.apply_oceanmask,
        ),
        0.0,
        atol = 10^(-2.5),
    )
    @test isapprox(
        ClimaAnalysis.global_rmse(
            ocean_var,
            zero_var,
            mask = ClimaAnalysis.apply_landmask,
        ),
        0.0,
        atol = 10^(-2.5),
    )
end

@testset "Replace" begin
    times = collect(range(0.0, 100, 2 * 180))
    data = ones(length(times))
    data[1:5] .= NaN
    dims = OrderedDict(["time" => times])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["time" => Dict("units" => "s")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_no_nan = ClimaAnalysis.replace(var, NaN => 0.0)
    @test var_no_nan.dims == var.dims
    @test var_no_nan.data == vcat(zeros(5), ones(355))
    @test var_no_nan.attributes == var.attributes
    @test var_no_nan.dim_attributes == var.dim_attributes

    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    data[42:47] .= NaN
    data[32042:32047] .= NaN
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_no_nan = ClimaAnalysis.replace(var, NaN => 1.0)
    @test var_no_nan.dims == var.dims
    @test var_no_nan.data == ones(length(lat), length(lon))
    @test var_no_nan.attributes == var.attributes
    @test var_no_nan.dim_attributes == var.dim_attributes

    lat = collect(range(-89.5, 89.5, 2))
    lon = collect(range(-179.5, 179.5, 2))
    data = [[missing, NaN] [NaN, missing]]
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_no_nan = ClimaAnalysis.replace(var, missing => 1.0, NaN => 2.0)
    @test var_no_nan.dims == var.dims
    @test var_no_nan.data == [[1.0, 2.0] [2.0, 1.0]]
    @test var_no_nan.attributes == var.attributes
    @test var_no_nan.dim_attributes == var.dim_attributes
end

@testset "Set units for dimension" begin
    # Units exist in dim_attribs
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    ClimaAnalysis.set_dim_units!(var, "lat", "degrees")
    @test ClimaAnalysis.dim_units(var, "lat") == "degrees"

    # Units do not exist in dim_attribs as a key
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs =
        OrderedDict(["lat" => Dict(), "lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    ClimaAnalysis.set_dim_units!(var, "lat", "degrees")
    @test ClimaAnalysis.dim_units(var, "lat") == "degrees"

    # Dimension is not present in dim_attribs
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    ClimaAnalysis.set_dim_units!(var, "lat", "degrees")
    @test ClimaAnalysis.dim_units(var, "lat") == "degrees"

    # Error handling
    @test_throws ErrorException ClimaAnalysis.set_dim_units!(
        var,
        "wacky",
        "idk",
    )
end

@testset "Generating masks" begin
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    data[1] = NaN
    data_ones = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    ones_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_ones)

    mask_fn = ClimaAnalysis.make_lonlat_mask(
        var,
        set_to_val = isnan,
        true_val = 0.0,
        false_val = 1.0,
    )
    var1 = mask_fn(ones_var)
    @test var1.data[1] == 0.0

    # Error handling
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lat))
    dims = OrderedDict(["lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.make_lonlat_mask(var)

    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lon))
    dims = OrderedDict(["lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.make_lonlat_mask(var)

    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    t = collect(range(1, 2, 2))
    data = ones(length(lat), length(lon), length(t))
    dims = OrderedDict(["lat" => lat, "lon" => lon, "time" => t])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.make_lonlat_mask(var)
end
