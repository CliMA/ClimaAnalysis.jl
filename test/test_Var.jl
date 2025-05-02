using Test
import ClimaAnalysis

import Statistics
import Interpolations as Intp
import NaNStatistics: nanmean, nansum, nanvar
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

    @test identity(longvar) == longvar

    # Unitful
    attribs = Dict("long_name" => "hi", "units" => "m/s")
    dim_attributes = OrderedDict(["long" => Dict("units" => "m")])

    var_with_unitful = ClimaAnalysis.OutputVar(
        attribs,
        Dict("long" => long),
        dim_attributes,
        data,
    )

    @test ClimaAnalysis.units(var_with_unitful) == "m/s"
    @test var_with_unitful.attributes["units"] == "m/s"

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
        Dict{String, Dict}(),
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
        Dict{String, Dict}(),
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

@testset "Shift longitudes" begin
    lon = collect(range(0.0, 360.0, 361))
    dims = OrderedDict("lon" => lon)
    data = collect(0.0:360.0)
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict("lon" => Dict("units" => "degrees"))
    var_0_360 = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # Shift longitudes from 0 to 360 degrees to -180 to 180 degrees
    var_neg_180_180 = ClimaAnalysis.shift_longitude(var_0_360, -180.0, 180.0)
    @test var_neg_180_180.dims["lon"] == collect(range(-179, 180, 360))
    @test var_neg_180_180.data == vcat(
        collect(range(181.0, 359.0, 179)),
        collect(range(0.0, 180.0, 181)),
    )

    # Shift longitudes from -180 to 180 degrees to 0 to 360 degrees
    var_0_360 = ClimaAnalysis.shift_longitude(var_neg_180_180, 0.0, 360.0)
    @test var_0_360.dims["lon"] == collect(range(0.0, 359.0, 360))
    @test var_0_360.data == collect(range(0.0, 359.0, 360))

    # Shift longitudes from 20 to 380 degrees to 40 to 400 degrees
    lon = [20.0, 30.0, 40.0, 360.0, 370.0]
    dims = OrderedDict("lon" => lon)
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict("lon" => Dict("units" => "degrees"))
    var_20_380 = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    var_40_400 = ClimaAnalysis.shift_longitude(var_20_380, 40.0, 400.0)
    @test var_40_400.dims["lon"] == [40.0, 360.0, 370.0, 380.0, 390.0]
    @test var_40_400.data == [3.0, 4.0, 5.0, 1.0, 2.0]

    # Center longitudes with 3D OutputVar
    lat = [-90.0, -30.0, 30.0, 90.0]
    lon = [-60.0, -30.0, 0.0, 30.0, 60.0]
    time = [0.0, 1.0, 5.0]
    n_elts = length(lat) * length(lon) * length(time)
    dims = OrderedDict("lat" => lat, "lon" => lon, "time" => time)
    size_of_data = (length(lat), length(lon), length(time))
    data = reshape(collect(1.0:n_elts), size_of_data...)
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(
        "lat" => Dict("units" => "degrees"),
        "lon" => Dict("units" => "degrees"),
        "time" => Dict("units" => "seconds"),
    )
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    var = ClimaAnalysis.shift_longitude(var, 0.0, 360.0)
    @test var.dims["lon"] == [0.0, 30.0, 60.0, 300.0, 330.0]
    @test var.data == data[:, [3, 4, 5, 1, 2], :]

    # Error handling
    @test_throws ErrorException ClimaAnalysis.shift_longitude(var, 0.0, 180.0)
end

@testset "Remake" begin
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lat" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    remake_var = ClimaAnalysis.remake(var)
    remake_var.attributes["test"] = "test1"
    remake_var.dims["z"] = [1.0, 2.0, 3.0]
    remake_var.dims["lat"][1] = 10000.0
    remake_var.dim_attributes["lat"]["cool"] = "cool1"
    remake_var.data[1, 1] = 2.0

    # Check that any mutation to remake_var does not change var
    @test var.attributes == Dict("long_name" => "hi")
    @test var.dims == OrderedDict([
        "lat" => collect(range(-89.5, 89.5, 180)),
        "lon" => collect(range(-179.5, 179.5, 360)),
    ])
    @test var.dim_attributes == OrderedDict(["lat" => Dict("units" => "deg")])
    @test var.data == ones(length(lat), length(lon))

    dims1 = OrderedDict(["lat1" => lat, "lon1" => lon])
    dim_attribs1 = OrderedDict(["lat1" => Dict("units" => "deg")])
    remake_var2 =
        ClimaAnalysis.remake(var, dims = dims1, dim_attributes = dim_attribs1)

    @test remake_var2.attributes == Dict("long_name" => "hi")
    @test remake_var2.dims == OrderedDict([
        "lat1" => collect(range(-89.5, 89.5, 180)),
        "lon1" => collect(range(-179.5, 179.5, 360)),
    ])
    @test remake_var2.dim_attributes ==
          OrderedDict(["lat1" => Dict("units" => "deg")])
    @test remake_var2.data == ones(length(lat), length(lon))
end

@testset "Interpolant boundary conditions and error handling" begin
    # Check boundary condtions for lon (equispaced and span), lat (equispaced and span), and
    # time
    lon = 0.5:1.0:359.5 |> collect
    lat = -89.5:1.0:89.5 |> collect
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test intp.et == (Intp.Periodic(), Intp.Flat(), Intp.Throw())

    # Not equispaced for lon and lat
    lon = 0.5:1.0:359.5 |> collect |> x -> push!(x, 42.0) |> sort
    lat = -89.5:1.0:89.5 |> collect |> x -> push!(x, 42.0) |> sort
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test intp.et == (Intp.Throw(), Intp.Throw(), Intp.Throw())

    # Does not span entire range for and lat
    lon = 0.5:1.0:350.5 |> collect
    lat = -89.5:1.0:80.5 |> collect
    time = 1.0:100 |> collect
    data = ones(length(lon), length(lat), length(time))
    dims = OrderedDict(["lon" => lon, "lat" => lat, "time" => time])
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test intp.et == (Intp.Throw(), Intp.Throw(), Intp.Throw())

    # Lon is exactly 360 degrees
    lon = 0.0:1.0:360.0 |> collect
    data = ones(length(lon))
    dims = OrderedDict(["lon" => lon])
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test intp.et == (Intp.Periodic(),)

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
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test isnothing(intp)

    # 2D dimensions
    arb_dim = reshape(collect(range(-89.5, 89.5, 16)), (4, 4))
    data = collect(1:16)
    dims = OrderedDict(["arb_dim" => arb_dim])
    intp = ClimaAnalysis.Var._make_interpolant(dims, data)
    @test isnothing(intp)
end

@testset "empty" begin
    dims = OrderedDict{String, Vector{Float64}}()
    data = Float64[]
    empty_var = ClimaAnalysis.OutputVar(dims, data)
    @test isempty(empty_var)

    dims = OrderedDict{String, Vector{Float64}}()
    data = fill(1.0)
    empty_var = ClimaAnalysis.OutputVar(dims, data)
    @test !isempty(empty_var)

    long = 0.0:180.0 |> collect
    dims = OrderedDict(["long" => long])
    data = ones(size(long))
    dim_attributes = OrderedDict(["lon" => Dict("b" => 2)])
    attribs = Dict("short_name" => "bob", "long_name" => "hi")
    not_empty_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)
    @test !isempty(not_empty_var)
end

@testset "Arithmetic and Mathematical Operations" begin
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data1 = Float64.(collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91))))

    dims = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict("units" => "s"),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict{String, Any}( # Specify Any for value type
        "short_name" => "bob",
        "long_name" => "hi",
        "start_date" => "2008",
        "cool" => "rad",
    )
    var1 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data1)

    # Incompatible dim attributes
    dim_attributes2 = OrderedDict([
        "time" => Dict("units" => "m"),
        "lon" => Dict("lol" => 2),
        "lat" => Dict("a" => 1),
    ])
    var2 = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes2, data1)

    # Compatible variable with different data and some different attributes
    data3 =
        Float64.(5.0 .+ collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91))))
    attribs3 = Dict{String, Any}(
        "long_name" => "bob", # Note: different from var1's long_name ("hi")
        "short_name" => "bula",
        "start_date" => "2008", # Same as var1
        "cool" => "not rad", # Different from var1
    )
    var3 = ClimaAnalysis.OutputVar(attribs3, dims, dim_attributes, data3) # Use same dim_attributes as var1

    @testset "Compatibility" begin
        @test !ClimaAnalysis.arecompatible(var1, var2) # Different dim_attributes
        @test ClimaAnalysis.arecompatible(var1, var3) # Same dims and dim_attributes
    end

    @testset "Binary Operations" begin
        # Addition (+)
        var1plus10 = var1 + 10.0 # Use Float for consistency
        @test var1plus10.data == data1 .+ 10.0
        @test ClimaAnalysis.short_name(var1plus10) == "bob + 10.0"
        @test ClimaAnalysis.long_name(var1plus10) == "hi + 10.0"
        # @test var1plus10((0.0, 0.0, 0.0)) == var1((0.0, 0.0, 0.0)) + 10.0 # Assuming interpolation works
        @test var1plus10.attributes == Dict(
            "short_name" => "bob + 10.0",
            "long_name" => "hi + 10.0",
            "start_date" => "2008",
        )
        @test var1plus10.dims == var1.dims # Check dims preserved
        @test var1plus10.dim_attributes == var1.dim_attributes # Check dim_attributes preserved


        tenplusvar1 = 10.0 + var1
        @test tenplusvar1.data == 10.0 .+ data1
        @test ClimaAnalysis.short_name(tenplusvar1) == "10.0 + bob"
        @test ClimaAnalysis.long_name(tenplusvar1) == "10.0 + hi"
        @test tenplusvar1((0.0, 0.0, 0.0)) == 10.0 + var1((0.0, 0.0, 0.0))
        @test tenplusvar1.attributes == Dict(
            "short_name" => "10.0 + bob",
            "long_name" => "10.0 + hi",
            "start_date" => "2008",
        )

        var1plusvar3 = var1 + var3
        @test var1plusvar3.data == data1 .+ data3
        @test ClimaAnalysis.short_name(var1plusvar3) == "bob + bula"
        @test ClimaAnalysis.long_name(var1plusvar3) == "hi + bob"
        # Only attributes common AND identical should be kept
        @test var1plusvar3.attributes == Dict(
            "short_name" => "bob + bula",
            "long_name" => "hi + bob",
            "start_date" => "2008", # Only start_date is common and identical
        )
        @test var1plusvar3.dims == var1.dims # Check dims preserved


        # Multiplication (*)
        var_times = var1 * var3
        @test var_times.data == var1.data .* var3.data
        @test ClimaAnalysis.short_name(var_times) == "bob * bula"
        @test ClimaAnalysis.long_name(var_times) == "hi * bob"
        @test var_times.attributes == Dict( # Check common attribute logic
            "short_name" => "bob * bula",
            "long_name" => "hi * bob",
            "start_date" => "2008",
        )

        var_times_real = var1 * 2.5
        @test var_times_real.data == var1.data .* 2.5
        @test ClimaAnalysis.short_name(var_times_real) == "bob * 2.5"

        real_times_var = 0.5 * var1
        @test real_times_var.data == 0.5 .* var1.data
        @test ClimaAnalysis.short_name(real_times_var) == "0.5 * bob"

        # Division (/)
        var_divide = var1 / var3
        @test var_divide.data == var1.data ./ var3.data
        @test ClimaAnalysis.short_name(var_divide) == "bob / bula"
        @test ClimaAnalysis.long_name(var_divide) == "hi / bob"
        @test var_divide.attributes == Dict( # Check common attribute logic
            "short_name" => "bob / bula",
            "long_name" => "hi / bob",
            "start_date" => "2008",
        )

        # Subtraction (-)
        var_minus = var1 - var3
        @test var_minus.data == var1.data .- var3.data
        @test ClimaAnalysis.short_name(var_minus) == "bob - bula"
        @test ClimaAnalysis.long_name(var_minus) == "hi - bob"
        @test var_minus.attributes == Dict( # Check common attribute logic
            "short_name" => "bob - bula",
            "long_name" => "hi - bob",
            "start_date" => "2008",
        )

        # Maximum (max)
        var1maxvar3 = max(var1, var3)
        @test var1maxvar3.data == max.(var1.data, var3.data) # Element-wise max
        @test ClimaAnalysis.short_name(var1maxvar3) == "bob max bula"
        @test ClimaAnalysis.long_name(var1maxvar3) == "hi max bob"
        @test var1maxvar3.attributes == Dict( # Common attribute logic
            "short_name" => "bob max bula",
            "long_name" => "hi max bob",
            "start_date" => "2008",
        )

        var1max1000 = max(var1, 1000.0)
        @test var1max1000.data == max.(var1.data, 1000.0)
        @test ClimaAnalysis.short_name(var1max1000) == "bob max 1000.0"

        thousandmaxvar1 = max(1000.0, var1)
        @test thousandmaxvar1.data == max.(1000.0, var1.data)
        @test ClimaAnalysis.short_name(thousandmaxvar1) == "1000.0 max bob"

        # Minimum (min)
        var1minvar3 = min(var1, var3)
        @test var1minvar3.data == min.(var1.data, var3.data) # Element-wise min
        @test ClimaAnalysis.short_name(var1minvar3) == "bob min bula"
        @test ClimaAnalysis.long_name(var1minvar3) == "hi min bob"
        @test var1minvar3.attributes == Dict( # Common attribute logic
            "short_name" => "bob min bula",
            "long_name" => "hi min bob",
            "start_date" => "2008",
        )

        var1min1000 = min(var1, 1000.0)
        @test var1min1000.data == min.(var1.data, 1000.0)
        @test ClimaAnalysis.short_name(var1min1000) == "bob min 1000.0"

        thousandminvar1 = min(1000.0, var1)
        @test thousandminvar1.data == min.(1000.0, var1.data)
        @test ClimaAnalysis.short_name(thousandminvar1) == "1000.0 min bob"

    end

    @testset "Unary Operations" begin
        # Logarithm (log)
        log_var1 = log(var1)
        @test log_var1.data == log.(var1.data)
        @test ClimaAnalysis.short_name(log_var1) == "log(bob)"
        @test ClimaAnalysis.long_name(log_var1) == "log(hi)"
        @test log_var1.attributes == Dict(
            "short_name" => "log(bob)",
            "long_name" => "log(hi)",
            "start_date" => "2008",
        )
        @test log_var1.dims == var1.dims
        @test log_var1.dim_attributes == var1.dim_attributes

        # Exponential (exp)
        exp_var1 = exp(var1)
        @test exp_var1.data ≈ exp.(var1.data) # Use ≈ for potential Float inaccuracies
        @test ClimaAnalysis.short_name(exp_var1) == "exp(bob)"
        @test ClimaAnalysis.long_name(exp_var1) == "exp(hi)"

        # Sine (sin)
        sin_var1 = sin(var1)
        @test sin_var1.data == sin.(var1.data)
        @test ClimaAnalysis.short_name(sin_var1) == "sin(bob)"
        @test ClimaAnalysis.long_name(sin_var1) == "sin(hi)"

        # Cosine (cos)
        cos_var1 = cos(var1)
        @test cos_var1.data == cos.(var1.data)
        @test ClimaAnalysis.short_name(cos_var1) == "cos(bob)"
        @test ClimaAnalysis.long_name(cos_var1) == "cos(hi)"

        # Tangent (tan)
        tan_var1 = tan(var1)
        @test tan_var1.data == tan.(var1.data)
        @test ClimaAnalysis.short_name(tan_var1) == "tan(bob)"
        @test ClimaAnalysis.long_name(tan_var1) == "tan(hi)"

        # Square Root (sqrt)
        sqrt_var1 = sqrt(var1)
        @test sqrt_var1.data == sqrt.(var1.data)
        @test ClimaAnalysis.short_name(sqrt_var1) == "sqrt(bob)"
        @test ClimaAnalysis.long_name(sqrt_var1) == "sqrt(hi)"

        # Unary Minus (-)
        neg_var1 = -var1
        @test neg_var1.data == -(var1.data)
        @test ClimaAnalysis.short_name(neg_var1) == "-(bob)" # Macro uses string(op)
        @test ClimaAnalysis.long_name(neg_var1) == "-(hi)"
        @test neg_var1.attributes == Dict(
            "short_name" => "-(bob)",
            "long_name" => "-(hi)",
            "start_date" => "2008",
        )
        @test neg_var1.dims == var1.dims
    end
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

    # Test weighted average lat with Float32 and NaN
    # We are testing that no error is thrown and do not care about correctness
    # here
    f32_long = Float32.(0.0:180.0 |> collect)
    f32_lat = Float32.(0.0:90.0 |> collect)
    f32_time = Float32.(0.0:10.0 |> collect)
    data_with_NaN = collect(reshape(1.0:(91 * 181 * 11), (11, 181, 91)))
    data_with_NaN[1, 1, 1] = NaN
    dims =
        OrderedDict(["time" => f32_time, "lon" => f32_long, "lat" => f32_lat])
    wei_lat_avg =
        ClimaAnalysis.remake(var, dims = dims, data = data_with_NaN) |>
        ClimaAnalysis.weighted_average_lat

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

    x = [0.0, 1.0]
    y = [0.0, 1.0]
    data = [[NaN, 1.0] [2.0, 3.0]]
    dims = OrderedDict(["x" => x, "y" => y])
    dim_attributes = OrderedDict(["x" => Dict("b" => 2), "y" => Dict("a" => 1)])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)
    avg_var = ClimaAnalysis.average_xy(var)
    @test avg_var.data[] == 2.0
    avg_var = ClimaAnalysis.average_xy(var, ignore_nan = false)
    @test isnan(avg_var.data[])
end

@testset "Average over arbitrary dims" begin
    time = [0.0, 1.0, 2.0]
    lon = [0.0, 10.0, 20.0]
    lat = [0.0, 20.0, 40.0]
    data = collect(reshape(1.0:27.0, (3, 3, 3)))
    data[3, 3, 3] = NaN # replace 27
    data[2, 2, 2] = NaN # replace 14
    data[1, 1, 1] = NaN # replace 1
    dims = OrderedDict(["time" => time, "lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "time" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # ignore_nan = true, update_long_name = true
    avg_var = ClimaAnalysis.Var._average_dims(
        var,
        ("time", "lon", "lat"),
        ignore_nan = true,
        update_long_name = true,
    )
    @test avg_var.data[] ≈ (sum(i for i in 1:27) - 1.0 - 14.0 - 27.0) / 24.0
    @test avg_var.attributes["long_name"] ==
          "hi averaged over time (0.0 to 2.0deg), lon (0.0 to 20.0deg), and lat (0.0 to 40.0deg)"

    # ignore_nan = false, update_long_name = true
    avg_var = ClimaAnalysis.Var._average_dims(
        var,
        ("time", "lon", "lat"),
        ignore_nan = false,
        update_long_name = true,
    )
    @test isnan(avg_var.data[])

    # ignore_nan = false, update_long_name = false
    avg_var = ClimaAnalysis.Var._average_dims(
        var,
        ("time", "lon", "lat"),
        ignore_nan = false,
        update_long_name = false,
    )
    @test avg_var.attributes["long_name"] == "hi"
end

@testset "Average over lon and lat" begin
    # no nan var
    lon = [-20.0, -10.0, 0.0, 10.0, 20.0]
    lat = [-30.0, 15.0, 40.0]
    data = reshape(1.0:(length(lon) * length(lat)), (length(lon), length(lat)))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # nan var
    nan_data = collect(data)
    nan_data[1, 2] = NaN # replace 6
    nan_data[2, 3] = NaN # replace 12
    nan_data[4, 1] = NaN # replace 4
    nan_var = ClimaAnalysis.remake(var, data = nan_data)

    # 2D case
    # ignore_nan = true, weighted = false
    avg_var =
        ClimaAnalysis.average_lonlat(var; ignore_nan = true, weighted = false)
    @test avg_var.data[] ≈ sum(i for i in 1:15) / 15
    @test avg_var.attributes["long_name"] ==
          "hi averaged over lat (-30.0 to 40.0deg) and lon (-20.0 to 20.0deg)"
    @test avg_var.dim_attributes == OrderedDict{String, Dict{String, String}}()
    @test avg_var.dims == OrderedDict{String, Vector{Float64}}()

    avg_nan_var = ClimaAnalysis.average_lonlat(
        nan_var;
        ignore_nan = true,
        weighted = false,
    )
    @test avg_nan_var.data[] ≈ (sum(i for i in 1:15) - 6.0 - 12.0 - 4.0) / 12

    # ignore_nan = false, weighted = false
    avg_var =
        ClimaAnalysis.average_lonlat(var; ignore_nan = false, weighted = false)
    @test avg_var.data[] ≈ sum(i for i in 1:15) / 15

    avg_nan_var = ClimaAnalysis.average_lonlat(
        nan_var;
        ignore_nan = false,
        weighted = false,
    )
    @test isnan(avg_nan_var.data[])

    # ignore_nan = true, weighted = true
    avg_var = ClimaAnalysis.weighted_average_lonlat(var; ignore_nan = true)
    @test avg_var.attributes["long_name"] ==
          "hi weighted averaged over lat (-30.0 to 40.0deg) and lon (-20.0 to 20.0deg)"
    @test avg_var.data[] ≈
          sum(data .* reshape(cosd.(lat), (1, 3)) ./ sum(5.0 * cosd.(lat)))

    avg_nan_var =
        ClimaAnalysis.weighted_average_lonlat(nan_var; ignore_nan = true)
    @test avg_nan_var.data[] ≈ nansum(
        nan_data .* reshape(cosd.(lat), (1, 3)) ./
        (sum(5.0 * cosd.(lat)) - cosd(15) - cosd(40) - cosd(-30)),
    )

    # ignore_nan = false, weighted = true
    avg_var = ClimaAnalysis.weighted_average_lonlat(var; ignore_nan = false)
    @test avg_var.data[] ≈
          sum(data .* reshape(cosd.(lat), (1, 3)) ./ sum(5.0 * cosd.(lat)))

    avg_nan_var =
        ClimaAnalysis.weighted_average_lonlat(nan_var; ignore_nan = false)
    @test isnan(avg_nan_var.data[])

    # 3D case
    lat = [-30.0, 15.0, 40.0]
    time = [0.0, 42.0, 85.0]
    lon = [-20.0, -10.0, 0.0, 10.0, 20.0]
    data0 =
        reshape(reshape(1.0:(length(lon) * length(lat)), (5, 3))', (3, 1, 5))
    data1 = collect(
        reshape(reshape(1.0:(length(lon) * length(lat)), (5, 3))', (3, 1, 5)),
    )
    data1[2, 1, 1] = NaN # replace 6
    data1[3, 1, 2] = NaN # replace 12
    data1[1, 1, 4] = NaN # replace 4
    data2 = reshape(
        Float64[[42, 47, 7] [3, 8, 4] [5, 3, 2] [9, 11, 12] [15, 14, 13]],
        (3, 1, 5),
    )
    data3d = cat(data0, data1, data2, dims = 2)
    dims = OrderedDict(["lat" => lat, "time" => time, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
        "lon" => Dict("units" => "deg"),
    ])
    var_3d = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data3d)

    # Precomputed values with ignore_nan = true
    # data0
    avg0 = sum(i for i in 1:15) / 15
    avg0_weighted =
        sum(data .* reshape(cosd.(lat), (1, 3)) ./ sum(5.0 * cosd.(lat)))

    # data1
    avg1 = (sum(i for i in 1:15) - 6.0 - 12.0 - 4.0) / 12
    avg1_weighted = nansum(
        nan_data .* reshape(cosd.(lat), (1, 3)) ./
        (sum(5.0 * cosd.(lat)) - cosd(15) - cosd(40) - cosd(-30)),
    )

    # data2
    avg2 = sum(data2) / 15
    avg2_weighted =
        sum(data2 .* reshape(cosd.(lat), (3, 1, 1)) ./ sum(5.0 * cosd.(lat)))

    # ignore_nan = true, weighted = false
    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_3d;
        ignore_nan = true,
        weighted = false,
    )
    @test isapprox(avg_var_3d.data, [avg0, avg1, avg2])
    @test avg_var_3d.attributes["long_name"] ==
          "hi averaged over lat (-30.0 to 40.0deg) and lon (-20.0 to 20.0deg)"
    @test avg_var_3d.dim_attributes ==
          OrderedDict(["time" => Dict("units" => "s")])
    @test avg_var_3d.dims == OrderedDict(["time" => time])

    # ignore_nan = false, weighted = false
    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_3d;
        ignore_nan = false,
        weighted = false,
    )
    @test isapprox(avg_var_3d.data, [avg0, NaN, avg2], nans = true)

    # ignore_nan = true, weighted = true
    avg_var_3d =
        ClimaAnalysis.average_lonlat(var_3d; ignore_nan = true, weighted = true)
    @test isapprox(
        avg_var_3d.data,
        [avg0_weighted, avg1_weighted, avg2_weighted],
    )

    # ignore_nan = false, weighted = true
    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_3d;
        ignore_nan = false,
        weighted = true,
    )
    @test isapprox(
        avg_var_3d.data,
        [avg0_weighted, NaN, avg2_weighted],
        nans = true,
    )

    # Test order of time, lon, lat
    var_reordered_3d = permutedims(var_3d, ("time", "lon", "lat"))

    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_reordered_3d,
        ignore_nan = true,
        weighted = false,
    )
    @test isapprox(avg_var_3d.data, [avg0, avg1, avg2])

    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_reordered_3d,
        ignore_nan = false,
        weighted = false,
    )
    @test isapprox(avg_var_3d.data, [avg0, NaN, avg2], nans = true)

    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_reordered_3d,
        ignore_nan = true,
        weighted = true,
    )
    @test isapprox(
        avg_var_3d.data,
        [avg0_weighted, avg1_weighted, avg2_weighted],
    )

    avg_var_3d = ClimaAnalysis.average_lonlat(
        var_reordered_3d,
        ignore_nan = false,
        weighted = true,
    )
    @test isapprox(
        avg_var_3d.data,
        [avg0_weighted, NaN, avg2_weighted],
        nans = true,
    )


    # Error handling
    lon = [-20.0, -10.0, 0.0, 10.0]
    lat = [-30.0, 1.0]
    data = reshape(1.0:(length(lon) * length(lat)), (length(lon), length(lat)))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_logs (
        :warn,
        "Detected latitudes are small. If units are radians, results will be wrong",
    ) match_mode = :any ClimaAnalysis.weighted_average_lonlat(var)
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
    @test ClimaAnalysis.conventional_dim_name("pfull") == "pressure"
    @test ClimaAnalysis.find_corresponding_dim_name(
        "lat",
        ["latitude", "longitude"],
    ) == "latitude"
    @test ClimaAnalysis.find_corresponding_dim_name(
        "lat",
        ["longitude", "latitude"],
    ) == "latitude"
    @test ClimaAnalysis.find_corresponding_dim_name(
        "hi",
        ["longitude", "hi"],
    ) == "hi"
    @test_throws ErrorException ClimaAnalysis.find_corresponding_dim_name(
        "lat",
        ["longitude", "pfull"],
    )

    # Pressure dim
    pressure = 0:100.0 |> collect
    data = ones(length(pressure))

    attribs = Dict("short_name" => "K")
    dim_attribs = Dict{String, Dict}()
    pressure_var = ClimaAnalysis.OutputVar(
        attribs,
        Dict("pfull" => pressure),
        dim_attribs,
        data,
    )

    @test ClimaAnalysis.has_pressure(pressure_var)
    @test ClimaAnalysis.pressures(pressure_var) == pressure
    @test ClimaAnalysis.pressure_name(pressure_var) == "pfull"
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
        "x" => Dict("units" => u"km"),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    @test ClimaAnalysis.dim_units(var, "y") == ""
    @test ClimaAnalysis.dim_units(var, "time") == "seconds"
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

    # Test _update_long_name_generic! manually
    # One dimension
    ClimaAnalysis.Var._update_long_name_generic!(var, var, "time", "reduced")
    @test var.attributes["long_name"] ==
          "hi reduced over time (0.0 to 10.0seconds)"

    # Two dimensions
    var.attributes["long_name"] = "hi"
    ClimaAnalysis.Var._update_long_name_generic!(
        var,
        var,
        ("time", "x"),
        "reduced",
    )
    @test var.attributes["long_name"] ==
          "hi reduced over time (0.0 to 10.0seconds) and x (0.0 to 180.0km)"

    # Three dimensions
    var.attributes["long_name"] = "hi"
    ClimaAnalysis.Var._update_long_name_generic!(
        var,
        var,
        ("time", "x", "y"),
        "reduced",
    )
    @test var.attributes["long_name"] ==
          "hi reduced over time (0.0 to 10.0seconds), x (0.0 to 180.0km), and y (0.0 to 90.0km)"
end

@testset "Variance" begin
    time = 0.0:3.0 |> collect
    lon = 0.0:4.0 |> collect
    lat = 0.0:5.0 |> collect

    data = collect(reshape(range(0.0, 1000.0, length = 4 * 5 * 6), (4, 5, 6)))
    data[1, 1, 1] = -10.0
    data[2, 2, 2] = -100.0
    data[3, 3, 3] = -1000.0

    dims = OrderedDict(["time" => time, "lon" => lon, "lat" => lat])
    dim_attributes = OrderedDict([
        "time" => Dict(),
        "lon" => Dict("b" => 2),
        "lat" => Dict("a" => 1),
    ])
    attribs = Dict("long_name" => "hi")
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attributes, data)

    # No NaNs
    var_time = ClimaAnalysis.variance_time(var, ignore_nan = false)
    var_lat = ClimaAnalysis.variance_lat(var, ignore_nan = false)
    var_lon = ClimaAnalysis.variance_lon(var, ignore_nan = false)

    # Note: nanvar and var produces slightly different results, so we use
    # Statistics.var instead of nanvar even though there are NaNs
    @test var_time.data == dropdims(Statistics.var(data, dims = 1), dims = 1)
    @test var_time.attributes["long_name"] ==
          "hi variance over time (0.0 to 3.0)"
    @test var_time.dims == OrderedDict(["lon" => lon, "lat" => lat])

    @test var_lon.data == dropdims(Statistics.var(data, dims = 2), dims = 2)
    @test var_lon.attributes["long_name"] == "hi variance over lon (0.0 to 4.0)"
    @test var_lon.dims == OrderedDict(["time" => time, "lat" => lat])

    @test var_lat.data == dropdims(Statistics.var(data, dims = 3), dims = 3)
    @test var_lat.attributes["long_name"] == "hi variance over lat (0.0 to 5.0)"
    @test var_lat.dims == OrderedDict(["time" => time, "lon" => lon])

    # NaNs
    var.data[1, 1, 1] = NaN
    var.data[2, 2, 2] = NaN
    var.data[3, 3, 3] = NaN

    var_time = ClimaAnalysis.variance_time(var)
    var_lat = ClimaAnalysis.variance_lat(var)
    var_lon = ClimaAnalysis.variance_lon(var)
    @test var_time.data == dropdims(nanvar(data, dims = 1), dims = 1)
    @test var_lon.data == dropdims(nanvar(data, dims = 2), dims = 2)
    @test var_lat.data == dropdims(nanvar(data, dims = 3), dims = 3)

    # Error handling
    # Missing dimension
    @test_throws ErrorException ClimaAnalysis.variance_time(var_time)
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
    x_dims = OrderedDict(["t" => x_lat, "pfull" => x_long])
    x_dim_attribs = OrderedDict([
        "t" => Dict("units" => "test_units1"),
        "pfull" => Dict("units" => "test_units2"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)
    @test_throws "Dimensions do not agree between x (Set([\"time\", \"pressure\"])) and y (Set([\"latitude\", \"longitude\"]))" ClimaAnalysis.Var._check_dims_consistent(
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

    # Check with keyword argument
    x_long = 0.0:180.0 |> collect
    x_lat = 0.0:90.0 |> collect
    x_data = reshape(1.0:(181 * 91), (91, 181))
    x_dims = OrderedDict(["lat" => x_lat, "long" => x_long])
    x_attribs = Dict("long_name" => "hi")
    x_dim_attribs = OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => "test_units1"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)

    y_long = 0.0:180.0 |> collect
    y_pfull = 0.0:2.0 |> collect
    y_time = 0.0:3.0 |> collect
    y_data = ones(length(y_long), length(y_pfull), length(y_time))
    y_dims = OrderedDict(["lon" => y_long, "pfull" => y_pfull, "t" => y_time])
    y_attribs = Dict("long_name" => "hello")
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "pfull" => Dict("units" => "something"),
        "time" => Dict("units" => "idk"),
    ])
    y_var = ClimaAnalysis.OutputVar(y_attribs, y_dims, y_dim_attribs, y_data)
    @test_nowarn ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["longitude"],
    )
    @test_nowarn ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = "longitude",
    )

    # Test if dimension is not present in x or y
    @test_throws "Cannot find space in the dimension names of x ([\"latitude\", \"longitude\"])" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["space"],
    )
    # Test if dimension is present in only one of them
    @test_throws "Cannot find latitude in the dimension names of y ([\"longitude\", \"pressure\", \"time\"])" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = "lat",
    )

    # Test if units are consistent between dimensions
    x_dim_attribs = OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => "this should not"),
    ])
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "be the same"),
        "pfull" => Dict("units" => "something"),
        "time" => Dict("units" => "idk"),
    ])
    x_var = ClimaAnalysis.remake(x_var, dim_attributes = x_dim_attribs)
    y_var = ClimaAnalysis.remake(y_var, dim_attributes = y_dim_attribs)
    @test_throws "Units for dimensions [\"long\"] in x is not consistent with units for dimensions [\"lon\"] in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = "longitude",
    )

    # Test if units are missing from any of the dimensions
    x_dim_attribs = OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => ""),
    ])
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => ""),
        "pfull" => Dict("units" => "something"),
        "time" => Dict("units" => "idk"),
    ])
    x_var = ClimaAnalysis.remake(x_var, dim_attributes = x_dim_attribs)
    y_var = ClimaAnalysis.remake(y_var, dim_attributes = y_dim_attribs)
    @test_throws "Units for dimensions [\"long\"] are missing in x and units for dimensions [\"lon\"] are missing in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = "longitude",
    )

    # Test with more than one dimension supplied for keyword argument
    x_pfull = 0.0:10.0 |> collect
    x_lat = 0.0:90.0 |> collect
    x_long = 0.0:180.0 |> collect
    x_data = reshape(1.0:(11 * 181 * 91), (11, 91, 181))
    x_dims = OrderedDict([
        "pressure_level" => x_pfull,
        "lat" => x_lat,
        "long" => x_long,
    ])
    x_attribs = Dict("long_name" => "hi")
    x_dim_attribs = OrderedDict([
        "pressure_level" => Dict("units" => "something"),
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => "test_units1"),
    ])
    x_var = ClimaAnalysis.OutputVar(x_attribs, x_dims, x_dim_attribs, x_data)

    y_long = 0.0:180.0 |> collect
    y_pfull = 0.0:2.0 |> collect
    y_time = 0.0:3.0 |> collect
    y_data = ones(length(y_long), length(y_pfull), length(y_time))
    y_dims = OrderedDict(["lon" => y_long, "pfull" => y_pfull, "t" => y_time])
    y_attribs = Dict("long_name" => "hello")
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "pfull" => Dict("units" => "something"),
        "time" => Dict("units" => "idk"),
    ])
    y_var = ClimaAnalysis.OutputVar(y_attribs, y_dims, y_dim_attribs, y_data)
    @test_nowarn ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["longitude", "pfull"],
    )

    # Test if dimension is not present in x or y
    @test_throws "Cannot find no in the dimension names of x ([\"pressure\", \"latitude\", \"longitude\"])" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["no", "dim"],
    )

    # Test if dimension is present in only one of them
    @test_throws "Cannot find latitude in the dimension names of y ([\"longitude\", \"pressure\", \"time\"])" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["latitude", "lon"],
    )

    # Test if units are consistent between dimensions
    x_dim_attribs = OrderedDict([
        "pressure_level" => Dict("units" => "not something"),
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => "test_units1"),
    ])
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "pfull" => Dict("units" => "something"),
        "time" => Dict("units" => "idk"),
    ])
    x_var = ClimaAnalysis.remake(x_var, dim_attributes = x_dim_attribs)
    y_var = ClimaAnalysis.remake(y_var, dim_attributes = y_dim_attribs)
    @test_throws "Units for dimensions [\"pressure_level\"] in x is not consistent with units for dimensions [\"pfull\"] in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["longitude", "pfull"],
    )

    # Test if units are missing from any of the dimensions
    x_dim_attribs = OrderedDict([
        "pressure_level" => Dict("units" => "something"),
        "lat" => Dict("units" => "test_units2"),
        "long" => Dict("units" => ""),
    ])
    y_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "pfull" => Dict("units" => ""),
        "time" => Dict("units" => "idk"),
    ])
    x_var = ClimaAnalysis.remake(x_var, dim_attributes = x_dim_attribs)
    y_var = ClimaAnalysis.remake(y_var, dim_attributes = y_dim_attribs)
    @test_throws "Units for dimensions [\"long\"] are missing in x and units for dimensions [\"pfull\"] are missing in y" ClimaAnalysis.Var._check_dims_consistent(
        x_var,
        y_var,
        dim_names = ["longitude", "pfull"],
    )
end

@testset "Permuting dims" begin
    # Reordering the dimensions of a var to match itself
    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    data = ones(length(long), length(lat))
    dims = OrderedDict(["long" => long, "lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # Testing with different types for perm
    perms = [
        ("latitude", "lon"),
        ["lat", "longitude"],
        keys(OrderedDict("latitude" => 1, "long" => 2)),
    ]
    for perm in perms
        permuted_var = permutedims(var, perm)
        @test permuted_var.attributes == var.attributes
        @test permuted_var.dims == OrderedDict(["lat" => lat, "long" => long])
        @test permuted_var.data == permutedims(data, (2, 1))
        @test permuted_var.dim_attributes == OrderedDict([
            "lat" => Dict("units" => "test_units2"),
            "long" => Dict("units" => "test_units1"),
        ])
    end
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
        "extra_info" => Dict("hi" => "bye"),
        "lat" => Dict("units" => "test_units2"),
    ])
    src_var_one =
        ClimaAnalysis.remake(src_var, dim_attributes = src_dim_attribs_one)
    src_var_empty =
        ClimaAnalysis.remake(src_var, dim_attributes = src_dim_attribs_empty)
    src_var_extra =
        ClimaAnalysis.remake(src_var, dim_attributes = src_dim_attribs_extra)
    reordered_var = ClimaAnalysis.reordered_as(src_var_one, dest_var)
    @test reordered_var.dim_attributes == src_dim_attribs_one
    reordered_var = ClimaAnalysis.reordered_as(src_var_empty, dest_var)
    @test reordered_var.dim_attributes == src_dim_attribs_empty
    reordered_var = ClimaAnalysis.reordered_as(src_var_extra, dest_var)
    @test reordered_var.dim_attributes == OrderedDict([
        "lat" => Dict("units" => "test_units2"),
        "extra_info" => Dict("hi" => "bye"),
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

@testset "Resampling over all dimensions" begin
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
    dest_var = ClimaAnalysis.remake(src_var, data = dest_data, dims = dest_dims)

    @test src_var.data == ClimaAnalysis.resampled_as(src_var, src_var).data
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data == reshape(1.0:(181 * 91), (181, 91))[1:91, 1:46]
    @test_throws BoundsError ClimaAnalysis.resampled_as(dest_var, src_var)

    # Test if reordering is automatically done
    src_var_transpose = permutedims(src_var, ("latitude", "longitude"))
    resampled_var = ClimaAnalysis.resampled_as(src_var_transpose, dest_var)
    @test resampled_var.data == reshape(1.0:(181 * 91), (181, 91))[1:91, 1:46]

    # Test with ordered iterable for dims
    resampled_var =
        ClimaAnalysis.resampled_as(src_var, long = dest_long, lat = dest_lat)
    @test resampled_var.data == reshape(1.0:(181 * 91), (181, 91))[1:91, 1:46]
    @test_throws BoundsError ClimaAnalysis.resampled_as(
        dest_var,
        long = src_long,
        lat = src_lat,
    )

    # BoundsError check
    src_long = 90.0:120.0 |> collect
    src_lat = 45.0:90.0 |> collect
    src_data = zeros(length(src_long), length(src_lat))
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_var = ClimaAnalysis.remake(src_var, dims = src_dims, data = src_data)

    dest_long = 85.0:115.0 |> collect
    dest_lat = 50.0:85.0 |> collect
    dest_data = zeros(length(dest_long), length(dest_lat))
    dest_dims = OrderedDict(["long" => dest_long, "lat" => dest_lat])
    dest_var =
        ClimaAnalysis.remake(dest_var, data = dest_data, dims = dest_dims)

    @test_throws BoundsError ClimaAnalysis.resampled_as(src_var, dest_var)

    # Error handling with ordered iterable for dims
    @test_throws ErrorException ClimaAnalysis.resampled_as(
        src_var,
        pfull = [1.0, 2.0],
    )
end

@testset "Resampling ongrid and oncenter" begin
    # Checking oncell (1d)
    src_long = 0.0:359.0 |> collect
    src_data = reshape(collect(1.0:(360.0)), (360,))
    src_dims = OrderedDict(["long" => src_long])
    src_attribs = Dict("long_name" => "hi")
    src_dim_attribs = OrderedDict(["long" => Dict("units" => "test_units1")])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )

    dest_long = [-2.0, -1.5, -1.0, -0.5, 0.0, 1.0, 2.0]
    dest_data = ones(size(dest_long))
    dest_dims = OrderedDict(["long" => dest_long])
    dest_var = ClimaAnalysis.remake(src_var, data = dest_data, dims = dest_dims)
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data == [359.0, 359.5, 360.0, 180.5, 1.0, 2.0, 3.0]

    # Checking ongrid (1d)
    src_long = 0.0:360.0 |> collect
    src_data = reshape(collect(1.0:(361.0)), (361,))
    src_dims = OrderedDict(["long" => src_long])
    src_attribs = Dict("long_name" => "hi")
    src_dim_attribs = OrderedDict(["long" => Dict("units" => "test_units1")])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data == [359.0, 359.5, 360.0, 360.5, 1.0, 2.0, 3.0]

    # Checking oncell (2d)
    src_long = 0.0:359.0 |> collect
    src_lat = 0.0:2.0 |> collect
    src_data = reshape(collect(1.0:(360.0 * 3)), (360, 3))
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
    dest_long = [-2.0, -1.5, -1.0, 1.0]
    dest_lat = [1.0, 1.5]
    dest_data = ones(length(dest_long), length(dest_lat))
    dest_dims = OrderedDict(["long" => dest_long, "lat" => dest_lat])
    dest_var = ClimaAnalysis.remake(src_var, data = dest_data, dims = dest_dims)
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data == [[719, 719.5, 720, 362] [899, 899.5, 900, 542]]

    # Checking ongrid (2d)
    src_long = 0.0:360.0 |> collect
    src_lat = 0.0:2.0 |> collect
    src_data = reshape(collect(1.0:(361.0 * 3)), (361, 3))
    src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
    src_var = ClimaAnalysis.remake(src_var, data = src_data, dims = src_dims)
    resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_var.data ==
          [[720, 720.5, 721, 363] [900.5, 901.0, 901.5, 543.5]]
end

@testset "Resampling with dim_names keyword" begin
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
    dest_var = ClimaAnalysis.remake(src_var, data = dest_data, dims = dest_dims)

    # Resampling over all dimensions with dim_names (should be the same as resampled_as)
    resampled_src_var1 = ClimaAnalysis.resampled_as(
        src_var,
        dest_var,
        dim_names = ["latitude", "longitude"],
    )
    resampled_src_var2 = ClimaAnalysis.resampled_as(src_var, dest_var)
    @test resampled_src_var1.data == resampled_src_var2.data
    @test resampled_src_var1.dims == resampled_src_var2.dims

    # Resampling an OutputVar with small number of dimensions (src_var) on an OutputVar with a large number of dimensions (dest_var)
    dest_time = 0.0:3.0 |> collect
    dest_long = 0.0:42.0 |> collect
    dest_pfull = 0.0:10.0 |> collect
    dest_data = reshape(
        1.0:(length(dest_time) * length(dest_long) * length(dest_pfull)),
        (length(dest_time), length(dest_long), length(dest_pfull)),
    )
    dest_dims = OrderedDict([
        "time" => dest_time,
        "lon" => dest_long,
        "pfull" => dest_pfull,
    ])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict([
        "lon" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
    ])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )
    resampled_src_var =
        ClimaAnalysis.resampled_as(src_var, dest_var, dim_names = "longitude")
    @test resampled_src_var.data == src_var.data[1:43, :]
    @test resampled_src_var.dims["long"] == dest_var.dims["lon"]

    # Resample OutputVar with 2 dimensions (src_var) on an OutputVar
    # with 3 dimensions (dest_var) with two dimensions being
    # resampled and in different order in both OutputVars
    dest_lat = 0.0:3.0 |> collect
    dest_long = 0.0:42.0 |> collect
    dest_pull = 0.0:10.0 |> collect
    dest_data = reshape(
        1.0:(length(dest_lat) * length(dest_long) * length(dest_pfull)),
        (length(dest_lat), length(dest_long), length(dest_pfull)),
    )
    dest_dims = OrderedDict([
        "latitude" => dest_lat,
        "longitude" => dest_long,
        "pfull" => dest_pfull,
    ])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict([
        "latitude" => Dict("units" => "test_units2"),
        "longitude" => Dict("units" => "test_units1"),
    ])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )
    resampled_src_var = ClimaAnalysis.resampled_as(
        src_var,
        dest_var,
        dim_names = ["lon", "latitude"],
    )
    @test resampled_src_var.data == src_var.data[1:43, 1:4]
    @test resampled_src_var.dims["lat"] == dest_var.dims["latitude"]
    @test resampled_src_var.dims["long"] == dest_var.dims["longitude"]

    # Resample OutputVar with 3 dimensions (src_var) on an OutputVar
    # with 2 dimensions (dest_var) with two dimensions being
    # resampled and in different order in both OutputVars
    src_long = 0.0:60.0 |> collect
    src_lat = 0.0:30.0 |> collect
    src_time = 0.0:5.0 |> collect
    src_data = reshape(
        1.0:(length(src_long) * length(src_lat) * length(src_time)),
        (length(src_long), length(src_lat), length(src_time)),
    )
    src_dims =
        OrderedDict(["long" => src_long, "lat" => src_lat, "time" => src_time])
    src_attribs = Dict("long_name" => "hi")
    src_dim_attribs = OrderedDict([
        "long" => Dict("units" => "test_units1"),
        "lat" => Dict("units" => "test_units2"),
        "time" => Dict("units" => "seconds"),
    ])
    src_var = ClimaAnalysis.OutputVar(
        src_attribs,
        src_dims,
        src_dim_attribs,
        src_data,
    )

    dest_t = 0.0:3.0 |> collect
    dest_longitude = 0.0:45.0 |> collect
    dest_data = reshape(
        1.0:(length(dest_t) * length(dest_longitude)),
        (length(dest_t), length(dest_longitude)),
    )
    dest_dims = OrderedDict(["t" => dest_t, "longitude" => dest_longitude])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict([
        "t" => Dict("units" => "seconds"),
        "longitude" => Dict("units" => "test_units1"),
    ])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )
    resampled_src_var =
        ClimaAnalysis.resampled_as(src_var, dest_var, dim_names = ["long", "t"])
    @test resampled_src_var.data == src_var.data[1:46, :, 1:4]
    @test resampled_src_var.dims["time"] == dest_var.dims["t"]
    @test resampled_src_var.dims["long"] == dest_var.dims["longitude"]

    # Testing resampling with only one OutputVar
    resampled_var =
        ClimaAnalysis.resampled_as(src_var, t = dest_t, long = dest_longitude)
    @test resampled_var.dims == resampled_src_var.dims
    @test resampled_var.dim_attributes == resampled_src_var.dim_attributes
    @test resampled_var.data == resampled_src_var.data
    @test resampled_var.attributes == resampled_src_var.attributes

    # Resampling with out one OutputVar, but no units for the dimensions
    src_dim_attribs = OrderedDict([
        "long" => Dict("units" => ""),
        "lat" => Dict("units" => ""),
        "time" => Dict("no_units" => ""),
    ])
    src_var_no_units =
        ClimaAnalysis.remake(src_var, dim_attributes = src_dim_attribs)
    resampled_var = ClimaAnalysis.resampled_as(
        src_var_no_units,
        t = dest_t,
        long = dest_longitude,
    )
    @test resampled_var.dims == resampled_src_var.dims
    @test resampled_var.dim_attributes == src_dim_attribs
    @test resampled_var.data == resampled_src_var.data
    @test resampled_var.attributes == resampled_src_var.attributes

    # Error handling
    # Out of bound errors
    dest_long = 0.0:200.0 |> collect
    dest_data = ones(length(dest_long))
    dest_dims = OrderedDict(["long" => dest_long])
    dest_attribs = Dict("long_name" => "hi")
    dest_dim_attribs = OrderedDict(["long" => Dict("units" => "test_units1")])
    dest_var = ClimaAnalysis.OutputVar(
        dest_attribs,
        dest_dims,
        dest_dim_attribs,
        dest_data,
    )
    @test_throws BoundsError ClimaAnalysis.resampled_as(
        src_var,
        dest_var,
        dim_names = "longitude",
    )
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
    var_without_unitful = ClimaAnalysis.remake(
        var_with_unitful,
        attributes = Dict{String, Any}(),
        dims = Dict("long" => long),
    )

    var_empty_unit = ClimaAnalysis.remake(
        var_with_unitful,
        attributes = Dict{String, Any}("units" => ""),
        dims = Dict("long" => long),
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

@testset "_check_time_dim" begin
    # Satisfies all conditions
    time = 0.0:10.0 |> collect
    data = ones(length(time))
    dims = OrderedDict(["time" => time])
    attribs = Dict("start_date" => Dates.DateTime(2010))
    dim_attribs = OrderedDict(["time" => Dict("units" => "s")])
    time_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test time_var |> ClimaAnalysis.Var._check_time_dim |> isnothing

    # Time dimension does not exist
    no_time_var = ClimaAnalysis.average_time(time_var)
    @test_throws ErrorException ClimaAnalysis.Var._check_time_dim(no_time_var)

    # Unit for the time dimension is not second
    dim_attribs = OrderedDict(["time" => Dict("units" => "min")])
    minute_var = ClimaAnalysis.remake(time_var, dim_attributes = dim_attribs)
    @test_throws ErrorException ClimaAnalysis.Var._check_time_dim(no_time_var)

    # Start date is not present
    attribs = Dict("no start date" => "idk")
    no_date_var = ClimaAnalysis.remake(time_var, attributes = attribs)
    @test_throws ErrorException ClimaAnalysis.Var._check_time_dim(no_date_var)
end

@testset "split along dim" begin
    lat = collect(range(-89.5, 89.5, 10))
    lon = collect(range(-179.5, 179.5, 10))
    time = collect(range(0.0, 10.0, 11))
    data = reshape(
        1.0:1.0:(length(lat) * length(lon) * length(time)),
        (length(lat), length(lon), length(time)),
    )
    dims = OrderedDict(["lat" => lat, "lon" => lon, "time" => time])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # Empty case
    empty_var = ClimaAnalysis.Var._split_along_dim(var, "lon", [[]]) |> first
    @test isempty(empty_var)

    # No splitting
    no_split_var =
        ClimaAnalysis.Var._split_along_dim(var, "time", [time]) |> first
    @test var.dims == no_split_var.dims
    @test var.data == no_split_var.data
    @test var.attributes == no_split_var.attributes
    @test var.dim_attributes == no_split_var.dim_attributes

    # Split along time dimension
    time_split_vars = ClimaAnalysis.Var._split_along_dim(
        var,
        "time",
        [[0.0, 5.0], [2.0, 7.0], []],
    )
    first_var = time_split_vars[begin]
    second_var = time_split_vars[begin + 1]
    third_var = time_split_vars[end]
    @test isempty(third_var)
    @test first_var.dims["time"] == [0.0, 5.0]
    @test second_var.dims["time"] == [2.0, 7.0]
    @test first_var.dims["lon"] == lon
    @test second_var.dims["lon"] == lon
    @test first_var.dims["lat"] == lat
    @test second_var.dims["lat"] == lat
    @test first_var.data == data[:, :, [1, 6]]
    @test second_var.data == data[:, :, [3, 8]]
    @test first_var.dim_attributes == var.dim_attributes
    @test second_var.dim_attributes == var.dim_attributes

    # Split along latitude dimension
    lat_split_vars = ClimaAnalysis.Var._split_along_dim(
        var,
        "lat",
        [lat[begin:5], lat[6:end]],
    )
    first_var = lat_split_vars[begin]
    second_var = lat_split_vars[end]
    @test first_var.dims["lat"] == lat[begin:5]
    @test second_var.dims["lat"] == lat[6:end]
    @test first_var.dims["lon"] == lon
    @test second_var.dims["lon"] == lon
    @test first_var.dims["time"] == time
    @test second_var.dims["time"] == time
    @test first_var.data == data[begin:5, :, :]
    @test second_var.data == data[6:end, :, :]
    @test first_var.dim_attributes == var.dim_attributes
    @test second_var.dim_attributes == var.dim_attributes
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

    # Check season is available as an attributes
    @test MAM.attributes["season"] == "MAM"
    @test JJA.attributes["season"] == "JJA"
    @test DJF.attributes["season"] == "DJF"

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

@testset "split_by_season_across_time" begin
    lon = collect(range(-179.5, 179.5, 36))
    lat = collect(range(-89.5, 89.5, 18))
    time = [0.0]
    push!(time, 2_678_400.0) # correspond to 2024-2-1
    push!(time, 5_184_000.0) # correspond to 2024-3-1
    push!(time, 7_862_400.0) # correspond to 2024-4-1
    push!(time, 10_454_400.0) # correspond to 2024-5-1
    push!(time, 13_132_800.0) # correspond to 2024-6-1
    push!(time, 15_724_800.0) # correspond to 2024-7-1
    push!(time, 18_403_200.0) # correspond to 2024-8-1
    push!(time, 21_081_600.0) # correspond to 2024-9-1
    push!(time, 36_720_000.0) # correspond to 2025-3-1
    data = reshape(
        1.0:1.0:(length(lat) * length(time) * length(lon)),
        (length(lat), length(time), length(lon)),
    )
    dims = OrderedDict(["lat" => lat, "time" => time, "lon" => lon])
    attribs = Dict("long_name" => "hi", "start_date" => "2024-1-1")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    seasonal_vars = ClimaAnalysis.split_by_season_across_time(var)
    @test length(seasonal_vars) == 6
    @test isempty(seasonal_vars[5])

    # Cheak season
    @test seasonal_vars[1].attributes["season"] == "DJF"
    @test seasonal_vars[2].attributes["season"] == "MAM"
    @test seasonal_vars[3].attributes["season"] == "JJA"
    @test seasonal_vars[4].attributes["season"] == "SON"
    @test seasonal_vars[6].attributes["season"] == "MAM"

    # Check arrays for time is correct
    @test seasonal_vars[1].dims["time"] == [0.0, 2_678_400.0]
    @test seasonal_vars[2].dims["time"] ==
          [5_184_000.0, 7_862_400.0, 10_454_400.0]
    @test seasonal_vars[3].dims["time"] ==
          [13_132_800.0, 15_724_800.0, 18_403_200.0]
    @test seasonal_vars[4].dims["time"] == [21_081_600.0]
    @test seasonal_vars[6].dims["time"] == [36_720_000.0]

    # Check data itself
    @test seasonal_vars[1].data == data[:, 1:2, :]
    @test seasonal_vars[2].data == data[:, 3:5, :]
    @test seasonal_vars[3].data == data[:, 6:8, :]
    @test seasonal_vars[4].data == data[:, 9:9, :]
    @test seasonal_vars[6].data == data[:, 10:10, :]
end

@testset "split_by_month" begin
    lon = collect(range(-179.5, 179.5, 36))
    lat = collect(range(-89.5, 89.5, 18))
    dates = [
        Dates.DateTime(2010, 1),
        Dates.DateTime(2011, 1),
        Dates.DateTime(2010, 2),
        Dates.DateTime(2010, 12),
    ]
    t = ClimaAnalysis.Utils.date_to_time.(Dates.DateTime(2010, 1), dates)

    data = reshape(
        1.0:1.0:(length(lat) * length(t) * length(lon)),
        (length(lat), length(t), length(lon)),
    )
    dims = OrderedDict(["lat" => lat, "time" => t, "lon" => lon])
    attribs = Dict("long_name" => "hi", "start_date" => "2010-1-1")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "time" => Dict("units" => "s"),
        "lon" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    monthly_vars = ClimaAnalysis.split_by_month(var)
    @test length(monthly_vars) == 12
    for i in 3:11
        @test isempty(monthly_vars[i])
    end

    # Cheak month
    for i in (1, 2, 12)
        @test monthly_vars[i].attributes["month"] == Dates.monthname(i)
    end

    # Check arrays for time is correct
    @test monthly_vars[1].dims["time"] == [0.0, 31536000.0]
    @test monthly_vars[2].dims["time"] == [2678400.0]
    @test monthly_vars[12].dims["time"] == [28857600.0]

    # Check data itself
    @test monthly_vars[1].data == data[:, 1:2, :]
    # Reshape to add extra dimension back
    @test monthly_vars[2].data == reshape(data[:, 3, :], (18, 1, 36))
    @test monthly_vars[12].data == reshape(data[:, 4, :], (18, 1, 36))
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
    @test ClimaAnalysis.units(squared_error_var) == "(kg)^2"

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
    @test ClimaAnalysis.units(var_unitful) == "(kg^2/m)^2"
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

    # Check that `dates` converts
    @test ClimaAnalysis.dates(var_s) == time_arr
    # Check that `dates` throws an error where there's no time
    @test_throws ErrorException ClimaAnalysis.dates(
        ClimaAnalysis.OutputVar(
            Dict("long" => collect(1.0:10.0)),
            collect(1.0:10.0),
        ),
    )

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
    var_lonlat = ClimaAnalysis.remake(var_latlon, dims = dims, data = data)

    land_var_lonlat = ClimaAnalysis.apply_landmask(var_lonlat)
    ocean_var_lonlat = ClimaAnalysis.apply_oceanmask(var_lonlat)
    land_var_latlon = ClimaAnalysis.apply_landmask(var_latlon)
    ocean_var_latlon = ClimaAnalysis.apply_oceanmask(var_latlon)
    @test isequal(land_var_lonlat.data, land_var_latlon.data')
    @test isequal(ocean_var_lonlat.data, ocean_var_latlon.data')

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
    @test isequal(land_var.data |> transpose, land_var_lonlat.data)
    @test isequal(ocean_var.data |> transpose, ocean_var_lonlat.data)

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
    land_var = ClimaAnalysis.Var.LAND_MASK
    land_var = ClimaAnalysis.replace(land_var, NaN => 0.0)
    ocean_var = ClimaAnalysis.Var.OCEAN_MASK
    ocean_var = ClimaAnalysis.replace(ocean_var, NaN => 0.0)
    land_var = ClimaAnalysis.set_units(land_var, "idk")
    ocean_var = ClimaAnalysis.set_units(ocean_var, "idk")
    data_zero = zeros(land_var.data |> size)
    zero_var = ClimaAnalysis.OutputVar(
        land_var.attributes,
        land_var.dims,
        land_var.dim_attributes,
        data_zero,
    )

    # Trim data because periodic boundary condition on the edges
    @test isequal(
        ClimaAnalysis.bias(
            land_var,
            zero_var,
            mask = ClimaAnalysis.apply_oceanmask,
        ).data[
            begin:(end - 1),
            :,
        ] |> A -> replace(A, NaN => 0.0),
        data_zero[begin:(end - 1), :],
    )
    @test isequal(
        ClimaAnalysis.bias(
            ocean_var,
            zero_var,
            mask = ClimaAnalysis.apply_landmask,
        ).data[
            begin:(end - 1),
            :,
        ] |> A -> replace(A, NaN => 0.0),
        data_zero[begin:(end - 1), :],
    )

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
    ] |> A -> replace(A, NaN => 0.0) == data_zero[begin:(end - 1), :]
    @test ClimaAnalysis.squared_error(
        ocean_var,
        zero_var,
        mask = ClimaAnalysis.apply_landmask,
    ).data[
        begin:(end - 1),
        :,
    ] |> A -> replace(A, NaN => 0.0) == data_zero[begin:(end - 1), :]

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

@testset "Concatenate OutputVars" begin
    # Initialize 3D OutputVar (lon, lat, time)
    times = [
        Dates.value(
            Dates.Second(Dates.DateTime(2010, i) - Dates.DateTime(2010, 1)),
        ) for i in 1:12
    ]
    lon = [-180.0, 0.0, 180.0]
    lat = [-90.0, -45.0, 0.0, 45.0, 90.0]
    data = cat(
        (i * ones(1, length(lon), length(lat)) for i in eachindex(times))...,
        dims = 1,
    )
    dims = OrderedDict(["time" => times, "lon" => lon, "lat" => lat])
    attribs = Dict(
        "long_name" => "LONG_NAME",
        "short_name" => "shrt_nm",
        "start_date" => "2010-1-1",
        "test_key" => "test_val",
        "units" => "kg",
    )
    dim_attribs = OrderedDict(["time" => Dict("units" => "s")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    # Test with common situations where you would need to concat
    # Split by season across time
    seasons = ClimaAnalysis.split_by_season_across_time(var)
    MAM = cat(seasons[2], seasons[5], dims = "time")
    @test MAM.dim_attributes == dim_attribs
    @test MAM.data == cat(data[3:5, :, :], data[12:12, :, :], dims = 1)
    @test MAM.dims == OrderedDict([
        "time" => [times[i] for i in (3, 4, 5, 12)],
        "lon" => lon,
        "lat" => lat,
    ])
    @test MAM.attributes == Dict(
        "long_name" => "LONG_NAME concatenated LONG_NAME",
        "short_name" => "shrt_nm",
        "start_date" => "2010-1-1",
        "units" => "kg",
    )

    # Window
    window_12 = ClimaAnalysis.window(
        var,
        "time",
        left = Dates.DateTime(2010, 1),
        right = Dates.DateTime(2010, 2),
    )
    window_45 = ClimaAnalysis.window(
        var,
        "time",
        left = Dates.DateTime(2010, 4),
        right = Dates.DateTime(2010, 5),
    )
    window_1245 = cat(window_12, window_45, dims = "time")
    @test window_1245.dim_attributes == dim_attribs
    @test window_1245.data == cat(data[1:2, :, :], data[4:5, :, :], dims = 1)
    @test window_1245.dims == OrderedDict([
        "time" => [times[i] for i in (1, 2, 4, 5)],
        "lon" => lon,
        "lat" => lat,
    ])
    @test window_1245.attributes == Dict(
        "long_name" => "LONG_NAME concatenated LONG_NAME",
        "short_name" => "shrt_nm",
        "start_date" => "2010-1-1",
        "units" => "kg",
    )

    # Error handling
    # Dimension does not exist
    @test_throws ErrorException cat(var, dims = "pfull")

    # Number of dimensions are not the same
    time = 0.0:10.0 |> collect
    data = ones(length(time))
    dims = OrderedDict(["time" => time])
    attribs = Dict("start_date" => Dates.DateTime(2010))
    dim_attribs = OrderedDict(["time" => Dict("units" => "s")])
    time_var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException cat(time_var, var, dims = "time")

    # Order of dimensions are not the same
    diff_order_var = permutedims(var, ("long", "lat", "t"))
    @test_throws ErrorException cat(diff_order_var, var, dims = "time")

    # Dimension units are not the same
    dim_attribs = OrderedDict(["time" => Dict("units" => "min")])
    diff_dim_units_var = ClimaAnalysis.remake(var, dim_attributes = dim_attribs)
    @test_throws ErrorException cat(diff_dim_units_var, var, dims = "time")

    # Values of dimensions are not the same
    dims = OrderedDict([
        "time" => times,
        "lon" => [-180.0, 10.0, 180.0],
        "lat" => lat,
    ])
    diff_vals_var = ClimaAnalysis.remake(var, dims = dims)
    @test_throws ErrorException cat(diff_vals_var, var, dims = "time")

    # Units of the data are not the same
    attribs = Dict(
        "long_name" => "LONG_NAME",
        "short_name" => "shrt_nm",
        "start_date" => "2010-1-1",
        "test_key" => "test_val",
        "units" => "kg^-1",
    )
    diff_units_var = ClimaAnalysis.remake(var, attributes = attribs)
    @test_throws ErrorException cat(diff_units_var, var, dims = "time")

    # Short name is different
    attribs = Dict(
        "long_name" => "LONG_NAME",
        "short_name" => "diff_shrt_nm",
        "start_date" => "2010-1-1",
        "test_key" => "test_val",
        "units" => "kg",
    )
    diff_shortname_var = ClimaAnalysis.remake(var, attributes = attribs)
    @test_throws ErrorException cat(diff_shortname_var, var, dims = "time")
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
    dim_attribs =
        OrderedDict(["lat" => Dict(), "lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.remake(var, dim_attributes = dim_attribs)
    ClimaAnalysis.set_dim_units!(var, "lat", "degrees")
    @test ClimaAnalysis.dim_units(var, "lat") == "degrees"

    # Dimension is not present in dim_attribs
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.remake(var, dim_attributes = dim_attribs)
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
    ones_var = ClimaAnalysis.remake(var, data = data_ones)

    mask_fn = ClimaAnalysis.make_lonlat_mask(
        var,
        set_to_val = isnan,
        true_val = 0.0,
        false_val = 1.0,
    )
    var1 = mask_fn(ones_var)
    @test var1.data[1] == 0.0

    # Check with a different longitude name
    ones_var = ClimaAnalysis.remake(
        ones_var,
        dims = OrderedDict(["latitude" => lat, "longitude" => lon]),
        dim_attributes = OrderedDict(["longitude" => Dict("units" => "deg")]),
    )
    # ones_var.dims = OrderedDict(["latitude" => lat, "longitude" => lon])
    # ones_var.dim_attributes = OrderedDict(["longitude" => Dict("units" => "deg")])
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
    var = ClimaAnalysis.remake(var, data = data, dims = dims)
    @test_throws ErrorException ClimaAnalysis.make_lonlat_mask(var)

    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    t = collect(range(1, 2, 2))
    data = ones(length(lat), length(lon), length(t))
    dims = OrderedDict(["lat" => lat, "lon" => lon, "time" => t])
    var = ClimaAnalysis.remake(var, data = data, dims = dims)
    @test_throws ErrorException ClimaAnalysis.make_lonlat_mask(var)
end

@testset "Reverse dimensions" begin
    lat = reverse(collect(range(-89.5, 89.5, 180)))
    lon = collect(range(-179.5, 179.5, 360))
    data = reshape(collect(1:(360 * 180)), (180, 360))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lon" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test isnothing(ClimaAnalysis.Var._make_interpolant(dims, data))

    # Test reverse_dim
    reverse_var = ClimaAnalysis.reverse_dim(var, "lat")
    @test reverse_var.dims["lat"] == reverse(lat)
    @test reverse_var.data == reverse(data, dims = 1)

    # Test reverse_dim!
    ClimaAnalysis.reverse_dim!(var, "lat")
    @test var.dims["lat"] == collect(range(-89.5, 89.5, 180))
    @test var.data ==
          reverse(reshape(collect(1:(360 * 180)), (180, 360)), dims = 1)

    # Error handling
    @test_throws ErrorException ClimaAnalysis.reverse_dim(var, "pressure_level")
    @test_throws ErrorException ClimaAnalysis.reverse_dim!(
        var,
        "pressure_level",
    )

    arb_dim = reshape(collect(range(-89.5, 89.5, 16)), (4, 4))
    data = collect(1:16)
    dims = OrderedDict(["arb_dim" => arb_dim])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["arb_dim" => Dict("units" => "idk")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test_throws ErrorException ClimaAnalysis.reverse_dim(var, "arb_dim")
    @test_throws ErrorException ClimaAnalysis.reverse_dim!(var, "arb_dim")
end

@testset "Convert units of dimensions" begin
    # Convert units of one dimension
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lat" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    var_units = ClimaAnalysis.convert_dim_units(
        var,
        "lat",
        "rads",
        conversion_function = x -> x * π / 180.0,
    )
    @test ClimaAnalysis.latitudes(var_units) == lat .* π ./ 180.0
    @test ClimaAnalysis.dim_units(var_units, "lat") == "rads"

    # Error handling
    # Dimension doesn't exists
    @test_throws ErrorException ClimaAnalysis.convert_dim_units(
        var,
        "dim do not exists",
        "idk",
        conversion_function = identity,
    )

    # Units do not exists for the dimension
    @test_throws ErrorException ClimaAnalysis.convert_dim_units(
        var,
        "lon",
        "idk",
        conversion_function = identity,
    )

    # Conversion function is not supplied
    @test_throws ErrorException ClimaAnalysis.convert_dim_units(
        var,
        "lat",
        "rads",
    )
end

@testset "Show" begin
    lat = collect(range(-89.5, 89.5, 180))
    lon = collect(range(-179.5, 179.5, 360))
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lat" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\nData defined over:\n  lat with 180 elements (-89.5 to 89.5)\n  lon with 360 elements (-179.5 to 179.5)"

    # Reverse lat
    lat = reverse(lat)
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    var = ClimaAnalysis.remake(var, dims = dims)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\nData defined over:\n  lat with 180 elements (89.5 to -89.5)\n  lon with 360 elements (-179.5 to 179.5)"

    # Unsorted lat
    lat[1] = -1000.0
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    var = ClimaAnalysis.remake(var, dims = dims)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\nData defined over:\n  lat with 180 elements (not sorted)\n  lon with 360 elements (-179.5 to 179.5)"

    lat = [1.0]
    lon = [2.0]
    data = ones(length(lat), length(lon))
    dims = OrderedDict(["lat" => lat, "lon" => lon])
    var = ClimaAnalysis.remake(var, data = data, dims = dims)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\nData defined over:\n  lat with 1 element (1.0)\n  lon with 1 element (2.0)"

    # Empty OutputVar
    dims = OrderedDict{String, Vector{Float64}}()
    data = Float64[]
    empty_var = ClimaAnalysis.OutputVar(dims, data)
    @test sprint(show, empty_var) ==
          "Attributes:\nDimension attributes:\nData defined over:"

    # OutputVar with string in dim_attributes as the value as opposed to another dictionary
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lat))
    dims = OrderedDict(["lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict([
        "lat" => Dict(
            "units" => "deg",
            "wacky" => Dict("something" => "another thing"),
        ),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\n    wacky => Dict(\"something\" => \"another thing\")\nData defined over:\n  lat with 180 elements (-89.5 to 89.5)"

    # OutputVar with 0 elements in one of the dimenisions
    lat = Float64[]
    data = ones(length(lat))
    dims = OrderedDict(["lat" => lat])
    attribs = Dict("long_name" => "hi")
    dim_attribs = OrderedDict(["lat" => Dict("units" => "deg")])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    @test sprint(show, var) ==
          "Attributes:\n  long_name => hi\nDimension attributes:\n  lat:\n    units => deg\nData defined over:\n  lat with 0 element"
end
