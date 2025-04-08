using Test

import ClimaAnalysis: Template
import OrderedCollections: OrderedDict

@testset "Lazy evaluation" begin
    lazy_eval = Template.make_lazy((x, y, z) -> x + y * z, 1, 2, 3)
    @test lazy_eval() == 7

    lazy_eval = Template.make_lazy((x, y, z) -> x + y * z, 1, 2)
    @test lazy_eval(3) == 7

    lazy_eval = Template.make_lazy((; a = 1, b = 2) -> a + b, a = 100)
    @test lazy_eval() == 102

    lazy_eval = Template.make_lazy((; a = 1, b = 2) -> a + b)
    @test lazy_eval() == 3

    lazy_eval = Template.make_lazy(
        (x, y, z; a = 1, b = 2) -> x * y + z * a + b,
        4,
        6;
        b = 10,
    )
    @test lazy_eval(12, a = 50) == 634
end

@testset "Template" begin
    var =
        Template.TemplateVar() |>
        Template.add_attribs(long_name = "hi") |>
        Template.add_time_dim(dim = collect(1.0:30.0)) |>
        Template.add_lon_dim(units = "deg") |>
        Template.add_lat_dim(units = "deg") |>
        Template.one_to_n_data(collected = false) |>
        Template.initialize

    @test var.attributes == Dict("long_name" => "hi")
    @test var.dim_attributes == OrderedDict(
        "time" => Dict("units" => "seconds"),
        "longitude" => Dict("units" => "deg"),
        "latitude" => Dict("units" => "deg"),
    )
    @test var.data == reshape(1.0:(30.0 * 361.0 * 181.0), (30, 361, 181))
    @test var.dims == OrderedDict(
        "time" => collect(1.0:30.0),
        "longitude" => collect(range(-180.0, 180.0, 361)),
        "latitude" => collect(range(-90.0, 90.0, 181)),
    )

    var =
        Template.TemplateVar() |>
        Template.add_attribs(long_name = "Test", short_name = "test") |>
        Template.add_time_dim(
            dim = [0.0, 1.0, 2.0],
            units = "ab",
            cool = "idk1",
        ) |>
        Template.add_lon_dim(dim = [3.0, 4.0], units = "cd", rad = "idk2") |>
        Template.add_lat_dim(
            dim = [5.0, 6.0, 7.0],
            units = "ef",
            wa = "idk3",
        ) |>
        Template.add_pfull_dim(
            dim = [8.0, 9.0, 10.0],
            units = "gh",
            sd = "idk4",
        ) |>
        Template.add_z_dim(dim = [11.0, 12.0], units = "ij", ef = "idk5") |>
        Template.zeros_data(data_type = Float32) |>
        Template.initialize
    @test var.data == zeros(Float32, 3, 2, 3, 3, 2)
    @test var.dims == OrderedDict(
        "time" => [0.0, 1.0, 2.0],
        "longitude" => [3.0, 4.0],
        "latitude" => [5.0, 6.0, 7.0],
        "pfull" => [8.0, 9.0, 10.0],
        "z" => [11.0, 12.0],
    )
    @test var.attributes == Dict("short_name" => "test", "long_name" => "Test")
    @test var.dim_attributes == OrderedDict(
        "time" => Dict("units" => "ab", "cool" => "idk1"),
        "longitude" => Dict("units" => "cd", "rad" => "idk2"),
        "latitude" => Dict("units" => "ef", "wa" => "idk3"),
        "pfull" => Dict("units" => "gh", "sd" => "idk4"),
        "z" => Dict("units" => "ij", "ef" => "idk5"),
    )


    # Test with add_dim
    var =
        Template.TemplateVar() |>
        Template.add_dim("time", [0.0, 1.0]) |>
        Template.initialize
    @test var.dims["time"] == [0.0, 1.0]
    @test isempty(var.attributes)
    @test var.dim_attributes["time"] == Dict()
    @test var.data == [1.0, 2.0]

    # TemplateVar with ones data
    var =
        Template.TemplateVar() |>
        Template.add_time_dim(dim = collect(1.0:5.0)) |>
        Template.ones_data() |>
        Template.initialize
    @test var.data == ones(5)

    # TemplateVar with ones to n data with collected = true
    var =
        Template.TemplateVar() |>
        Template.add_time_dim(dim = collect(1.0:5.0)) |>
        Template.one_to_n_data(collected = true) |>
        Template.initialize
    @test var.data == collect(1.0:5.0)

    # Empty TemplateVar
    @test Template.TemplateVar() |> Template.initialize |> isempty

    # TemplateVar with no data function
    var =
        Template.TemplateVar() |>
        Template.add_time_dim(dim = collect(1.0:5.0)) |>
        Template.initialize
    @test var.data == 1:5
end

@testset "Template factory" begin
    var =
        Template.make_template_var("t", "long", "lat", "pfull", "z") |>
        Template.initialize
    @test isempty(var.attributes)
    dim_sizes = (11, 361, 181, 11, 11)
    @test var.data == reshape(1:prod(dim_sizes), dim_sizes...)
    @test var.dims == OrderedDict(
        "time" => collect(0.0:10.0),
        "lon" => collect(-180.0:180.0),
        "lat" => collect(-90.0:90.0),
        "pfull" => collect(0.0:10.0),
        "z" => collect(0.0:10.0),
    )
    @test var.dim_attributes == OrderedDict(
        "time" => Dict("units" => "seconds"),
        "lon" => Dict("units" => "degrees_north"),
        "lat" => Dict("units" => "degrees_east"),
        "pfull" => Dict("units" => "Pa"),
        "z" => Dict("units" => "m"),
    )

    # Check adding the same dimensions does not permute the dimensions
    var =
        Template.make_template_var("t", "long", "lat"; abc = "def") |>
        Template.add_lon_dim(name = "lon", dim = [1.0, 2.0]) |>
        Template.add_time_dim(name = "t", dim = [3.0, 4.0]) |>
        Template.add_lat_dim(name = "lat", dim = [5.0, 6.0]) |>
        Template.initialize

    @test var.dims == OrderedDict(
        "t" => [3.0, 4.0],
        "lon" => [1.0, 2.0],
        "lat" => [5.0, 6.0],
    )
    dim_sizes = (2, 2, 2)
    @test var.data == reshape(1:prod(dim_sizes), dim_sizes...)
    @test var.dim_attributes == OrderedDict(
        "t" => Dict("units" => "seconds"),
        "lon" => Dict("units" => "degrees_north"),
        "lat" => Dict("units" => "degrees_east"),
    )
    @test var.attributes == Dict("abc" => "def")
end
