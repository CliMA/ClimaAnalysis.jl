using Test
import ClimaAnalysis
import CairoMakie

using OrderedCollections

@testset "CairoMakieExt" begin

    tmp_dir = mktempdir(cleanup = false)
    @info "Tempdir", tmp_dir

    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data3D = reshape(1.0:(91 * 181 * 10), (10, 181, 91))
    dims3D = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    attribs = Dict([
        "long_name" => "My name",
        "short_name" => "name",
        "units" => "bob",
    ])
    path = "a/b/c"
    dim_attributes3D = OrderedDict([
        "time" => Dict(["units" => "s"]),
        "lon" => Dict(["units" => "degrees"]),
        "lat" => Dict(["units" => "degrees"]),
    ])
    var3D =
        ClimaAnalysis.OutputVar(attribs, dims3D, dim_attributes3D, data3D, path)

    fig = CairoMakie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.contour_plot2D!(
        fig,
        var3D,
    )

    data2D = reshape(1.0:(91 * 181), (181, 91))
    dims2D = OrderedDict(["lon" => long, "lat" => lat])
    attribs = Dict([
        "long_name" => "My name",
        "short_name" => "name",
        "units" => "bob",
    ])
    path = "a/b/c"
    dim_attributes2D = OrderedDict([
        "lon" => Dict(["units" => "degrees"]),
        "lat" => Dict(["units" => "degrees"]),
    ])
    var2D =
        ClimaAnalysis.OutputVar(attribs, dims2D, dim_attributes2D, data2D, path)

    ClimaAnalysis.Visualize.contour_plot2D!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D.png")
    CairoMakie.save(output_name, fig)

    # New figure
    fig = CairoMakie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.sliced_contour_plot!(
        fig,
        var3D,
    )

    cut = Dict("time" => 1)
    ClimaAnalysis.Visualize.sliced_contour_plot!(fig, var3D, cut)
    output_name = joinpath(tmp_dir, "test3D_sliced.png")
    CairoMakie.save(output_name, fig)

    # New figure
    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.contour_plot!(fig, var3D; time = 1)
    output_name = joinpath(tmp_dir, "test3D_sliced_kwargs.png")
    CairoMakie.save(output_name, fig)
end
