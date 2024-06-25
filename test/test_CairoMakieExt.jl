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
    var3D = ClimaAnalysis.OutputVar(attribs, dims3D, dim_attributes3D, data3D)

    fig = CairoMakie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.heatmap2D!(fig, var3D)

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
    var2D = ClimaAnalysis.OutputVar(attribs, dims2D, dim_attributes2D, data2D)

    ClimaAnalysis.Visualize.heatmap2D!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D.png")
    CairoMakie.save(output_name, fig)

    # Test with a GridLayout
    fig = CairoMakie.Figure()
    layout = fig[1, 1] = CairoMakie.GridLayout()
    ClimaAnalysis.Visualize.heatmap2D!(layout, var2D)

    output_name = joinpath(tmp_dir, "test2D_gd.png")
    CairoMakie.save(output_name, fig)


    # New figure
    fig = CairoMakie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.sliced_heatmap!(
        fig,
        var3D,
    )

    cut = Dict("time" => 1)
    ClimaAnalysis.Visualize.sliced_heatmap!(fig, var3D, cut)
    output_name = joinpath(tmp_dir, "test3D_sliced.png")
    CairoMakie.save(output_name, fig)

    # New figure
    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.heatmap!(fig, var3D; time = 1)
    output_name = joinpath(tmp_dir, "test3D_sliced_kwargs.png")
    CairoMakie.save(output_name, fig)

    @test_throws ErrorException ClimaAnalysis.Visualize.line_plot1D!(fig, var3D)

    data1D = reshape(1.0:(91), (91))
    dims1D = OrderedDict(["lat" => lat])
    attribs = Dict([
        "long_name" => "My name",
        "short_name" => "name",
        "units" => "bob",
    ])
    path = "a/b/c"
    dim_attributes1D = OrderedDict(["lat" => Dict(["units" => "degrees"])])
    var1D = ClimaAnalysis.OutputVar(attribs, dims1D, dim_attributes1D, data1D)

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.line_plot1D!(fig, var1D)
    output_name = joinpath(tmp_dir, "test1D.png")
    CairoMakie.save(output_name, fig)

    fig = CairoMakie.Figure()
    cut = Dict("lon" => 30)
    ClimaAnalysis.Visualize.sliced_line_plot!(fig, var2D, cut)
    output_name = joinpath(tmp_dir, "test2D_sliced.png")
    CairoMakie.save(output_name, fig)

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.line_plot!(fig, var2D; lon = 30)
    output_name = joinpath(tmp_dir, "test2D_sliced_kwargs.png")
    CairoMakie.save(output_name, fig)

    # Test plot!
    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var2D; lon = 30)
    output_name = joinpath(tmp_dir, "test_plot2D_sliced.png")
    CairoMakie.save(output_name, fig)

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var3D; lon = 30, time = 1)
    output_name = joinpath(tmp_dir, "test_plot3D_sliced.png")
    CairoMakie.save(output_name, fig)

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var3D; time = 1)
    output_name = joinpath(tmp_dir, "test_plot3D_sliced_once.png")
    CairoMakie.save(output_name, fig)

    # Test passing more_kwargs

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(
        fig,
        var3D;
        time = 1,
        more_kwargs = Dict(:cb => [:vertical => :false]),
    )

    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(
        fig,
        var3D;
        time = 1,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = :inferno),
        ),
    )
    output_name = joinpath(tmp_dir, "test_plot3D_sliced_once.png")
    CairoMakie.save(output_name, fig)

    # Test dim_on_y
    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.plot!(
        fig,
        var3D;
        time = 1,
        lon = 30,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = :inferno),
            :axis => ClimaAnalysis.Utils.kwargs(dim_on_y = true),
        ),
    )
    output_name = joinpath(tmp_dir, "test_plot3D_sliced_swapped.png")
    CairoMakie.save(output_name, fig)

    # Test overriding title, xlabel, and ylabel
    fig = CairoMakie.Figure()
    ClimaAnalysis.Visualize.heatmap2D!(
        fig,
        var2D,
        more_kwargs = Dict(
            :axis => ClimaAnalysis.Utils.kwargs(
                title = "My title",
                xlabel = "My xlabel",
                ylabel = "My ylabel",
            ),
        ),
    )

    output_name = joinpath(tmp_dir, "test2D_title.png")
    CairoMakie.save(output_name, fig)

end
