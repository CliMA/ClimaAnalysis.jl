using Test
import ClimaAnalysis
import Makie
import GeoMakie

using OrderedCollections

@testset "MakieExt" begin
    tmp_dir = mktempdir(cleanup = false)
    @info "Tempdir", tmp_dir

    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect

    data2D = reshape(1.0:(91 * 181), (181, 91))
    dims2D = OrderedDict(["lon" => long, "lat" => lat])
    attribs = Dict([
        "long_name" => "The quick brown fox jumps over the lazy dog. The quick brown fox.",
        "short_name" => "name",
        "units" => "bob",
    ])
    path = "a/b/c"
    dim_attributes2D = OrderedDict([
        "lon" => Dict(["units" => "degrees"]),
        "lat" => Dict(["units" => "degrees"]),
    ])
    var2D = ClimaAnalysis.OutputVar(attribs, dims2D, dim_attributes2D, data2D)

    fig = Makie.Figure()

    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D_globe.png")
    Makie.save(output_name, fig)

    fig2 = Makie.Figure()

    ClimaAnalysis.Visualize.contour2D_on_globe!(
        fig2,
        var2D,
        more_kwargs = Dict(
            :plot =>
                ClimaAnalysis.Utils.kwargs(colormap = Makie.colorschemes[:vik]),
        ),
    )

    output_name = joinpath(tmp_dir, "test_contours2D_globe.png")
    Makie.save(output_name, fig2)

    # Test cmap 
    MakieExt = Base.get_extension(ClimaAnalysis, :MakieExt)
    test_cmap = MakieExt._constrained_cmap(
        Makie.colorschemes[:vik],
        0.0,
        15000.0 + (5000.0 / 3.0),
        mid = 5000.0,
        categorical = true,
    )

    fig3 = Makie.Figure()

    ClimaAnalysis.Visualize.contour2D_on_globe!(
        fig3,
        var2D,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = test_cmap),
        ),
    )

    output_name = joinpath(tmp_dir, "test_contours2D_globe_with_test_cmap.png")
    Makie.save(output_name, fig3)

    test_cmap = MakieExt._constrained_cmap(
        range(Makie.colorant"red", stop = Makie.colorant"green", length = 15),
        0.0,
        15000.0 + (5000.0 / 3.0),
    )

    fig4 = Makie.Figure()

    ClimaAnalysis.Visualize.contour2D_on_globe!(
        fig4,
        var2D,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = test_cmap),
        ),
    )

    output_name = joinpath(tmp_dir, "test_contours2D_globe_with_test_cmap2.png")
    Makie.save(output_name, fig4)

end
