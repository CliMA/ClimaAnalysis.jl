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

    fig = Makie.Figure()

    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D_globe.png")
    Makie.save(output_name, fig)

    fig2 = Makie.Figure()

    ClimaAnalysis.Visualize.contour2D_on_globe!(fig2, var2D)

    output_name = joinpath(tmp_dir, "test_contours2D_globe.png")
    Makie.save(output_name, fig2)

end
