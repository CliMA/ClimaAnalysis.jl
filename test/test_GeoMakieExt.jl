using Test
import ClimaAnalysis
import Makie
import CairoMakie
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

    # Intialize another figure to see if plotting with multiple figures initialized is
    # possible
    fig = Makie.Figure()
    fig1 = Makie.Figure()

    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D_globe.png")
    Makie.save(output_name, fig)

    fig2 = Makie.Figure()

    ClimaAnalysis.Visualize.contour2D_on_globe!(
        fig2,
        var2D,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = :vik),
        ),
    )

    output_name = joinpath(tmp_dir, "test_contours2D_globe.png")
    Makie.save(output_name, fig2)

    # Test cmap
    test_cmap = ClimaAnalysis.Visualize._constrained_cmap(
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

    test_cmap = ClimaAnalysis.Visualize._constrained_cmap(
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

    # Test with oceanmask
    fig5 = Makie.Figure()

    ClimaAnalysis.Visualize.heatmap2D_on_globe!(
        fig5,
        var2D,
        mask = ClimaAnalysis.Visualize.oceanmask(),
        more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :blue)),
    )

    output_name = joinpath(tmp_dir, "test_contours2D_globe_with_oceanmask.png")
    Makie.save(output_name, fig5)

    # Test plot_bias
    fig6 = Makie.Figure()

    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = collect(reshape(-32400:32399, (360, 180))) ./ (32399.0 / 5.0)
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    data_zero = zeros(length(lon), length(lat))
    var_zero = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_zero)

    ClimaAnalysis.Visualize.plot_bias_on_globe!(fig6, var, var_zero)
    output_name = joinpath(tmp_dir, "plot_bias.png")
    Makie.save(output_name, fig6)

    # Test plot bias with keyword arguments
    fig7 = Makie.Figure()
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig7,
        var,
        var_zero,
        more_kwargs = Dict(
            :axis => Dict(:title => "no title"),
            :plot => Dict(:extendhigh => nothing),
        ),
    )
    output_name = joinpath(tmp_dir, "plot_bias_kwargs.png")
    Makie.save(output_name, fig7)

    # Test plots with apply_landmask and apply_oceanmask
    fig8 = Makie.Figure()

    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = ones(length(lon), length(lat))
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
    land_var = var |> ClimaAnalysis.apply_landmask
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig8, land_var)
    output_name = joinpath(tmp_dir, "plot_apply_land_mask.png")
    Makie.save(output_name, fig8)

    fig9 = Makie.Figure()
    ocean_var = var |> ClimaAnalysis.apply_oceanmask
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig9, ocean_var)
    output_name = joinpath(tmp_dir, "plot_apply_ocean_mask.png")
    Makie.save(output_name, fig9)

    # Test bias plots with apply_landmask and apply_oceanmask
    fig10 = Makie.Figure()

    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = collect(reshape(-32400:32399, (360, 180))) ./ (32399.0 / 5.0)
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "idk", "short_name" => "ta", "units" => "K")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    data_zero = zeros(length(lon), length(lat))
    var_zero = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_zero)
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig10,
        var,
        var_zero;
        mask = ClimaAnalysis.Visualize.landmask(),
        more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :red)),
    )
    output_name = joinpath(tmp_dir, "plot_bias_landmask.png")
    Makie.save(output_name, fig10)

    fig11 = Makie.Figure()
    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig11,
        var,
        var_zero;
        mask = ClimaAnalysis.Visualize.oceanmask(),
        more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :blue)),
    )
    output_name = joinpath(tmp_dir, "plot_bias_oceanmask.png")
    Makie.save(output_name, fig11)

    # Test bias plots with a mask that is not landmask() or oceanmask()
    fig12 = Makie.Figure()
    land_mask_modified = first(ClimaAnalysis.Visualize.landmask(), 50)

    ClimaAnalysis.Visualize.heatmap2D_on_globe!(
        fig12,
        var2D,
        mask = land_mask_modified,
        more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :red)),
    )
    output_name = joinpath(tmp_dir, "plot_modified_mask.png")
    Makie.save(output_name, fig12)

    # # Make plot with mask generating function
    fig13 = Makie.Figure()
    ncpath = joinpath(@__DIR__, "sample_nc/test_gpp_mask.nc")
    var = ClimaAnalysis.OutputVar(ncpath, "gpp")
    var.attributes["short_name"] = "gpp"
    var = ClimaAnalysis.replace(var, missing => NaN)
    var = ClimaAnalysis.slice(var, time = ClimaAnalysis.times(var) |> first)
    mask_fn = ClimaAnalysis.make_lonlat_mask(
        var,
        set_to_val = x -> !isnan(x),
        true_val = 0.0,
    )
    var = mask_fn(var)
    ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig13, var)
    output_name = joinpath(tmp_dir, "plot_custom_mask.png")
    Makie.save(output_name, fig13)

    # Make bias with mask generating function
    fig14 = Makie.Figure()

    lon = collect(range(-179.5, 179.5, 360))
    lat = collect(range(-89.5, 89.5, 180))
    data = collect(reshape(-32400:32399, (360, 180))) ./ (32399.0 / 5.0)
    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs =
        Dict("long_name" => "idk", "short_name" => "short", "units" => "kg")
    dim_attribs = OrderedDict([
        "lon" => Dict("units" => "deg"),
        "lat" => Dict("units" => "deg"),
    ])
    var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

    data_zero = zeros(length(lon), length(lat))
    var_zero = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_zero)

    var_mask = ClimaAnalysis.OutputVar(ncpath, "gpp")
    var_mask.attributes["short_name"] = "gpp"
    var_mask = ClimaAnalysis.replace(var_mask, missing => NaN)
    var_mask = ClimaAnalysis.slice(
        var_mask,
        time = ClimaAnalysis.times(var_mask) |> first,
    )
    mask_fn = ClimaAnalysis.make_lonlat_mask(var_mask; set_to_val = isnan)

    ClimaAnalysis.Visualize.plot_bias_on_globe!(
        fig14,
        var,
        var_zero,
        mask = mask_fn,
        cmap_extrema = (-50.0, 50.0),
        # The keyword `nan_color` do not work right now for CairoMakie.
        # See https://github.com/MakieOrg/Makie.jl/issues/4524
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(nan_color = :red),
        ),
    )
    output_name = joinpath(tmp_dir, "plot_bias_with_custom_mask.png")
    Makie.save(output_name, fig14)
end
