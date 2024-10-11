using Test
import ClimaAnalysis
import Makie
import CairoMakie

using OrderedCollections

@testset "MakieExt" begin

    tmp_dir = mktempdir(cleanup = false)
    @info "Tempdir", tmp_dir

    long = 0.0:180.0 |> collect
    lat = 0.0:90.0 |> collect
    time = 0.0:10.0 |> collect

    data3D = reshape(1.0:(91 * 181 * 11), (11, 181, 91))
    dims3D = OrderedDict(["time" => time, "lon" => long, "lat" => lat])
    attribs = Dict([
        "long_name" => "The quick brown fox jumps over the lazy dog. The quick brown fox.",
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

    # Intialize another figure to see if plotting with multiple figures initialized is
    # possible
    fig = Makie.Figure()
    fig1 = Makie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.heatmap2D!(fig, var3D)

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

    ClimaAnalysis.Visualize.heatmap2D!(fig, var2D)

    output_name = joinpath(tmp_dir, "test2D.png")
    Makie.save(output_name, fig)

    # Test with a GridLayout
    fig = Makie.Figure()
    layout = fig[1, 1] = Makie.GridLayout()
    ClimaAnalysis.Visualize.heatmap2D!(layout, var2D)

    output_name = joinpath(tmp_dir, "test2D_gd.png")
    Makie.save(output_name, fig)


    # New figure
    fig = Makie.Figure()
    @test_throws ErrorException ClimaAnalysis.Visualize.sliced_heatmap!(
        fig,
        var3D,
    )

    cut = Dict("time" => 1)
    ClimaAnalysis.Visualize.sliced_heatmap!(fig, var3D, cut)
    output_name = joinpath(tmp_dir, "test3D_sliced.png")
    Makie.save(output_name, fig)

    # New figure
    fig = Makie.Figure()
    ClimaAnalysis.Visualize.heatmap!(fig, var3D; time = 1)
    output_name = joinpath(tmp_dir, "test3D_sliced_kwargs.png")
    Makie.save(output_name, fig)

    @test_throws ErrorException ClimaAnalysis.Visualize.line_plot1D!(fig, var3D)

    data1D = reshape(1.0:(91), (91))
    dims1D = OrderedDict(["lat" => lat])
    attribs = Dict([
        "long_name" => "The quick brown fox jumps over the lazy dog. The quick brown fox.",
        "short_name" => "name",
        "units" => "bob",
    ])
    path = "a/b/c"
    dim_attributes1D = OrderedDict(["lat" => Dict(["units" => "degrees"])])
    var1D = ClimaAnalysis.OutputVar(attribs, dims1D, dim_attributes1D, data1D)

    # Intialize another figure to see if plotting with multiple figures initialized is
    # possible
    fig = Makie.Figure()
    fig2 = Makie.Figure()
    ClimaAnalysis.Visualize.line_plot1D!(fig, var1D)
    output_name = joinpath(tmp_dir, "test1D.png")
    Makie.save(output_name, fig)

    fig = Makie.Figure()
    cut = Dict("lon" => 30)
    ClimaAnalysis.Visualize.sliced_line_plot!(fig, var2D, cut)
    output_name = joinpath(tmp_dir, "test2D_sliced.png")
    Makie.save(output_name, fig)

    fig = Makie.Figure()
    ClimaAnalysis.Visualize.line_plot!(fig, var2D; lon = 30)
    output_name = joinpath(tmp_dir, "test2D_sliced_kwargs.png")
    Makie.save(output_name, fig)

    # Test plot!
    fig = Makie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var2D; lon = 30)
    output_name = joinpath(tmp_dir, "test_plot2D_sliced.png")
    Makie.save(output_name, fig)

    fig = Makie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var3D; lon = 30, time = 1)
    output_name = joinpath(tmp_dir, "test_plot3D_sliced.png")
    Makie.save(output_name, fig)

    fig = Makie.Figure()
    ClimaAnalysis.Visualize.plot!(fig, var3D; time = 1)
    output_name = joinpath(tmp_dir, "test_plot3D_sliced_once.png")
    Makie.save(output_name, fig)

    # Test passing more_kwargs

    fig = Makie.Figure()
    ClimaAnalysis.Visualize.plot!(
        fig,
        var3D;
        time = 1,
        more_kwargs = Dict(:cb => [:vertical => :false]),
    )

    fig = Makie.Figure()
    ClimaAnalysis.Visualize.plot!(
        fig,
        var3D;
        time = 1,
        more_kwargs = Dict(
            :plot => ClimaAnalysis.Utils.kwargs(colormap = :inferno),
        ),
    )
    output_name = joinpath(tmp_dir, "test_plot3D_sliced_once.png")
    Makie.save(output_name, fig)

    # Test dim_on_y
    fig = Makie.Figure()
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
    Makie.save(output_name, fig)

    # Test overriding title, xlabel, and ylabel
    fig = Makie.Figure()
    ClimaAnalysis.Visualize.heatmap2D!(
        fig,
        var2D,
        more_kwargs = Dict(
            :axis => ClimaAnalysis.Utils.kwargs(
                title = "My title: The quick brown fox jumps over the lazy dog. The quick brown fox.",
                xlabel = "My xlabel",
                ylabel = "My ylabel",
            ),
        ),
    )

    output_name = joinpath(tmp_dir, "test2D_title.png")
    Makie.save(output_name, fig)

    # Plotting box plot
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "CliMA")
    rmse_var["CliMA", :] = [12.0, 12.0, 11.0, 14.0, 6.0]
    ClimaAnalysis.add_unit!(
        rmse_var,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
        ),
    )
    rmse_var[2, 5] = 4.0

    rmse_var1 = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var1 = ClimaAnalysis.add_model(rmse_var1, "CliMA", "test_model")
    rmse_var1["CliMA", :] = [12.0, 12.0, 11.0, 14.0, 6.0]
    rmse_var1["test_model", :] = [12.0, 12.0, 11.0, 14.0, 6.0]
    ClimaAnalysis.add_unit!(
        rmse_var1,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
            "test_model" => "units",
        ),
    )
    rmse_var1[2, 5] = 4.0

    fig = Makie.Figure(; size = (800, 300 * 3 + 400), fontsize = 20)
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_var,
        model_names = ["CliMA"],
        best_and_worst_category_name = "ANN",
    )
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_var,
        model_names = ["CliMA"],
        ploc = (2, 1),
        best_and_worst_category_name = "ANN",
    )
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_var1,
        model_names = ["CliMA", "ACCESS-ESM1-5"],
        ploc = (3, 1),
        best_and_worst_category_name = "ANN",
    )
    output_name = joinpath(tmp_dir, "test_boxplots.png")
    Makie.save(output_name, fig)

    # Plotting leaderboard
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "CliMA")
    rmse_var[:, :] = [
        [10.0 11.0 12.0 13.0 14.0]
        [36.0 37.0 38.0 39.0 30.0]
        [11.0 12.0 13.0 14.0 15.0]
    ]
    ClimaAnalysis.add_unit!(
        rmse_var,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
        ),
    )

    rmse_var1 = ClimaAnalysis.read_rmses(csv_file_path, "ta1")
    rmse_var1 = ClimaAnalysis.add_model(rmse_var1, "CliMA")
    rmse_var1[:, :] = [
        [6.0 7.0 8.0 9.0 10.0]
        [11.0 12.0 13.0 14.0 15.0]
        [1.0 2.0 3.0 4.0 5.0]
    ]
    ClimaAnalysis.add_unit!(
        rmse_var1,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
        ),
    )

    rmse_var2 = ClimaAnalysis.read_rmses(csv_file_path, "ta2")
    rmse_var2 = ClimaAnalysis.add_model(rmse_var2, "CliMA")
    rmse_var2[:, :] = [
        [0.5 1.0 1.5 2.0 2.5]
        [6.0 7.0 8.0 9.0 10.0]
        [11.0 12.0 13.0 14.0 15.0]
    ]
    ClimaAnalysis.add_unit!(
        rmse_var2,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
        ),
    )


    # Normalized RMSEs should improve going from ta to ta1 to ta2 for CliMA model
    # Colors for ta of both models should be similar (close to 1)
    # Colors for ta2 of best model should be greener (improve) from ta1 to ta2
    # for CliMA, the normalized RMSEs from greatest to least should be ta1, ta, and ta2
    fig = Makie.Figure(; fontsize = 20)
    ClimaAnalysis.Visualize.plot_leaderboard!(
        fig,
        rmse_var,
        rmse_var1,
        rmse_var2,
        best_category_name = "ANN",
    )
    output_name = joinpath(tmp_dir, "test_leaderboard.png")
    Makie.save(output_name, fig)

    # Plot box plots and leaderboard in one plot
    rmse_vars = (rmse_var, rmse_var1, rmse_var2)
    fig = Makie.Figure(; size = (800, 300 * 3 + 400), fontsize = 20)
    for i in 1:3
        ClimaAnalysis.Visualize.plot_boxplot!(
            fig,
            rmse_vars[i],
            ploc = (i, 1),
            best_and_worst_category_name = "ANN",
        )
    end
    ClimaAnalysis.Visualize.plot_leaderboard!(
        fig,
        rmse_vars...,
        best_category_name = "ANN",
        ploc = (4, 1),
    )
    output_name = joinpath(tmp_dir, "test_boxplot_and_leaderboard.png")
    Makie.save(output_name, fig)

    # Plotting box plot with NaN
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "CliMA")
    ClimaAnalysis.add_unit!(
        rmse_var,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "CliMA" => "units",
        ),
    )
    rmse_var[2, 5] = 10.0
    fig = Makie.Figure(; fontsize = 20)
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_var,
        model_names = ["CliMA"],
        best_and_worst_category_name = "ANN",
    )
    output_name = joinpath(tmp_dir, "test_boxplot_nan.png")
    Makie.save(output_name, fig)

    fig = Makie.Figure(; fontsize = 20)
    ClimaAnalysis.Visualize.plot_leaderboard!(
        fig,
        rmse_var,
        model_names = ["CliMA"],
        best_category_name = "ANN",
    )
    output_name = joinpath(tmp_dir, "test_leaderboard_nan.png")
    Makie.save(output_name, fig)

    # Testing with long name
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(
        rmse_var,
        "long_name_name_name_name_name_name_name",
    )
    ClimaAnalysis.add_unit!(
        rmse_var,
        Dict(
            "ACCESS-ESM1-5" => "units",
            "ACCESS-CM2" => "units",
            "long_name_name_name_name_name_name_name" => "units",
        ),
    )
    rmse_var[2, 5] = 10.0
    fig = Makie.Figure(; fontsize = 20)
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_var,
        model_names = ["long_name_name_name_name_name_name_name"],
        best_and_worst_category_name = "ANN",
    )
    output_name = joinpath(tmp_dir, "test_boxplot_long_name.png")
    Makie.save(output_name, fig)

    # Test error handling for plot_leaderboard
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var1 = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var1 = ClimaAnalysis.add_category(rmse_var1, "hi")
    rmse_var2 = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var2 = ClimaAnalysis.add_category(rmse_var2, "hello")
    @test_throws ErrorException ClimaAnalysis.Visualize.plot_leaderboard!(
        fig,
        rmse_var1,
        rmse_var2,
        model_names = ["CliMA"],
        best_category_name = "ANN",
    )
end
