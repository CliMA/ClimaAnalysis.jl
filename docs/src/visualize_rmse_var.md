# Visualizing `RMSEVariable`s

Instead of computing summary statistics, it may be more helpful to plot a box plot or a
heatmap. `ClimaAnalysis` provides the functions `plot_boxplot!` and `plot_leaderboard!`
to help visualize the root mean squared errors (RMSEs) in a `RMSEVariable`.

The function [`Visualize.plot_boxplot!`](@ref) makes a box plot for each category in the
`RMSEVariable`. The best model and worst model and any other models in `model_names` are
plotted. When finding the best and worst single models, any models in `model_names` will be
excluded. The category to find the best and worst model defaults to "ANN", but can be
changed using the parameter `best_and_worst_category_name`.

The function [`Visualize.plot_leaderboard!`](@ref) makes a heatmap of the
RMSEs between the variables of interest and the categories. The best model for each variable
of interest and the models in `model_names` are shown in the heatmap. Similar to
`plot_boxplot!`, the category to find the best model defaults to "ANN", but can be changed
using the parameter `best_category_name`. The values of the heatmap are normalized by
dividing over the median model's RMSEs for each variable.

```@setup plotting
import ClimaAnalysis
import CairoMakie

csv_file_path = "./data/test_csv.csv"
rmse_var_ta = ClimaAnalysis.read_rmses(csv_file_path, "ta")
rmse_var_ta = ClimaAnalysis.add_model(rmse_var_ta, "CliMA", "test1", "test2")
rmse_var_ta[:, :] = [
    [10.0 11.0 12.0 13.0 14.0]
    [36.0 37.0 38.0 39.0 30.0]
    [11.0 12.0 13.0 14.0 15.0]
    [13.0 13.0 13.0 13.0 15.0]
    [24.0 24.0 24.0 24.0 24.0]
]
ClimaAnalysis.add_unit!(
    rmse_var_ta,
    Dict(
        "ACCESS-ESM1-5" => "K",
        "ACCESS-CM2" => "K",
        "CliMA" => "K",
        "test1" => "K",
        "test2" => "K",
    ),
)

rmse_var_pr = ClimaAnalysis.read_rmses(csv_file_path, "pr")
rmse_var_pr = ClimaAnalysis.add_model(rmse_var_pr, "CliMA")
rmse_var_pr[:, :] = [
    [6.0 7.0 8.0 9.0 10.0]
    [11.0 12.0 13.0 14.0 15.0]
    [1.0 2.0 3.0 4.0 11.0]
]
ClimaAnalysis.add_unit!(
    rmse_var_pr,
    Dict(
        "ACCESS-ESM1-5" => "kg m^-2 s^-1",
        "ACCESS-CM2" => "kg m^-2 s^-1",
        "CliMA" => "kg m^-2 s^-1",
    ),
)

rmse_var_ha = ClimaAnalysis.read_rmses(csv_file_path, "ha")
rmse_var_ha = ClimaAnalysis.add_model(rmse_var_ha, "CliMA")
rmse_var_ha[:, :] = [
    [0.5 1.0 1.5 2.0 2.5]
    [6.0 7.0 8.0 9.0 10.0]
    [11.0 12.0 13.0 14.0 7.0]
]
ClimaAnalysis.add_unit!(
    rmse_var_ha,
    Dict(
        "ACCESS-ESM1-5" => "m^2 s^-2",
        "ACCESS-CM2" => "m^2 s^-2",
        "CliMA" => "m^2 s^-2",
    ),
)
```

```@example plotting
import ClimaAnalysis
import CairoMakie

# Plot box plots
rmse_vars = (rmse_var_ta, rmse_var_pr, rmse_var_ha)
fig = CairoMakie.Figure(; size = (800, 300 * 3 + 400), fontsize = 20)
for i in 1:3
    ClimaAnalysis.Visualize.plot_boxplot!(
        fig,
        rmse_vars[i],
        ploc = (i, 1),
        best_and_worst_category_name = "ANN",
    )
end

# Plot leaderboard
ClimaAnalysis.Visualize.plot_leaderboard!(
    fig,
    rmse_vars...,
    best_category_name = "ANN",
    ploc = (4, 1),
)
CairoMakie.save("./assets/boxplot_and_leaderboard.png", fig)

nothing # hide
```

![box plot](./assets/boxplot_and_leaderboard.png)
