# Visualizing `OutputVar`s

This page is under construction, in the meantime, consult [`Visualize`](@ref).

### Masking part of the output in `GeoMakie`

When performing ocean or land simulations, it is often convenient to hide the
other component (e.g., hide the ocean and focus on the continents). For
`GeoMakie` plots, there is a direct way to accomplish this. In this section, we
discuss this feature.

The main `GeoMakie` plots are [`Visualize.contour2D_on_globe!`](@ref) and
[`Visualize.heatmap2D_on_globe!`](@ref). Both these functions take a `mask` argument. By
default, `mask` is `nothing`, meaning that the entire output is displayed on the
globe. Alternatively, `mask` can be a collection of polygons that can be plotted
with [`Makie.poly`](https://docs.makie.org/v0.21/reference/plots/poly).
`ClimaAnalysis` comes with the most important ones [`Visualize.oceanmask`](@ref) and
[`Visualize.landmask`](@ref), to hide ocean and continents respectively.

For example, suppose `var` it the variable we want to plot with a ocean mask
```julia
import ClimaAnalysis.Visualize: contour2D_on_globe!, oceanmask
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import GeoMakie
import CairoMakie

fig = CairoMakie.Figure()

contour2D_on_globe!(fig,
                    var,
                    mask = oceanmask(),
                    more_kwargs = Dict(:mask => ca_kwargs(color = :blue)),
                   )

CairoMakie.save("myfigure.pdf", fig)
```

In this example, we plotted `var` on the globe and overplotted a blue ocean.
`ca_kwargs` (`Utils.kwargs`) is a convenience function to pass keyword arguments
more easily.

The output might look something like:

![oceanmask](./assets/oceanmask.png)

### Plotting bias

After [computing the bias](@ref bias) between observational and simulation data, you may
want to plot the bias and display information such as the root mean squared error (RMSE) and
the global bias in the plot. To do this, you use the function [`plot_bias_on_globe!(fig, sim,
obs)`](@ref Visualize.plot_bias_on_globe!). In the example below, we plot the bias between our
simulation and some observations stored in `ta_1d_average.nc`.

```julia
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!
import GeoMakie
import CairoMakie

obs_var = ClimaAnalysis.OutputVar("ta_1d_average.nc")
sim_var = ClimaAnalysis.get(ClimaAnalysis.simdir("simulation_output"), "ta")

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig, sim_var, obs_var)
CairoMakie.save("myfigure.pdf", fig)
```

The output produces something like:

![biasplot](./assets/bias_plot.png)

We can also plot the bias using an ocean mask. This also means we compute the bias only
over land.

```julia
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!, oceanmask
import GeoMakie
import CairoMakie

obs_var = ClimaAnalysis.OutputVar("ta_1d_average.nc")
sim_var = ClimaAnalysis.get(ClimaAnalysis.simdir("simulation_output"), "ta")

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig,
                    var,
                    mask = oceanmask(),
                    more_kwargs = Dict(:mask => ca_kwargs(color = :blue)),
                   )
CairoMakie.save("myfigure.pdf", fig)
```

The output produces something like:

![biasplot_oceanmask](./assets/bias_plot_oceanmask.png)

We can also plot the bias using a custom mask generated from `make_lonlat_mask`.

!!! note "Passing a masking function for `mask`"
    ClimaAnalysis do not support mask keyword arguments for masking functions. If you want
    the values of the mask to not show in a plot, then pass `true_val = NaN` as a keyword
    argument to `make_lonlat_mask`. The color of `NaN` is controlled by the keyword
    `nan_color` which can be passed for the plotting function (`:plot`).

    Note that if the backend is CairoMakie, then the keyword `nan_color` does nothing. See
    this [issue](https://github.com/MakieOrg/Makie.jl/issues/4524).

```julia
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!, oceanmask
import GeoMakie
import CairoMakie

mask_var = ClimaAnalysis.OutputVar("ocean_mask.nc")
mask_fn = ClimaAnalysis.generate_lonlat_mask(mask_var, NaN, 1.0)

obs_var = ClimaAnalysis.OutputVar("ta_1d_average.nc")
sim_var = ClimaAnalysis.get(ClimaAnalysis.simdir("simulation_output"), "ta")

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig, var, mask = mask_fn)
CairoMakie.save("myfigure.pdf", fig)
```

The output produces something like:

![bias_with_custom_mask_plot](./assets/plot_bias_with_custom_mask.png)

### Makie integration

```@setup makie_integration
import ClimaAnalysis.Template:
    TemplateVar,
    make_template_var,
    add_attribs,
    add_dim,
    add_data,
    one_to_n_data,
    initialize

lat = -90.0:1.0:90.0
lon = -180:1.0:180.0
data2d = [cos(x / 90) * sin(y / 180) for x in lat, y in lon]
var2d =
    TemplateVar() |>
    add_dim("lat", lat, units = "degrees") |>
    add_dim("lon", lon, units = "degrees") |>
    add_attribs(; long_name = "Air temperature", short_name = "ta", units = "K") |>
    add_data(; data = data2d) |>
    initialize
```

In versions of ClimaAnalysis after v0.5.22, Makie functions can be used natively
with `OutputVar`s. In the example below, we use the `CairoMakie` backend and the
`CairoMakie.plot` function to plot `var`. The type of the plot is automatically
determined by the number of dimensions of the `OutputVar`. This supports
`OutputVar` representing one dimensional or two dimensional data.

```@example makie_integration
import ClimaAnalysis
import CairoMakie

CairoMakie.plot(var2d; colormap = :vik10)
```

You can also use other plotting functions in `Makie` such as `heatmap` and
`lines`. In this example, we use `heatmap` to plot a `OutputVar` with the
longitude and latitude dimensions and `lines` to plot a weighted
latitude-averaged `OutputVar`. You can expect most plotting functions from
`Makie` to work out of the box.

```@example makie_integration
fig = CairoMakie.Figure(; size = (1000, 400))
# If no axis is provided, then a default axis is added
CairoMakie.heatmap(fig[1,1], var2d)

# You can make your own axis for more customization
lat_weighted_var = ClimaAnalysis.weighted_average_lat(var2d)
long_name = ClimaAnalysis.long_name(lat_weighted_var)
lon_units = ClimaAnalysis.dim_units(lat_weighted_var, "longitude")
short_name = ClimaAnalysis.short_name(lat_weighted_var)
var_units = ClimaAnalysis.units(lat_weighted_var)
ax = CairoMakie.Axis(fig[1, 2], title = long_name, xlabel = "Longitude ($lon_units)", ylabel = "$short_name ($var_units)")
CairoMakie.lines!(ax, lat_weighted_var; color = :blue)
fig
```
