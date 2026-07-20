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

For example, suppose `var` is the variable we want to plot with an ocean mask

```@setup visualize_mask
import ClimaAnalysis
import ClimaAnalysis.Template:
    TemplateVar, add_attribs, add_dim, add_data, initialize

lon = collect(-180.0:1.0:180.0)
lat = collect(-90.0:1.0:90.0)
data = [288.0 - 40.0 * sind(la)^2 + 5.0 * cosd(lo) for lo in lon, la in lat]
var =
    TemplateVar() |>
    add_dim("lon", lon, units = "degrees_east") |>
    add_dim("lat", lat, units = "degrees_north") |>
    add_attribs(long_name = "Air temperature", short_name = "ta", units = "K") |>
    add_data(data = data) |>
    initialize
```

```@example visualize_mask
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

fig
```

In this example, we plotted `var` on the globe and overplotted a blue ocean.
`ca_kwargs` (`Utils.kwargs`) is a convenience function to pass keyword arguments
more easily.

### Plotting bias

After [computing the bias](@ref bias) between observational and simulation data, you may
want to plot the bias and display information such as the root mean squared error (RMSE) and
the global bias in the plot. To do this, you use the function [`plot_bias_on_globe!(fig, sim,
obs)`](@ref Visualize.plot_bias_on_globe!). In the example below, we plot the bias between our
simulation data (`sim_var`) and observational data (`obs_var`), where both are
`OutputVar`s defined over longitude and latitude.

```@setup bias_plots
import ClimaAnalysis
import ClimaAnalysis.Template:
    TemplateVar, add_attribs, add_dim, add_data, initialize

lon = collect(range(-179.5, 179.5, 360))
lat = collect(range(-89.5, 89.5, 180))
sim_data = [280.0 + 20.0 * cosd(la) + 5.0 * sind(lo) for lo in lon, la in lat]
obs_data = sim_data .+ [3.0 * cosd(2 * la) * sind(3 * lo) for lo in lon, la in lat]

sim_var =
    TemplateVar() |>
    add_dim("lon", lon, units = "degrees_east") |>
    add_dim("lat", lat, units = "degrees_north") |>
    add_attribs(long_name = "Air temperature", short_name = "ta", units = "K") |>
    add_data(data = sim_data) |>
    initialize

obs_var =
    TemplateVar() |>
    add_dim("lon", lon, units = "degrees_east") |>
    add_dim("lat", lat, units = "degrees_north") |>
    add_attribs(long_name = "Air temperature", short_name = "ta", units = "K") |>
    add_data(data = obs_data) |>
    initialize

mask_var = ClimaAnalysis.Var.LANDSEA_MASK
```

```@example bias_plots
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!
import GeoMakie
import CairoMakie

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig, sim_var, obs_var)
fig
```

We can also plot the bias using an ocean mask. This also means we compute the bias only
over land.

```@example bias_plots
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!, oceanmask
import ClimaAnalysis.Utils: kwargs as ca_kwargs
import GeoMakie
import CairoMakie

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig,
                    sim_var,
                    obs_var,
                    mask = oceanmask(),
                    more_kwargs = Dict(:mask => ca_kwargs(color = :blue)),
                   )
fig
```

We can also plot the bias using a custom mask generated from `generate_lonlat_mask`. In
the example below, `mask_var` is a `OutputVar` defined over longitude and latitude, whose
data contains only zeros and ones, where the ones indicate land and the zeros indicate
the ocean. The mask generated from
`ClimaAnalysis.generate_lonlat_mask(mask_var, NaN, 1.0)` replaces values with `NaN`
wherever the mask is zero and keeps values the same wherever the mask is one. As a
result, the bias is plotted only over land.

!!! note "Passing a masking function for `mask`"
    ClimaAnalysis does not support mask keyword arguments for masking functions. If you
    want the values of the mask to not show in a plot, then pass `NaN` for the `zero_to`
    or `one_to` positional arguments of `generate_lonlat_mask`. The color of `NaN` is
    controlled by the keyword `nan_color` which can be passed for the plotting function
    (`:plot`).

    Note that if the backend is CairoMakie, then the keyword `nan_color` does nothing. See
    this [issue](https://github.com/MakieOrg/Makie.jl/issues/4524).

```@example bias_plots
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!
import GeoMakie
import CairoMakie

mask_fn = ClimaAnalysis.generate_lonlat_mask(mask_var, NaN, 1.0)

fig = CairoMakie.Figure()
# Wrap mask_fn in a function so that it is recognized as a masking function
plot_bias_on_globe!(fig, sim_var, obs_var, mask = v -> mask_fn(v))
fig
```

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
data2d = [288.0 - 40.0 * sind(x)^2 + 5.0 * cosd(y) for x in lat, y in lon]
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
`CairoMakie.plot` function to plot `var2d`. The type of the plot is automatically
determined by the number of dimensions of the `OutputVar`. This supports
`OutputVar` representing one dimensional or two dimensional data.

```@example makie_integration
import ClimaAnalysis
import CairoMakie

CairoMakie.plot(var2d; colormap = :thermal)
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
