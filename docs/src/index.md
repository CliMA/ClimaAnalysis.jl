# ClimaAnalysis

`ClimaAnalysis`, your one-stop-shop for all your CliMA analysis needs.

## Quick start

### `SimDir`

Suppose you ran a `ClimaAtmos` simulation, and the output is saved in the folder
`simulation_output`. The first step in using `ClimaAnalysis` is to instantiate a
`SimDir`:
``` julia
import ClimaAnalysis

simdir = ClimaAnalysis.SimDir("simulation_output")
```
`ClimaAnalysis.SimDir` scans the `simulation_output`, finds all the output
files, and organizes them. Consider the following directories:
```
global_diagnostics/
├── output_0001
│   ├── ta_3.0h_average.nc
│   ├── ta_3.0h_max.nc
│   ├── ta_3.0h_min.nc
│   ├── ta_4.0h_max.nc
│   ├── ts_1.0h_max.nc
│   ├── ua_6.0h_average.nc
│   └── va_2.0h_average.nc
├── output_0002
│   ├── ua_6.0h_average.nc
│   └── va_2.0h_average.nc
└── output_0003
    ├── ts_1.0h_max.nc
    ├── ua_6.0h_average.nc
    └── va_2.0h_average.nc
```

As of version `0.1.0`, `ClimaAnalysis` uses file names to identify files and
variables. In this, `ClimaAnalysis` assumes that the default names for outputs
are used in `ClimaAtmos` (i.e.,
`<short_name>_<reduction_time>_<reduction_type>.nc`, as in `ta_1h_max.nc`, or
`<short_name>_1d_inst.nc`).

Once you have a `SimDir`, you can inspect the output. For example, to find what
variables are available:
``` julia-repl
julia> println(summary(simdir))
Output directory: simulation_output
Variables:
- va
    average (2.0h)
- ua
    average (6.0h)
- orog
    inst (1.0d)
- ta
    average (3.0h)
    max (4.0h, 3.0h)
    min (3.0h)
- ts
    max (1.0h)
```
Now, you can access any given variable
``` julia
ts_max = get(simdir; short_name = "ts", reduction = "max", period = "3.0h")
```

`ts_max` is a ` OutputVar`, a type that contains the variable as well as some
metadata. When there is only one combination `short_name/reduction/period`, the
function `get` can be used with `get(simdir, short_name)` (e.g., `get(simdir,
"orog")` in the previous example). In this case, you can get multiple variables
at the same time with
```julia
orog, ua, va = get(simdir, "orog", "ua", "va")
```

If there are more files with the same combination of short name, reduction, and
period, then the function `get` automatically stitch the `.nc` files together
along the time dimension.
``` julia
# Stitch `ua_6.0h_average.nc` in output_0001, output_0002, and output_0003
ua_average = get(simdir; short_name = "ua", reduction = "average", period = "6.0h")

# Stitch `ts_1.0h_max.nc` in output_0001 and output_0003
ts_average = get(simdir; short_name = "ts", reduction = "max", period = "1.0h")
```

The order of files in the directory tree traversed top-down determines the order of the
files when stitching them. Stitching the datasets is not possible when the start dates are
not the same, the names of the time dimension are not the same across the datasets, the time
dimension does not exist, or the times are not in sequential order.

Let us learn about ` OutputVar`s

### `OutputVar`

` OutputVar`s contain the raw data (in `.data`), the attributes read from the
file, and the information regarding the dimension over which the variable is
defined.
``` julia-repl
julia> ts_max.dims
OrderedCollections.OrderedDict{String, Vector{Float32}} with 4 entries:
  "time" => [10800.0, 21600.0, 32400.0, 43200.0]
  "lon"  => [-180.0, -177.989, -175.978, -173.966, -171.955, -169.944, -167.933, -165.922…
  "lat"  => [-80.0, -77.9747, -75.9494, -73.924, -71.8987, -69.8734, -67.8481, -65.8228, …
  "z"    => [0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, …
```
Here we have the dimensions and their values. The dimensions are ordered as in
the file, so that the first index of `.data` is `time`, and so on.

We can find the attributes of the dimensions in `.attributes`:
``` julia-repl
julia> ts_max.dim_attributes["lon"]
  "lon"  => Dict("units"=>"degrees_east")
```

Some of the attributes are exposed with function calls. For example
``` julia-repl
julia> long_name(ts_max)
  Surface Temperature, max within 1.0 Hour(s)
```
These function use the attributes in the NetCDF files. When not available, empty strings are returned.

Given an `OutputVar`, we can perform manipulations. For instance, we can take
the average over latitudes:
``` julia
ts_max_lat_averaged = ClimaAnalysis.average_lat(ts_max)
```
Now,
``` julia
ts_max_lat_averaged.dims =
OrderedCollections.OrderedDict{String, Vector{Float32}} with 3 entries:
  "time" => [10800.0, 21600.0, 32400.0, 43200.0]
  "lon"  => [-180.0, -177.989, -175.978, -173.966, -171.955, -169.944, -167.933, -165.922…
  "z"    => [0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, …
```
We can also take a time/altitude slice, for example, the plane with altitude closest to 8000 meters.
``` julia
ts_max_lat_averaged_sliced = ClimaAnalysis.slice(ts_max_lat_averaged, 8_000)
```
Alternatively, you can also call `ClimaAnalysis.slice(ts_max_lat_averaged_sliced, time = 8_000)`.
Now,
``` julia
ts_max_lat_averaged_sliced.dims =
OrderedCollections.OrderedDict{String, Vector{Float32}} with 2 entries:
  "time" => [10800.0, 21600.0, 32400.0, 43200.0]
  "lon"  => [-180.0, -177.989, -175.978, -173.966, -171.955, -169.944, -167.933, -165.922…
```

You can get the dimensions from standard names, for example, to find the
`times`, simply run
``` julia
times(ts_max_lat_averaged_sliced) =
4-element Vector{Float32}:
 10800.0
 21600.0
 32400.0
 43200.0
```

`OutputVar`s can be evaluated on arbitrary points. For instance
``` julia-repl
julia> ts_max([12000., 23., 45., 1200.])
```
will return the value of the maximum temperature at time 12000, longitude 23,
latitude 45, and altitude 1200. This can be used to interpolate `OutputVar`s
onto new grids.

#### Mathematical operations

`OutputVar`s support the usual mathematical operations. For instance, if
`ts_max` is an `OutputVar`, `2 * ts_max` will be an `OutputVar` with doubled
values.

For binary operations (e.g., `+, -, *, /`), `ClimaAnalysis` will check if the
operation is well defined (i.e., the two variables are defined on the physical
space). Binary operations do remove some attribute information. If two
`OutputVar`s share the same start date, then the start date will remain in the
resulting `OutputVar`s after performing binary operations on them.

Unary operations are supported to, for example `log(max(var, 1e-8))` returns a
new `OutputVar` with the functions applied to the data.

#### `Visualize`

We can directly visualize `OutputVar`s.

If [`Makie`](https://docs.makie.org/stable/) is available, `ClimaAnalysis` can
be used for plotting. Importing `Makie` and `ClimaAnalysis` in the same session
automatically loads the necessary `ClimaAnalysis` plotting modules.

If we want to make a heatmap for `ts_max` at time of 100 s at altitude `z` of 30000 meters:

``` julia
import CairoMakie
import ClimaAnalysis.Visualize as viz

fig = CairoMakie.Figure(size = (400, 600))

viz.plot!(
  fig,
  ts_max,
  time = 100.0,
  z = 30_000.0
)

CairoMakie.save("ts_max.png", fig)
```

If we want to have a line plot, we can simply add another argument (e.g., `lat =
30`), to slice through that value.

If you want to customize some of the properties, you can pass the `more_kwargs`
to the `plot!` function. `more_kwargs` is a dictionary that can contain
additional arguments to the `Axis` (`:axis`), `plot` (`:plot`), and `Colorbar`
(`:cb`) functions. `more_kwargs` is a Dictionary that maps the symbols `:axis`,
`:plot`, and `:cb` to their additional arguments. For instance, to choose the
alpha value of the plot, the label of the colorbar, and the subtitle, you can
do the following:
``` julia
viz.plot!(
    fig,
    ts_max,
    time = 100.0,
    z = 30_000.0,
    more_kwargs = Dict(
        :plot => Dict(:alpha => 0.5),
        :cb => Dict(:label => "My label"),
        :axis => Dict(:subtitle => "My subtitle"),
    ),
)
```
Note the `Symbol` in plot, cb, and axis!. `:plot`, `:cb`, and `:axis` have to be
a mapping of `Symbol`s and values. `ClimaAnalysis` has a convenience function
`kwargs` to more easily pass down the keyword arguments avoiding this step. With
that, the above example becomes
``` julia
import ClimaAnalysis.Utils : kwargs as ca_kwargs
viz.plot!(
  fig,
  ts_max,
  time = 100.0,
  z = 30_000.0,
  more_kwargs = Dict(
      :plot => ca_kwargs(alpha = 0.5),
      :cb => ca_kwargs(label = "My label"),
      :axis => ca_kwargs(subtitle = "My subtitle"),
  ),
)
```
With `Utils.kwargs`, you can just pass the arguments as you would pass them to
the constructor.

If you need more control over the placement of plots, you can pass
`Makie.GridLayout` objects to the plotting functions. For example,
``` julia
using CairoMakie

fig = Figure()
layout = fig[1, 2] = GridLayout()

viz.plot!(
  layout,
  ts_max,
  time = 100.0,
  z = 30_000.0,
  more_kwargs = Dict(
      :plot => ca_kwargs(alpha = 0.5),
      :cb => ca_kwargs(label = "My label"),
      :axis => ca_kwargs(subtitle = "My subtitle"),
  ),
)
```

When you pass a `GridLayout`, the optional argument `p_loc` refers to the
placement within the layout. When you pass a `Figure`, it refers to the
placement within the figure.

If you have `GeoMakie` and are working on a variable defined on a `long-lat`
grid, you can directly plot on a projected global surface. For that, load
`GeoMakie` and use the `heatmap2D_on_globe!` function.
