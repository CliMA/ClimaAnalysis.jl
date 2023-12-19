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
files, and organizes them.

As of version `0.1.0`, `ClimaAnalysis` uses file names to identify files and
variables. In this, `ClimaAnalysis` assumes that the default names for outputs
are used in `ClimaAtmos` (i.e.,
`<short_name>_<reduction_time>_<reduction_type>.nc`, as in `ta_1h_max.nc`, or
`<short_name>_inst.nc`).

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
    inst
- ta
    average (3.0h)
    max (4.0h, 3.0h)
    min (3.0h)
- ts
    max (1.0h)
```
Now, you can access any given variable
``` julia
ta_max = get(simdir; short_name = "t12", reduction = "max", period = "3.0h")
```
`ta_max` is a ` OutputVar`, a type that contains the variable as
well as some metadata.

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
ts_max_lat_averaged_sliced = ClimaAnalysis.slice_z(ts_max_lat_averaged, 8_000)
```
Alternatively, you can also call `ClimaAnalysis.slice(ts_max_lat_averaged_sliced, time = 8_000)`.
Now,
``` julia
ts_max_lat_averaged_sliced.dims =
OrderedCollections.OrderedDict{String, Vector{Float32}} with 2 entries:
  "time" => [10800.0, 21600.0, 32400.0, 43200.0]
  "lon"  => [-180.0, -177.989, -175.978, -173.966, -171.955, -169.944, -167.933, -165.922…
```

We can directly visualize these quantities.

### `Visualize`

If `CairoMakie` is available, `ClimaAnalysis` can be used for plotting.
Importing `CairoMakie` and `ClimaAnalysis` in the same session automatically
loads the necessary `ClimaAnalysis` plotting modules.

If we want to make a heatmap for `ta_max` at time of 100 s at altitude `z` of 30000 meters:

``` julia
import CairoMakie
import ClimaAnalysis: Visualize

fig = CairoMakie.Figure(resolution = (400, 600))

viz.plot!(
  fig,
  ta_max,
  time = 100.0,
  z = 30_000.0
)

CairoMakie.save("ta_max.png", fig)
```

If we want to have a line plot, we can simply add another argument (e.g., `lat =
30`), to slice through that value.

If you want to customize some of the properties, you can pass the `more_kwargs`
 to the `plot!` function. `more_kwargs` is a dictionary that can contain
 additional arguments to the `Axis` (`:axis`), `plot` (`:plot`), and `Colorbar`
 (`:cb`) functions. `more_kwargs` is a Dictionary that maps the symbols `:axis`,
 `:plot`, and `:cb` to their additional arguments. For instance, to choose the
 colormap for the colorbar to viridis
``` julia
viz.plot!(
  fig,
  ta_max,
  time = 100.0,
  z = 30_000.0,
  more_kwargs = Dict(:cb => [:colormap => :viridis])
)
```
Note the `Symbol` in colormap!. `:cb` has to be a mapping of `Symbol`s and
values. `ClimaAnalysis` has a convenience function `kwargs` to more easily pass
down the keyword arguments avoiding this step. With that, the above example becomes
``` julia
import ClimaAnalysis: Utils
viz.plot!(
  fig,
  ta_max,
  time = 100.0,
  z = 30_000.0,
  plot_kwargs = Dict(:cb => Utils.kwargs(colormap = :inferno))
)
```
With `Utils.kwargs`, you can just pass the arguments as you would pass them to
the constructor.


