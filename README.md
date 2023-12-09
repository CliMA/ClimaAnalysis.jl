<h1 align="center">
  ClimaAnalysis.jl
</h1>
<p align="center">
  <strong>Analyzing and visualizing ClimaAtmos simulations</strong>
</p>

`ClimaAnalysis.jl` is a Julia library to post-process and visualize `ClimaAtmos`
simulations.

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

> [!NOTE]
> As of version `0.1.0`, `ClimaAnalysis` uses file names to identify files and
> variables. In this, `ClimaAnalysis` assumes that the default names for outputs
> are used in `ClimaAtmos` (i.e.,
> `<short_name>_<reduction_time>_<reduction_type>.nc`, as in `ta_1h_max.nc`, or
> `<short_name>_inst.nc`).

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
``` julia-repl
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
``` julia-repl
ts_max_lat_averaged = ClimaAnalysis.average_lat(ts_max)
```
Now,
``` julia-repl
ts_max_lat_averaged.dims =
OrderedCollections.OrderedDict{String, Vector{Float32}} with 4 entries:
  "time" => [10800.0, 21600.0, 32400.0, 43200.0]
  "lon"  => [-180.0, -177.989, -175.978, -173.966, -171.955, -169.944, -167.933, -165.922…
  "z"    => [0.0, 5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, …
```

We can directly visualize these quantities.

### `Visualize`

If `CairoMakie` is available, `ClimaAnalysis` can be used for plotting.
Importing `CairoMakie` and `ClimaAnalysis` in the same session automatically
loads the necessary `ClimaAnalysis` plotting modules.

If we want to make a contour plot for `ts_max`:

``` julia
import CairoMakie
import ClimaAnalysis: Visualize

time = ta_max.dims["time"][end]
lon = ta_max.dims["lon"]
lat = ta_max.dims["lat"]
ZINDEX = 3
z = ta_max.dims[z][ZINDEX]

lat_units = ta_max.dim_attributes["lat"]["units"]
lon_units = ta_max.dim_attributes["lon"]["units"]
z_units = ta_max.dim_attributes["z"]["units"]
units = ta_max.attributes["units"]
long_name = ta_max.attributes["long_name"]

fig = CairoMakie.Figure(resolution = (400, 600))
title = "$long_name z = $z $z_units, t = $time $t_units"
viz.contour_plot!(
  fig;
  X = lon,
  Y = lat,
  Z = Z,
  title,
  xlabel = "lon [lon_units]",
  ylabel = "lat [lat_units]",
  colorbar_label = "ta_max [$units]",
  p_loc = (1, 1),
)

CairoMakie.save("ta_max.png", fig)
```


## Features

- Read, organize, and process NetCDF files
- Visualize contour plots with `CairoMakie`
- Apply averages and other reductions to the output variables.
