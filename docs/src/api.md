# API

```@meta
CurrentModule = ClimaAnalysis
```

## SimDir

```@docs
SimDir
Base.get
available_vars
available_reductions
available_periods
```

## OutputVar

```@docs
OutputVar
read_var
is_z_1D
short_name
long_name
units
slice
average_lat
average_lon
average_time
window
arecompatible
```


## Utilities

For development and not

```@docs
Utils.match_nc_filename
Utils.squeeze
Utils.nearest_index
Utils.kwargs
Utils.seconds_to_prettystr
```

## CairoMakie

```@docs
Visualize.heatmap2D!
Visualize.sliced_heatmap!
Visualize.heatmap!
Visualize.line_plot1D!
Visualize.sliced_line_plot!
Visualize.sliced_plot!
Visualize.plot!
```

## CairoMakie

```@docs
Visualize.heatmap2D_on_globe!
```
