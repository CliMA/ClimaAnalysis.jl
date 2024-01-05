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
slice
slice_general
slice_x
slice_y
slice_z
slice_lon
slice_lat
slice_time
average_lat
average_lon
average_time
window
```


## Utilities

For development and not

```@docs
Utils.match_nc_filename
Utils.squeeze
Utils.nearest_index
Utils.kwargs
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
