# API

```@meta
CurrentModule = ClimaAnalysis
```

## Sim

```@docs
Sim.SimDir
Base.get
Sim.available_vars
Sim.available_reductions
Sim.available_periods
```

## Var

```@docs
Var.OutputVar
Var.read_var
Var.is_z_1D
Var.short_name
Var.long_name
Var.units
Var.slice
Var.average_lat
Var.weighted_average_lat
Var.average_lon
Var.average_x
Var.average_y
Var.average_xy
Var.average_time
Var.window
Var.arecompatible
Var.center_longitude!
Var.times
Var.dates
Var.longitudes
Var.latitudes
Var.altitudes
Var.time_name
Var.date_name
Var.longitude_name
Var.latitude_name
Var.altitude_name
Var.has_time
Var.has_date
Var.has_longitude
Var.has_latitude
Var.has_altitude
Var.conventional_dim_name
Var.dim_units
Var.range_dim
Var.resampled_as
```

## Utilities

For development and not

```@docs
Utils.match_nc_filename
Utils.squeeze
Utils.nearest_index
Utils.kwargs
Utils.seconds_to_prettystr
Utils.warp_string
Utils.split_by_season
```

## Atmos

```@docs
Atmos
Atmos.to_pressure_coordinates
```

## Makie

```@docs
Visualize.heatmap2D!
Visualize.sliced_heatmap!
Visualize.heatmap!
Visualize.line_plot1D!
Visualize.sliced_line_plot!
Visualize.sliced_plot!
Visualize.plot!
```

## GeoMakie

```@docs
Visualize.contour2D_on_globe!
Visualize.heatmap2D_on_globe!
```
