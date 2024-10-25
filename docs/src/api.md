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
Base.isempty(var::OutputVar)
Var.short_name
Var.long_name
Var.units
Var.has_units
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
Var.reordered_as
Var.resampled_as
Var.convert_units
Var.set_units
Var.set_dim_units!
Var.integrate_lonlat
Var.integrate_lat
Var.integrate_lon
Var.split_by_season(var::OutputVar)
Var.bias
Var.global_bias
Var.squared_error
Var.global_mse
Var.global_rmse
Var.shift_to_start_of_previous_month
Var.apply_landmask
Var.apply_oceanmask
Var.make_lonlat_mask
Base.replace(var::OutputVar, old_new::Pair...)
```

## Leaderboard
```@docs
Leaderboard.RMSEVariable
Leaderboard.RMSEVariable(short_name, model_names::Vector{String})
Leaderboard.RMSEVariable(short_name, model_names::Vector{String}, units::Dict)
Leaderboard.RMSEVariable(short_name, model_names::Vector{String}, category_names::Vector{String}, units::Dict)
Leaderboard.RMSEVariable(short_name::String, model_names::Vector{String}, category_names::Vector{String}, RMSEs, units::Dict)
Leaderboard.RMSEVariable(short_name, model_names::Vector{String}, units::String)
Leaderboard.RMSEVariable(short_name, model_names::Vector{String}, category_names::Vector{String}, units::String)
Leaderboard.RMSEVariable(short_name::String, model_names::Vector{String}, category_names::Vector{String}, RMSEs, units::String)
Leaderboard.model_names
Leaderboard.category_names
Leaderboard.rmse_units
Leaderboard.read_rmses
Base.getindex(rmse_var::RMSEVariable, model_name, category)
Base.getindex(rmse_var::RMSEVariable, model_name::String)
Base.setindex!(rmse_var::RMSEVariable, rmse, model_name, category)
Base.setindex!(rmse_var::RMSEVariable, rmse, model_name::String)
Leaderboard.add_category
Leaderboard.add_model
Leaderboard.add_unit!
Leaderboard.find_best_single_model
Leaderboard.find_worst_single_model
Leaderboard.median
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
Visualize.plot_boxplot!
Visualize.plot_leaderboard!
```

## GeoMakie

```@docs
Visualize.oceanmask
Visualize.landmask
Visualize.contour2D_on_globe!
Visualize.heatmap2D_on_globe!
Visualize.plot_bias_on_globe!
```
