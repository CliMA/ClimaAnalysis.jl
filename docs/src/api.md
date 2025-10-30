# API

```@meta
CurrentModule = ClimaAnalysis
```

## Sim

```@docs
Sim.SimDir
get(simdir::SimDir)
get(simdir::SimDir, short_names...)
available_vars(simdir::SimDir)
Sim.available_reductions
Sim.available_periods
Base.show(io::IO, simdir::SimDir)
```

## Catalog

```@docs
Catalog.NCCatalog
Catalog.NCCatalog()
Catalog.add_file!
get(catalog::NCCatalog, short_name; var_kwargs = ())
Catalog.available_vars(catalog::NCCatalog)
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
Var.remake
Var.slice
Var.average_lat
Var.weighted_average_lat
Var.average_lon
Var.average_lonlat
Var.weighted_average_lonlat
Var.average_x
Var.average_y
Var.average_xy
Var.average_time
Var.variance_lon
Var.variance_lat
Var.variance_time
Var.window
Var.select
Var.view_select
Var.arecompatible
Var.shift_longitude
Var.center_longitude!
Var.times
Var.dates
Var.longitudes
Var.latitudes
Var.altitudes
Var.pressures
Var.time_name
Var.date_name
Var.longitude_name
Var.latitude_name
Var.altitude_name
Var.pressure_name
Var.has_time
Var.has_date
Var.has_longitude
Var.has_latitude
Var.has_altitude
Var.has_pressure
Var.conventional_dim_name
Var.find_corresponding_dim_name
Var.find_corresponding_dim_name_in_var
Var.dim_units
Var.range_dim
Base.permutedims(var::OutputVar, perm)
Var.reordered_as
Var.resampled_as
Var.convert_units
Var.set_units
Var.convert_dim_units
Var.set_dim_units!
Var.set_short_name!
Var.integrate_lonlat
Var.integrate_lat
Var.integrate_lon
Var.split_by_season(var::OutputVar)
Var.split_by_season_across_time(var::OutputVar)
Var.average_season_across_time
Var.split_by_month(var::OutputVar)
Var.bias
Var.global_bias
Var.squared_error
Var.global_mse
Var.global_rmse
Var.shift_to_start_of_previous_month
Var.shift_to_previous_week
Var.shift_to_previous_day
Var.transform_dates
Var.LonLatMask
Var.apply_landmask
Var.apply_oceanmask
Var.generate_lonlat_mask
Var.generate_land_mask
Var.generate_ocean_mask
Var.make_lonlat_mask
Base.replace(var::OutputVar, old_new::Pair...; count::Integer)
Base.replace(new::Union{Function, Type}, var::OutputVar; count::Integer)
Base.replace!(var::OutputVar, old_new::Pair...; count::Integer)
Base.replace!(new::Union{Function, Type}, var::OutputVar; count::Integer)
Base.cat(vars::OutputVar...; dim::String)
Var.reverse_dim
Var.reverse_dim!
Base.show(io::IO, var::OutputVar)
```

## Selectors

```@docs
ClimaAnalysis.AbstractSelector
ClimaAnalysis.NearestValue
ClimaAnalysis.MatchValue
ClimaAnalysis.Index
ClimaAnalysis.get_index(var, dim_name, val, ::NearestValue)
ClimaAnalysis.get_index(var, dim_name, val, ::MatchValue)
ClimaAnalysis.get_index(var, dim_name, val, ::Index)
```

## FlatVar

!!! note "Accessor functions for versions of ClimaAnalysis after v0.5.18"
    For `Metadata` and `FlatVar`, you can use accessor functions (e.g. `has_time` and
    `times`) defined for `OutputVar`s on them.

```@docs
Var.Metadata
Var.FlatVar
Var.flatten
Var.unflatten
Var.flatten_dim_order
Var.flattened_length
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
Leaderboard.reorder_categories
Leaderboard.match_category_order
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
Atmos.global_rmse_pfull
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
