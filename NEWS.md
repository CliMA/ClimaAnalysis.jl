ClimaAnalysis.jl Release Notes
===============================
Main
-------

## Mask aware flatten

With this release, you can pass in a mask to `flatten` to do a mask aware
flatten, where all the values corresponding to zeros in the mask are dropped.
Any value of the `OutputVar` whose coordinates correspond to zeros on the mask are
excluded from the flattened data.

```julia
ocean_mask = ClimaAnalysis.generate_ocean_mask(NaN, 1.0; threshold = 0.5)
flat_masked_var = ClimaAnalysis.flatten(var, mask = ocean_mask)
```

## Accessors for FlatVar and Metadata

You can now access data about the dimensions from `FlatVar` and `Metadata` with
the functions you would use for `OutputVar`s. See the example below.

```julia
# var is a OutputVar
flat_var = ClimaAnalysis.flatten(var)
lon_exists = ClimaAnalysis.has_longitude(flat_var)
lon = ClimaAnalysis.longitudes(flat_var)
lon_exists = ClimaAnalysis.has_longitude(flat_var.metadata)
lon = ClimaAnalysis.longitudes(flat_var.metadata)
```

## Bug fixes
- Fixed conversion from dates to relative times for NetCDF files where the
  temporal dimension is typed as `Union{T, Missing}`, where `T` is a subtype of
  `Dates.AbstractDateTime`, but contains no actual `missing` values.


v0.5.19
-------
This release introduces the following features and bug fixes
- [Selectors](#selectors)
- [Masks](#masks)
- [Read multiple NetCDF files](#read-multiple-netcdf-files)
- [NCCatalog](#nccatalog)

## Selectors

For both `slice` and `window`, a new keyword argument `by` is now available. You
can use this argument to specify how the data is selected.

```julia
# Default is selecting by using the nearest value
ClimaAnalysis.slice(var, time = 2.2, by = ClimaAnalysis.NearestValue())
ClimaAnalysis.window(var, "time", left = 2.2, right = 9.5, by = ClimaAnalysis.NearestValue())

# You can also slice by approximately matching the value
ClimaAnalysis.slice(var, time = 2.0, by = ClimaAnalysis.MatchValue())
ClimaAnalysis.window(var, "time", left = 2.0, right = 9.1, by = ClimaAnalysis.MatchValue())

# You can use an index
ClimaAnalysis.slice(var, time = 2, by = ClimaAnalysis.Index())
ClimaAnalysis.window(var, "time", left = 1, right = 4, by = ClimaAnalysis.Index())
```

## Masks

The function `Var.make_lonlat_mask` is now deprecated in favor of
`Var.generate_lonlat_mask`. Additionally, the land and ocean masks have been improved to
provide better masking near the coastlines. Finally, the `threshold` keyword argument is
added for `Var.apply_landmask` and `Var.apply_oceanmask`.

```julia
# var is a OutputVar whose data consists of zeros and ones
# In this example. all zeros are replaced with NaNs and all ones are replaced with 1.0 when
# applying the mask
mask_fn = generate_lonlat_mask(var, NaN, 1.0, threshold = 0.99)

# Apply mask
masked_var = mask_fn(another_var)
```

## Read multiple NetCDF files

Reading multiple NetCDF files along the time dimension is now supported for the `OutputVar`
constructor. For example, you can now do

```julia
pr_var = ClimaAnalysis.OutputVar(["pr1.nc", "pr2.nc"])
```

to load the NetCDF files "pr1.nc" and "pr2.nc" along the time dimension in a
`OutputVar`.

## NCCatalog

You can now organize and easily load external NetCDF files with `ClimaAnalysis.NCCatalog`.
See the example below of adding a file and loading a variable from the file as a
`OutputVar`.

```julia
catalog = ClimaAnalysis.NCCatalog()
# The file "precip.nc" contains a variable named "precip", but you can retrieve it as "pr".
ClimaAnalysis.add_file!(catalog, "precip.nc", "precip" => "pr")
# The short name of pr_var is "pr".
pr_var = get(catalog, "pr")
```

v0.5.18
-------
This release introduces the following features and bug fixes
- [Concatenation of OutputVars](#concatenation-of-outputvars)
- [FlatVar](#flatvar)
- [Seasonal average](#seasonal-average)
- [New keyword argument for `split_by_season`](#new-keyword-argument-for-split_by_season)
- [TemplateVar](#templatevar)
- [Bug fixes](#bug-fixes)

## Concatenation of OutputVars

You can concatenate `OutputVar`s along a single dimension. This function is helpful if you
need to concatenate `OutputVar`s after applying `split_by_season_across_time` or `window`.

```julia
seasons = ClimaAnalysis.split_by_season_across_time(var);
DJF = cat(seasons[begin:4:end]..., dim = "time");
```

## FlatVar

You can now flatten `data` of `OutputVar`s using `ClimaAnalysis.flatten` and `ClimaAnalysis.unflatten`.
Flattening data gives `FlatVar` whose flattened data can be accessed with `flat_var.data` and the metadata
can be accessed with `flat_var.metadata`. Then, unflattening data reconstructs the original `OutputVar`.
This feature is helpful with preparing data for use in other pipelines.

```julia
flat_var = ClimaAnalysis.flatten(var)

# Flattened data and metadata can be accessed like this
flat_var.data
flat_var.metadata

# Reconstructing the original OutputVar
unflattened_var = ClimaAnalysis.unflatten(flat_var)
unflattened_var = ClimaAnalysis.unflatten(flat_var.metadata, flat_var.data)
```

## Seasonal average

You can now compute seasonal averages with `average_season_across_time`. This returns
a `OutputVar`, where the time dimension consists of the first date of each season, and the
years and seasons can be accessed by `var.attributes["season"]` and `var.attributes["year"]`.

```julia
seasonal_averages_var = average_season_across_time(var)

# Return a vector of seasons
seasonal_averages_var.attributes["season"]

# Return a vector of years
seasonal_averages_var.attributes["year"]
```

## New keyword argument for `split_by_season`

You can now specify the order of the seasons for `split_by_season`. See the
example below of specifying the order of the seasons being returned. This is
helpful if you only want one season, or you want the seasons in a particular
order.

```julia
DJF, SON = split_by_season(var, seasons = ("DJF", "SON"))
```

## TemplateVar

For testing and developing ClimaAnalysis, this release adds `TemplateVar`, which easily
allows the creation of `OutputVar`. See the example below.

```julia
var =
    Template.TemplateVar() |>
    Template.add_attribs(long_name = "Test", short_name = "test") |>
    Template.add_time_dim(name = "t", dim = collect(0.0:2.0)) |>
    Template.add_lon_dim(units = "degrees") |>
    Template.add_lat_dim(axis = "Y") |>
    Template.one_to_n_data(collected = true) |>
    Template.initialize
```

For more information, see the
[documentation](https://clima.github.io/ClimaAnalysis.jl/dev/developer/).

## Bug fixes

- Fixed support for reductions when dimensions have only one point.
- Exclude NaNStatistics v0.6.51 and v0.6.52 in compat, because of correctness issue with
`nanmean`. See this [issue](https://github.com/brenhinkeller/NaNStatistics.jl/issues/59) in
NaNStatistics.jl.

v0.5.17
-------

## Bug fixes

- If the points of a grid represent center of a cell instead of the edges of a cell, then
resampling may produce incorrect results when interpolating outside of the longitude range.
This is fixed in this release.

## Reconstructing dates

`ClimaAnalysis.dates` will now try returning a list of `Dates` even when `dates` is not
among the available dimensions. For this to work, `var.attributes` have to contain
`start_date`. The units of time are assumed to be seconds.

## Support for more binary and unary operations

`OutputVar`s can now be directly manipulated with additional operations,
including `sqrt`, `min`, `log`, `-`.

## Variadic `get` for `SimDir`

`get` can now be used to fetch multiple variables at the same time. For example
```julia
ta, va, ua = get(simdir, "ta", "va", "ua")
```
This works only if the variables are fully identified by their short name.

## More flexible `slice` and `window`

Now, `slice` and `window` will try to match the name of the dimension with a
database of known names. For example, they understand `lat` even if the data has
dimensions `latitude`.

`slice` and `window` can also now work directly with `Dates.DateTime`s (assuming
that the `OutputVar` contains enough information to convert from `times` to
`dates`, namely the `start_date` among the attributes).

Example
```julia
# Suppose `var` is an `OutputVar`
julia> var.attributes["start_date"]
2008-01-01T00:00:00

julia> time_name(var)
"t"

# Notice how we call `time` and not `t`. This was not previously possible!
julia> slice(var, time = Dates.DateTime(2009, 12, 15))
```

## Split by month

Similar to `split_by_season` and `split_by_season_across_time`, the function
`split_by_month` splits a `OutputVar` by month and return a vector of `OutputVar`s where
each `OutputVar` corresponds to a single month.

## More Makie types for plotting
Before this release, plotting was limited to `Makie.Figure` and `Makie.Layout`. With this
release, plotting is extended to `Makie.GridPosition` and `Makie.GridSubposition`. See the
example below.

```
var = get(sim_dir, "ta")
fig = Makie.Figure()
ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, 1], var)
ClimaAnalysis.Visualize.heatmap2D_on_globe!(fig[1, 2][1, 1], var)
```

## Computing variance

Similar to `average_lon`, `average_lat`, and `average_time`, there are now
`variance_lon`, `variance_lat`, and `variance_time`, which compute the variance
along the longitude, latitude, or time dimension respectively.

## Permute dimensions

In addition to `reordered_as`, there is also `permutedims(var, perm)` which permute the
dimensions of an `OutputVar` according to `perm`, an iterable of dimension names.

## Resample with keyword arguments of dimension name and array

You can now resample without needing another `OutputVar` to resample on. Note
that checking for units is not done. See the example below.

```julia
resampled_var = resampled_as(src_var, long = [0.0, 1.0], lat = [2.0, 3.0])
```


v0.5.16
-------

## Shift longitudes
The implementation of `center_longitude!` is incorrect and buggy. To resolve this, the
function `shift_longitude` is encouraged to be used instead. For example, to shift
longitudes from 0 to 360 degrees to -180 to 180 degrees, one can use
`shift_longitude(var, -180.0, 180.0)` and to shift longitudes from -180 to 180 degrees to 0
to 360 degrees, one can use `shift_longitude(var, 0.0, 360.0)`.

v0.5.15
-------

## Fix Var.weighted_lonlat
`Var.weighted_average_lonlat` can throw an error depending on the order of the dimensions.
This has been fixed in this release.

v0.5.14
-------
This release introduces the following features and bug fixes
- [Split by seasons across time](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#split-by-seasons-across-time)
- [Partial resampling](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#partial-resampling)
- [Add global average along longitude and latitude](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#add-global-average-along-longitude-and-latitude)
- [Show method defined for SimDir](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#show-method-defined-for-simdir)
- [In-place functions](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#in-place-functions)
- [Retain attributes when doing binary operations between `OutputVar`s](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#remove-interpolant-from-outputvar)

## Split by seasons across time
It may be the case that you want to split a `OutputVar` by season, while keeping each year
separate. This is different from `split_by_season`, which ignores that seasons can come from
different years. This can be done by using `split_by_season_across_time`. For example, if a
`OutputVar` contains times corresponding to 2010-01-01, 2010-03-01, 2010-06-01, 2010-09-01,
and 2010-12-01, then the result of `split_by_season_across_time` is five `OutputVar`s, each
corresponding to a distinct date. Even though 2010-01-01 and 2010-12-01 are in the same
season, there are two `OutputVar`s, because the dates do not belong in the same season and
year.

## Partial resampling
Sometimes, it is necessary to resample over a subset of all the dimensions of a variable
(e.g., resample the spatial dimensions while preserving the temporal ones). `resampled_as`
now includes a new keyword argument to allow this. When called with `dim_names`,
`resampled_as` will resample only the specified dimensions. For example, if one wants to
resample only over longitude and latitude, then one can use
`resampled_as(src_var, dest_var, dim_names = ["longitude, "latitude"])`.

## Add global average along longitude and latitude
Before, you could only compute an average over one dimension. With this release, you can use
`average_lonlat` to compute the average over the longitude and latitude dimensions and
`weighted_average_lonlat` to compute the global latitude-weighted average over the longitude
and latitude dimensions. Both functions ignore `NaN` by default.

## Show method defined for SimDir
Previously, a brief description of `SimDir` can be printed to the terminal using
`summary(simdir)`. With this release, the same information can also be viewed by using
`show(simdir)`.

## In-place functions
With this release, there is an in-place version of `reverse_dim` which is `reverse_dim!`.

## Retain attributes when doing binary operations between `OutputVar`s
With this release, if two `OutputVar`s share the same start date, then the start date will
remain in the resulting `OutputVar`s after performing binary operations on them.

v0.5.13
-------

## Remove interpolant from OutputVar
With this release, the field `interpolant` is removed from `OutputVar`. This was done to
prevent the largely unnecessary construction of an interpolant every time a `OutputVar` is
constructed. With this change, a `OutputVar` take up about 50% less memory.

However, functions like `resampled_as` and interpolating using a `OutputVar` will be slower
as an interpolant must be generated. This means repeated calls to these functions will be
slower compared to the previous versions of ClimaAnalysis.

## Better error messages
There is now error hints when using a function that requires another package such as Makie
or GeoMakie to be loaded as well. The error hint tells the user which package need to be
loaded in, so that the function can be used.

## Remake OutputVars
With this release, you can remake a `OutputVar` using an already existing `OutputVar`. This
is helpful if you need to construct a new `OutputVar` from an already existing one, but only
need to modify one field while leaving the other fields the same.

## Reordering categories
There is now support for reordering the categories in a `RMSEVariable`. This is helpful
when making plots using `Visualize.plot_boxplot!` and `Visualize.plot_leaderboard!`.

## Pretty print
`OutputVar`s now display information about the attributes, dimension attributes, and the
dimensions that the data is defined over in the terminal.

v0.5.12
-------

## Features
This release introduces the following features and bug fixes
- [Compute global root mean squared error for 3D `OutputVar`s](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#compute-global-root-mean-squared-error-for-3d-outputvars)
- [Reversing a dimension](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#reversing-a-dimension)
- [Converting the units of a dimension](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#converting-the-units-of-a-dimension)
- [Bug fixes](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#bug-fixes)

### Compute global root mean squared error for 3D `OutputVar`s
For 3D variables defined over longitude, latitude, and pressure, one can find the global
RMSE in pressure space using `ClimaAnalysis.global_rmse_pfull`. See an example of this
below, where global RMSE is computed between 3D variables in pressure space.

```julia
# Load in 3D temperature variable defined over longitude, latitude, and pressure
julia> obs_var = OutputVar("era5_pfull_ta_data.nc"); # load in observational data

# Load in 3D temperature variable defined over longitude, latitude, and z
julia> sim_var = get(simdir("simulation_output"), "ta"); # load in simulation data

# Load in 3D pressure variable defined over longitude, latitude, and z
julia> pressure_3D_var = get(simdir("simulation_output"), "pfull"); # load in simulation data

# This function will automatically converts to pressure coordinates for `sim_var` and
# `obs_var` if the keywords `sim_pressure` and `obs_pressure` are supplied respectively
julia> ClimaAnalysis.global_rmse_pfull(sim_var, obs_var, sim_pressure = pressure_3D_var)
3.4
```

### Reversing a dimension
Before, a `OutputVar` will not be created, if the dimensions are not in increasing order.
This is because an interpolant cannot be made in this case. With this release, the condition
is relaxed, so that a `OutputVar` can be made, but an interpolant is not made. The function
`Var.reverse_dim` is provided to reverse the order of a dimension by name. See the example
below.

```julia
# Reversing pressure levels so that an interpolant can be made
var_reversed = ClimaAnalysis.reverse_dim(var, "pressure_level")
```

### Converting the units of a dimension
Similar to converting the units of data, you can convert the units of a dimension with
`ClimaAnalysis.convert_dim_units`. See the example below.
```julia
new_var = ClimaAnalysis.convert_dim_units(
        var,
        "lat",
        "rads",
        conversion_function = x -> x * π / 180.0,
    )
```

## Bug fixes
- `Atmos.to_pressure_coordinates` now works with Unitful units.
- `Atmos.to_pressure_coordinates` now uses reasonable pressure values when `target_pressure`
  is not specified. In particular, the vertical dimension is mapped to pressure levels by z
  -> P0 * exp(-z / H_EARTH), where P0 = 10000 and H_EARTH = 7000.0, following a simple
  hydrostatic model for the atmosphere.
- `Var.dim_units` returns a string even if the type of the units is Unitful.
- Add periodic boundary condition when extrapolating on longitude dimension when it is
  exactly 360 degrees.
- Update land-sea mask to use the land-sea mask [here](https://github.com/CliMA/ClimaArtifacts/tree/main/landsea_mask).

v0.5.11
-------

## Features

### Applying a land/sea mask to `OutputVar`s

There is now support to applying a land or sea mask to `OutputVar`s though
`apply_landmask(var)` and `apply_oceanmask(var)` respectively. A land mask sets all the
coordinates corresponding to land to zeros in the data of the `OutputVar` and a sea mask
sets all the coordinates corresponding to ocean to zeros in the data of the `OutputVar`.
Furthermore, the parameter `mask` is added to the functions `bias`, `global_bias`,
`squared_error`,`global_mse`, and `global_rmse` which takes either `apply_landmask` and
`apply_oceanmask`. This is useful for computing these quantities only over the land or sea.

```julia
var_no_land = ClimaAnalysis.apply_landmask(var)
var_no_ocean = ClimaAnalysis.apply_oceanmask(var)
```

### Aggregating data over the time dimension using SimDir

Before, `SimDir` only search for `.nc` files in `simulation_path`. This update now allows
for `SimDir` to recursively search for `.nc` files in `simulation_path`. If the same
combination of short name, reduction, and period is found, then calling `get` on `SimDir`
aggregate the NetCDF files along the time dimension. This is useful if you need to restart
the simulation multiple times and the data is scattered across different `.nc` files.

See the example below of this functionality. The file directories might look like this:
```
global_diagnostics/
├── output_0001
│   ├── pfull_2.0d_inst.nc
│   ├── ts_1.0h_max.nc
│   └── ts_10d_average.nc
├── output_0002
│   ├── ts_1.0h_max.nc
│   └── ts_10d_average.nc
└── output_0003
    └── pfull_2.0d_inst.nc
```
The following code would work:
```julia
ts_max_var = ClimaAnalysis.get(simdir, short_name = "ts", reduction = "max", period = "1.0h")
pfull_var = ClimaAnalysis.get(simdir, short_name = "pfull", reduction = "2.0d", period = "inst")
```

### Replace values in data of a `OutputVar`
When dealing with land or ocean data, there can potentially be `missing` or `NaN` values in
the data. The function `replace` can be used to replace `missing` or `NaN` values in
`Var.data` with another value like 0.0. See the example below of this usage.
```julia
ClimaAnalysis.replace(var, NaN => 0.0, missing => 0.0)
```

### Set units for dimensions
Similar to `set_units`, there is the function `set_dim_units!` which one can use to set
the units of a dimension.

```julia
new_var = ClimaAnalysis.set_dim_units!(var, "lon", "degrees_east")
```

### Automatically generate a mask from a OutputVar
Masking function can automatically be generated from the function `make_lonlat_mask`. See
the example below of generating a masking function that mask out any data that is `NaN` in
`var`.

```julia
# Replace NaN with 0.0 and everything else with 1.0 for the mask
mask_fn = ClimaAnalysis.make_lonlat_mask(var; set_to_val = isnan, true_val = 0.0, false_val = 1.0)
another_masked_var = mask_fn(another_var)
```

### Using masking function when plotting
Masking functions can now be passed for the `mask` keyword for plotting functions. See the
example below of plotting with a masking function.

```julia
import ClimaAnalysis
import ClimaAnalysis.Visualize: plot_bias_on_globe!, oceanmask
import GeoMakie
import CairoMakie

mask_var = ClimaAnalysis.OutputVar("ocean_mask.nc")
mask_fn = ClimaAnalysis.make_lonlat_mask(mask_var; set_to_val = isnan)

obs_var = ClimaAnalysis.OutputVar("ta_1d_average.nc")
sim_var = ClimaAnalysis.get(ClimaAnalysis.simdir("simulation_output"), "ta")

fig = CairoMakie.Figure()
plot_bias_on_globe!(fig, var, mask = mask_fn)
CairoMakie.save("myfigure.pdf", fig)
```

## Bug fixes
- Masking now affects the colorbar.
- `Var.shift_to_start_of_previous_month` now checks for duplicate dates and throws an error
if duplicate dates are detected.
- Fix issue with plotting multiple figures at the same time.
- Improve shading for `Visualize.heatmap2D_on_globe!`.
- Add support for automatically converting CFTime.AbstractCFDateTime dates to seconds.
- Treat `NaN`s as zeros when integrating (`integrate_lon, integrate_lat, integrate_lonlat`).
- Fix duplicate values on colorbar for plotting bias.
- Fix `NaN`s for colorbar when plotting bias.

v0.5.10
-------

## Features

### Set units

You can now set units for a `OutputVar`. This is useful if you need to change the name of
the units or units are missing.

```julia
new_var = ClimaAnalysis.set_units(var, "kg m s^-1")
```

### Extrapolating `OutputVar` on longitude and latitude
Extrapolation is now possible for the longitude and latitude dimensions. If the dimension
arrays are equispaced and span the entire range, then a periodic boundary condition is added
for the longitude dimension and a flat boundary condition is added for the latitude
dimension.

### Preprocess dates and times
There is now support for preprocessing dates and times. The constructor for reading NetCDF
files now automatically converts dates to seconds in the time dimension. This is done
because `ClimaAnalysis` does not support interpolating on dates which mean functions that
rely on the interpolats, such as `resampled_as`, will not work otherwise.

Also, the constructor supports two additional parameters `new_start_date` and `shift_by`.
After converting from dates to seconds, the seconds are shifted to match `new_start_date`.
If preprocessing of dates is needed before shifting to `new_start_date`, then the parameter
`shift_by` can be used as it accepts a function that takes in `Dates.DateTime` elements and
return Dates.DateTime elements. This function is applied to each element of the time array.
```julia
# Shift the dates to first day of month, convert to seconds, and adjust seconds to
# match "1/1/2010"
shift_var = OutputVar(
        "test.nc",
        "pr",
        new_start_date = "1/1/2010", # or Dates.DateTime(2010, 1, 1)
        shift_by = Dates.firstdayofmonth,
    )
```

The function `Var.shift_to_start_of_previous_month` is added to shift the times in the time
dimension to the end of the previous month. This function is helpful in ensuring consistency
in dates between simulation and observational data. One example of this is when adjusting
monthly averaged data. For instance, data on 2010-02-01 in the `OutputVar` corresponds to
the monthly average for January. This function shifts the times so that 2010-01-01 will
correspond to the monthly average for January.
```julia
sim_var = shift_to_start_of_previous_month(sim_var)
```

## Bug fixes

- Interpolation is not possible with dates. When dates are detected in any dimension, an
  interpolat will not be made.
- Fix identifying variables with underscore in the short name (such as
  `net_toa_flux`). ([#109](https://github.com/CliMA/ClimaAnalysis.jl/pull/109
  "PR109"))

## Minor changes

- Var.arecompatible only check for the units of the dimensions instead of checking that the
  dimension attributes fully match.

v0.5.9
------

## Features

### Reorder dimensions
Before, resampling requires that the order of the dimensions is the same between the two
`OutputVar`s. This feature adds the functionality of reordering the dimensions in a
`OutputVar` to match the ordering of another `OutputVar`. The function `resampled_as` is
updated to use `reordered_as`. See the example below of this functionality.

```julia
julia> src_var.dims |> keys |> collect
2-element Vector{String}:
 "long"
 "lat"

julia> dest_var.dims |> keys |> collect
2-element Vector{String}:
 "lat"
 "long"

julia> reordered_var = ClimaAnalysis.reordered_as(src_var, dest_var);

julia> reordered_var.dims |> keys |> collect
2-element Vector{String}:
 "lat"
 "long"
```

## Bug fixes

- Fix models repeating in legend of box plots by not considering the models in `model_names`
  when finding the best and worst models.
- Fix legend from covering the box plot by adding the parameter `legend_text_width` which
  control the number of characters on each line of the legend of the box plot.
- Use default marker size instead of a marker size of 20 when plotting other models beside
  `CliMA` on the box plot.
- Fix support for `""` in units.

v0.5.8
------

## Features
This release introduces the following features
- [Directly reading NetCDF files](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#add-support-for-directly-reading-netcdf-files)
- [Resampling a OutputVar using the dimensions from another OutputVar](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#resampling-a-outputvar-using-the-dimensions-from-another-outputvar)
- [Add support for converting units](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#add-support-for-converting-units)
- [Applying a land/sea mask to GeoMakie plots](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#applying-a-landsea-mask-to-geomakie-plots)
- [Integrating OutputVar with respect to longitude or latitude](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#integrating-outputvar-with-respect-to-longitude-or-latitude)
- [Splitting OutputVar by season](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#splitting-outputvar-by-season)
- [Compute bias and squared error between OutputVar](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#compute-bias-and-squared-error-between-outputvar)
- [Represent RMSEs from various models](https://github.com/CliMA/ClimaAnalysis.jl/blob/main/NEWS.md#represent-rmses-from-various-models)

### Add support for directly reading NetCDF files

Starting version 0.5.8, `ClimaAnalysis` now supports NetCDF files that were not
generated with
[`ClimaDiagnostics`](https://github.com/CliMA/ClimaDiagnostics.jl) [0].

To load a `NetCDF` file into a `ClimaAnalysis.OutputVar`, just pass the path of
such file to the constructor
```julia
import ClimaAnalysis: OutputVar

myfile = OutputVar("my_netcdf_file.nc")
```
`ClimaAnalysis` will try to find a variable in this file. If multiple are available,
`ClimaAnalysis` picks the latest in alphabetical order. If you want to specify one,
pass it to the constructor:
```julia
import ClimaAnalysis: OutputVar

myfile = OutputVar("my_netcdf_file.nc", "myvar")
```

[0] Restrictions apply.

### Resampling a `OutputVar` using the dimensions from another `OutputVar`

You can use the `resampled_as(src_var, dest_var)` function where `src_var` is a
OutputVar with the data you want to resample using the dimensions in another
OutputVar `dest_var`. If resampling is possible, then a new `OutputVar` is
returned where the data in `src_var` is resampled using a linear interpolation
to fit the dimensions in `dest_var`. Resampling is not possible when the
dimensions in either `OutputVar`s are missing units, the dimensions between the
`OutputVar`s do not agree, or the data in `src_var` is not defined everywhere on
the dimensions in `dest_var`.

```julia
julia> src_var.data
3×4 Matrix{Float64}:
 1.0  4.0  7.0  10.0
 2.0  5.0  8.0  11.0
 3.0  6.0  9.0  12.0

julia> src_var.dims
OrderedDict{String, Vector{Float64}} with 2 entries:
  "lon"      => [0.0, 1.0, 2.0]
  "latitude" => [0.0, 1.0, 2.0, 3.0]

julia> dest_var.dims # dims that src_var.data should be resampled on
OrderedDict{String, Vector{Float64}} with 2 entries:
  "long" => [0.0, 1.0]
  "lat"  => [0.0, 1.0, 2.0]

julia> resampled_var = ClimaAnalysis.resampled_as(src_var, dest_var);

julia> resampled_var.data
2×3 Matrix{Float64}:
 1.0  4.0  7.0
 2.0  5.0  8.0

julia> resampled_var.dims # updated dims that are the same as the dims in dest_var
OrderedDict{String, Vector{Float64}} with 2 entries:
  "lon"      => [0.0, 1.0]
  "latitude" => [0.0, 1.0, 2.0]
```

### Add support for converting units

`ClimaAnalysis` now uses
[Unitful](https://painterqubits.github.io/Unitful.jl/stable) to handle variable
units, when possible.

When a `OutputVar` has `units` among its `attributes`, `ClimaAnalysis` will try
to use `Unitful` to parse it. If successful, `OutputVar` can be directly
converted to other compatible units. For example, if `var` has units of `m/s`,
```julia-repl
julia> ClimaAnalysis.convert_units(var, "cm/s")
```
will convert to `cm/s`.

Some units are not recognized by `Unitful`. Please, open an issue about that:
we can add more units.

In those cases, or when units are incompatible, you can also pass a
`conversion_function` that specify how to transform units.
```julia-repl
julia> ClimaAnalysis.convert_units(var, "kg/s", conversion_function = (x) - 1000x)
```

### Applying a land/sea mask to `GeoMakie` plots

When plotting with `GeoMakie` (ie, using the `contour2D_on_globe!` and
`heatmap2D_on_globe!` function), it is now possible to mask out a portion of the
output. The most common use cases are to hide the ocean or the continents.

To hide the ocean, you can now pass `mask = ClimaAnalysis.Visualize.oceanmask()`
to the globe plotting functions. You can customize how the mask is plotted by
passing the `:mask` extra keywords. For example:
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

### Integrating `OutputVar` with respect to longitude or latitude

You can use the `integrate_lon(var)`, `integrate_lat(var)`, or `integrate_lonlat(var)`
functions for integrating along longitude, latitude, or both respectively. The bounds of
integration are determined by the range of the dimensions longitude and latitude in `var`.
The unit of both longitude and latitude should be degree.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints for integrating
longitude and latitude. See the example of integrating over a sphere where the data is all
ones to find the surface area of a sphere.

```julia
julia> lon = collect(range(-179.5, 179.5, 360));

julia> lat = collect(range(-89.5, 89.5, 180));

julia> data = ones(length(lon), length(lat));

julia> dims = OrderedDict(["lon" => lon, "lat" => lat]);

julia> dim_attribs = OrderedDict([
           "lon" => Dict("units" => "degrees_east"),
           "lat" => Dict("units" => "degrees_north"),
       ]);

julia> attribs = Dict("long_name" => "f");

julia> var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data);

julia> integrated_var = integrate_lonlat(var);

julia> integrated_var.dims # no dimensions since longitude and latitude are integrated out
OrderedDict{String, Vector{Float64}}()

julia> integrated_var.data # approximately 4π (the surface area of a sphere)
0-dimensional Array{Float64, 0}:
12.566530113084296

julia> long_name(integrated_var) # updated long name to reflect the data being integrated
"f integrated over lon (-179.5 to 179.5degrees_east) and integrated over lat (-89.5 to 89.5degrees_north)"
```

### Splitting OutputVar by season
`OutputVar`s can be split by seasons using `split_by_season(var)` provided that a start date
can be found in `var.attributes["start_date"]` and time is a dimension in the `OutputVar`.
The unit of time is expected to be second. The function `split_by_season(var)` returns a
vector of four `OutputVar`s with each `OutputVar` corresponding to a season. The months of
the seasons are March to May, June to August, September to November, and December to
February. The order of the vector is MAM, JJA, SON, and DJF. If there are no dates found for
a season, then the `OutputVar` for that season will be an empty `OutputVar`.

```julia
julia> attribs = Dict("start_date" => "2024-1-1");

julia> time = [0.0, 5_184_000.0, 13_132_800.0]; # correspond to dates 2024-1-1, 2024-3-1, 2024-6-1

julia> dims = OrderedDict(["time" => time]);

julia> dim_attribs = OrderedDict(["time" => Dict("units" => "s")]); # unit is second

julia> data = [1.0, 2.0, 3.0];

julia> var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data);

julia> MAM, JJA, SON, DJF = ClimaAnalysis.split_by_season(var);

julia> isempty(SON) # empty OutputVar because no dates between September to November
true

julia> [MAM.dims["time"], JJA.dims["time"], DJF.dims["time"]]
3-element Vector{Vector{Float64}}:
 [5.184e6]
 [1.31328e7]
 [0.0]

julia> [MAM.data, JJA.data, DJF.data]
3-element Vector{Vector{Float64}}:
 [2.0]
 [3.0]
 [1.0]
```

### Compute bias and squared error between OutputVar
Bias and squared error can be computed from simulation data and observational data in
`OutputVar`s using `bias(sim, obs)` and `squared_error(sim, obs)`. The function `bias(sim,
obs)` returns a `OutputVar` whose data is the bias (`sim.data - obs.data`) and computes the
 global bias of `data` in `sim` and `obs` over longitude and latitude. The result is stored
in `var.attributes["global_bias"]`. The function `squared_error(sim, obs)` returns a
`OutputVar` whose data is the squared error (`(sim.data - obs.data)^2`) and computes the
global mean squared error (MSE) and the global root mean squared error (RMSE) of `data` in
`sim` and `obs` over longitude and latitude. The result is stored in
`var.attributes["global_mse"]` and `var.attributes["global_rmse"]`. Resampling is
automatically done by resampling `obs` on `sim`. If you are only interested in computing
global bias, MSE, or RMSE, you can use `global_bias(sim, obs)`, `global_mse(sim, obs)`, or
`global_rmse(sim, obs)`.

As of now, these functions are implemented for `OutputVar`s with only the dimensions
longitude and latitude. Furthermore, units must be supplied for data and dimensions in `sim`
and `obs` and the units for longitude and latitude should be degrees.

Consider the following example, where we compute the bias and RMSE between our simulation
and some observations stored in "ta\_1d\_average.nc".

```julia
julia> obs_var = OutputVar("ta_1d_average.nc"); # load in observational data

julia> sim_var = get(simdir("simulation_output"), "ta"); # load in simulation data

julia> ClimaAnalysis.short_name(sim_var)
"ta"

julia> bias_var = ClimaAnalysis.bias(sim_var, obs_var); # bias_var is a OutputVar that can be plotted

julia> global_bias(sim, obs)
2.0

julia> units(bias_var)
"K"

julia> se_var = ClimaAnalysis.squared_error(sim_var, obs_var); # can also be plotted

julia> global_mse(sim, obs)
4.0

julia> global_rmse(sim, obs)
2.0

julia> units(se_var)
"K^2"
```

### Plotting bias

Building upon the other features introduced in this release, you can now directly plot bias
and root mean squared error between two variables with the `plot_bias_on_globe!` function.
Typically, this is done to compare simulated data against observations.

In the example below, we plot the bias between our simulation and some observations stored
in `ta_1d_average.nc`.

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

### Represent RMSEs from various models

To facilitate analysis of root mean squared errors (RMSEs) over different models and
categories (e.g., seasons) for a single variable of interest, `RMSEVariable` is introduced in
this release. See the examples below of constructing a `RMSEVariable` using a short name, a
vector of model names, a vector of categories, and a dictionary mapping model names to units
or a string of the name of the unit.

```julia
import ClimaAnalysis

rmse_var = ClimaAnalysis.RMSEVariable("ta", ["ACCESS-CM2", "ACCESS-ESM1-5"])
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"),
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"),
)
# Convenience functions if models all share the same unit
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    "K",
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    "K",
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    ones(2, 5),
    "K",
)
```

A `RMSEVariable` can be inspected using `model_names`, `category_names`, and `rmse_units`
which provide the model names, the category names, and the units respectively.

#### Reading RMSEs from CSV file

A CSV file containing model names in the first column and root mean squared errors in the
subsequent columns with a header describing each category (i.e. seasons) can be read into
a `RMSEVariable`. See the example below on how to use this functionality.

```julia
rmse_var = ClimaAnalysis.read_rmses("./data/test_csv.csv", "ta")
rmse_var = ClimaAnalysis.read_rmses(
    "./data/test_csv.csv",
    "ta",
    units = Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"), # passing units as a dictionary
)
rmse_var = ClimaAnalysis.read_rmses(
    "./data/test_csv.csv",
    "ta",
    units = "K", # passing units as a string
)
```

#### Indexing
`RMSEVariable` supports indexing by integer or string. See the example for indexing into
a `RMSEVariable`.

```julia
rmse_var["ACCESS-CM2"]
rmse_var[:, "MAM"]
rmse_var["ACCESS-CM2", ["ANN", "DJF", "MAM"]]
rmse_var[2,5] = 11.2
rmse_var[:, :]
```

#### Adding categories, models, and units
Adding categories (e.g., seasons, months, years, etc.), models, and units to a `RMSEVariable`
can be done using `add_category`, `add_model`, and `add_units!`.

See the example below for how to use this functionality.

```julia
rmse_var2 = ClimaAnalysis.add_category(rmse_var, "Jan") # can take in mode than one category
rmse_var = ClimaAnalysis.add_model(rmse_var, "CliMA") # can take in more than one model name
ClimaAnalysis.add_unit!(rmse_var, "CliMA", "K")
ClimaAnalysis.add_unit!(rmse_var, Dict("CliMA" => "K")) # for adding multiple units
```

#### Summary statistics
Comparsion between models can be done using `find_best_single_model`,
`find_worst_single_model`, and `median`. The functions `find_best_single_model` and
`find_worst_single_model` default to the category "ANN" (corresponding to the annual mean),
but any category be considered using the parameter `category_name`. Furthermore, the model's
root mean squared errors (RMSEs) and name is returned. The function `median` only return the
model's RMSEs. Any `NaN` that appear in the data is ignored when computing the summary
statistics. See the example below on how to use this functionality.

```julia
ClimaAnalysis.find_best_single_model(rmse_var, category_name = "DJF")
ClimaAnalysis.find_worst_single_model(rmse_var, category_name = "DJF")
ClimaAnalysis.median(rmse_var)
```

#### Plotting RMSEVariable
`RMSEVariable` can be visualized as a box plot or heat map using `plot_boxplot!` and
`plot_leaderboard!`. The function `plot_boxplot!` makes a box plot for each category in the
`RMSEVariable` and plots any other models as specified by `model_names`. The function
`plot_leaderboard!` makes a heatmap of the RMSEs between the variables of interest and the
categories. The values of the heatmap are normalized by dividing over the median model's
RMSEs for each variable.

## Bug fixes

- Increased the default value for `warp_string` to 72.
- Binary operation between Real and OutputVar now update the interpolation of the resulting
  OutputVar

## New compat requirements

`ClimaAnalysis` 0.5.8 drops support for versions of `GeoMakie` prior to `0.7.3`.
This change is required to acquire land-sea mask data from `GeoMakie`. Version
`0.7.3` specifically is also required because it fixes a precompilation bug in
`GeoMakie`. As a result, the minimum version of `Makie` is now `0.21.5`.

- `GeoMakie` >= 0.7.3
- `Makie` >= 0.21.5
- `CairoMakie` >= 0.12.0

v0.5.7
------
- Add support for evaluating `OutputVar`s onto arbitrary target points (with
  multilinear interpolation).
- `average` operations now ignore `NaN`s by default.
- Add `has_*` methods to query whether a `Var` has a given dimension (e.g., `z`).
- Support `Makie` backends besides `CairoMakie`.
- Add methods to get the range and units of a given dimension in `Var`.

v0.5.6
------
- Fix finding variables with name like `clwup_1m_40s_inst.nc` (composed period).
- Add support for weighted averages in `average_lat`.

v0.5.5
------
- Fix reading `NetCDF` files with dimensions in incorrect order.

v0.5.4
------
- Added support for extraction dimension from functions, such as `times`.
- Reorganized internal modules so that each file is a module.

v0.5.3
------
- Add `Visualize.contour2D_on_globe!` for discrete contours.
