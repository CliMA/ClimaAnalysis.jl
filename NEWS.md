ClimaAnalysis.jl Release Notes
===============================
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
mask_fn = ClimaAnalysis.make_lonlat_mask(var; set_to_zero = isnan)
another_masked_var = mask_fn(another_var)
```

## Bug fixes
- Masking now affects the colorbar.
- `Var.shift_to_start_of_previous_month` now checks for duplicate dates and throws an error
if duplicate dates are detected.
- Fix issue with plotting multiple figures at the same time.
- Improve shading for `Visualize.heatmap2D_on_globe!`.
- Add support for automatically converting CFTime.AbstractCFDateTime dates to seconds.
- Treat `NaN`s as zeros when integrating (`integrate_lon, integrate_lat, integrate_lonlat`)

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

julia> ClimaAnalysis.isempty(SON) # empty OutputVar because no dates between September to November
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
