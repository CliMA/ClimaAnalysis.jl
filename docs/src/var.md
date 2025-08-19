# `OutputVar`s

`OutputVar`s are the heart of `ClimaAnalysis`. This page is under construction,
in the meantime, consult [`OutputVar`](@ref).

`OutputVar`s can be directly generated from most NetCDF files. Just pass the
path to the constructor:
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

## Physical units

`OutputVar`s can contain information about their physical units. For
`OutputVar`s read from NetCDF files, this is obtained from the `units` attribute
(and stored in the `attributes["units"]`).

When possible, `ClimaAnalysis` uses
[Unitful](https://painterqubits.github.io/Unitful.jl/stable) to handle units.
This enables automatic unit conversion for `OutputVar`s.

Consider the following example:
```julia
import ClimaAnalysis
values = 0:100.0 |> collect
data = copy(values)
attribs = Dict("long_name" => "speed", "units" => "m/s")
dim_attribs = Dict{String, Any}()
var = ClimaAnalysis.OutputVar(attribs, Dict("distance" => values), dim_attribs, data)

var_cms = ClimaAnalysis.convert_units(var, "cm/s")
```
In this example, we set up` var`, an `OutputVar` with units of meters per second.
Then, we called [`ClimaAnalysis.convert_units`](@ref) to convert the units to
centimeters per second.

Sometimes, this automatic unit conversion is not possible (e.g., when you want
to transform between incompatible units). In this case, you can pass a function
that specify how to apply this transformation. For example, in the previous
case, we can assume that we are talking about water and transform units into a
mass flux:
```julia
new_var = ClimaAnalysis.convert_units(var, "kg m/s", conversion_function = (x) -> 1000x)
```

!!! note "Unparseable units"
    If you find some unparseable units, please open an issue. We can fix them!


If units do not exist, or you want to change the name of the units, then one can use the
`set_units` function.
```julia
new_var = ClimaAnalysis.set_units(var, "kg m s^-1")
```

For converting the units of a dimension, you can use
[`ClimaAnalysis.convert_dim_units`](@ref). As of now, automatic conversion is not supported
which means you need to supply the conversion function. See the example below.
```julia
new_var = ClimaAnalysis.convert_dim_units(
        var,
        "lat",
        "rads",
        conversion_function = x -> x * π / 180.0,
    )
```

Similarly, to set the units of a dimension, you can use the `dim_set_units!` function.
```julia
new_var = ClimaAnalysis.set_dim_units!(var, "lon", "degrees_east")
```

!!! warning "Override existing units"
    If units already exist, this will override the units for data or the dimension in `var`.



## Interpolations and extrapolations

Interpolating a `OutputVar` onto coordinates can be done by doing the following:
```julia
var((0.0, 0.0)) # var is a two-dimensional OutputVar
```

A multilinear interpolation is used to determine the value at the coordinate (0, 0).
!!! warning "Interpolate on dates"
    If any of the dimensions contains `Dates.DateTime` elements, interpolation is not
    possible. `Interpolations.jl` does not support interpolating on dates.

Extrapolating is supported only on the longitude and latitude dimensions. For the longitude
and latitude dimensions, a periodic boundary condition and a flat boundary condition are
added, respectively, when the dimension array is equispaced and spans the entire range. For
all other cases, extrapolating beyond the domain of the dimension will throw an error.

## Preprocess dates and seconds
When loading a NetCDF file, dates in the time dimension are automatically converted to
seconds and a start date is added to the attributes of the `OutputVar`. This is done because
`ClimaAnalysis` does not support interpolating on dates which mean functions that rely on
the interpolats, such as `resampled_as`, will not work otherwise.

Two additional parameters are provided to help preprocess dates which are `new_start_date`
and `shift_by`. If `new_start_date` is provided, then dates in the time dimension will
automatically be converted with reference to the `new_start_date` rather than the first date
found in the NetCDF file. The parameter `new_start_date` can be any string parseable by the
[Dates](https://docs.julialang.org/en/v1/stdlib/Dates/) module or a `Dates.DateTime` object.
If additional preprocessing is needed, then one can provide a function that takes in and
returns a `Date.DateTime` object. This function is applied to each date before converting
each dates to seconds with reference with the start date.
```@julia dates_to_seconds
# Shift the dates to first day of month, convert to seconds, and adjust seconds to
# match the date 1/1/2010
obs_var = ClimaAnalysis.OutputVar(
                "pr.nc",
                "precip",
                new_start_date = "2010-01-01T00:00:00", # or Dates.DateTime(2010, 1, 1)
                shift_by = Dates.firstdayofmonth,
            )
```

Additionally, the function `shift_to_start_of_previous_month(var::OutputVar)` is provided to
help with preprocessing. This function shifts the times in the time dimension to the start
of the previous month. After applying this function, the start date in the attributes
corresponds to the first element in the time array.

```@julia beginning
sim_var = shift_to_start_of_previous_month(sim_var)
```

This function is helpful in ensuring consistency in dates between simulation and
observational data. One example of this is when adjusting monthly averaged data. For
instance, suppose that data on 2010-02-01 in `sim_var` corresponds to the monthly
average for January. This function shifts the times so that 2010-01-01 will correspond to
the monthly average for January.


## Integration

`OutputVar`s can be integrated with respect to longitude, latitude, or both using
`integrate_lon(var)`, `integrate_lat(var)`, or `integrate_lonlat(var)` respectively. The
bounds of integration are determined by the range of the dimensions longitude and latitude
in `var`. The unit of both longitude and latitude should be degree.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints for integrating
longitude and latitude.

See the example of integrating over a sphere where the data is all ones to find the surface
area of a sphere.

```@julia integrate_lonlat
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

## Split by season
`OutputVar`s can be split by seasons using `split_by_season(var)` provided that a start date
can be found in `start_date(var)` and time is a dimension in the `OutputVar`. The unit of
time is expected to be second. The function `split_by_season(var)` returns a vector of four
`OutputVar`s with each `OutputVar` corresponding to a season. The months of the seasons are
March to May, June to August, September to November, and December to February. The order of
the vector is MAM, JJA, SON, and DJF. If there are no dates found for a season, then the
`OutputVar` for that season will be an empty `OutputVar`.

```@setup split_by_season
import ClimaAnalysis
import OrderedCollections: OrderedDict

attribs = Dict("start_date" => "2024-1-1");
time = [0.0, 5_184_000.0, 13_132_800.0]; # correspond to dates 2024-1-1, 2024-3-1, 2024-6-1
dims = OrderedDict(["time" => time]);
dim_attribs = OrderedDict(["time" => Dict("units" => "s")]); # unit is second
data = [1.0, 2.0, 3.0];
var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data);
```

```@repl split_by_season
var.attributes
ClimaAnalysis.times(var) # correspond to dates 2024-1-1, 2024-3-1, 2024-6-1
var.data
MAM, JJA, SON, DJF = ClimaAnalysis.split_by_season(var);
isempty(SON) # empty OutputVar because no dates between September to November
[MAM.dims["time"], JJA.dims["time"], DJF.dims["time"]]
[MAM.data, JJA.data, DJF.data]
```

!!! note "`seasons` keyword argument"
    In ClimaAnalysis v0.5.18, the `seasons` keyword argument allows you to change the order
    of the seasons or select only certain seasons. For example, you can call
    `split_by_season(var, seasons = ("DJF", "SON"))` to extract only the seasons `DJF` and
    `SON`, in that order.

### Split by season and year

It may be the case that you want to split a `OutputVar` by season, while keeping each year
separate. This is different from `split_by_season`, which ignores that seasons can come from
different years. This can be done by using `split_by_season_across_time`. For example, if a
`OutputVar` contains times corresponding to 2010-01-01, 2010-03-01, 2010-06-01, 2010-09-01,
and 2010-12-01, then the result of `split_by_season_across_time` is five `OutputVar`s, each
corresponding to a distinct date. Even though 2010-01-01 and 2010-12-01 are in the same
season, there are two `OutputVar`s, because the dates do not belong in the same season and
year.

```@setup split_by_season_across_time
import ClimaAnalysis
import OrderedCollections: OrderedDict

lon = collect(range(-179.5, 179.5, 36))
lat = collect(range(-89.5, 89.5, 18))
time = [0.0]
push!(time, 5_097_600.0) # correspond to 2024-3-1
push!(time, 13_046_400.0) # correspond to 2024-6-1
push!(time, 20_995_200.0) # correspond to 2024-9-1
push!(time, 28_857_600.0) # correspond to 2024-12-1

data = reshape(
    1.0:1.0:(length(lat) * length(time) * length(lon)),
    (length(lat), length(time), length(lon)),
)
dims = OrderedDict(["lat" => lat, "time" => time, "lon" => lon])
attribs = Dict("long_name" => "hi", "start_date" => "2010-1-1")
dim_attribs = OrderedDict([
    "lat" => Dict("units" => "deg"),
    "time" => Dict("units" => "s"),
    "lon" => Dict("units" => "deg"),
])
var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)
```

```@repl split_by_season_across_time
ClimaAnalysis.start_date(var)
ClimaAnalysis.times(var) # dates from the first of January, March, June, August, and December
split_var = ClimaAnalysis.split_by_season_across_time(var);
length(split_var) # months span over 5 seasons
ClimaAnalysis.times(split_var[1]) # correspond to 1/1 (middle of DJF)
ClimaAnalysis.times(split_var[2]) # correspond to 3/1 (start of MAM)
ClimaAnalysis.times(split_var[3]) # correspond to 6/1 (start of JJA)
ClimaAnalysis.times(split_var[4]) # correspond to 9/1 (start of SON)
ClimaAnalysis.times(split_var[5]) # correspond to 12/1 (start of DJF)
```

## Bias and squared error
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
and some observations stored in `ta_1d_average.nc`.

```@julia bias_and_mse
julia> obs_var = OutputVar("ta_1d_average.nc"); # load in observational data

julia> sim_var = get(simdir("simulation_output"), "ta"); # load in simulation data

julia> ClimaAnalysis.short_name(sim_var)
"ta"

julia> bias_var = ClimaAnalysis.bias(sim_var, obs_var); # bias_var is a OutputVar that can be plotted

julia> ClimaAnalysis.global_bias(sim, obs)
2.0

julia> ClimaAnalysis.units(bias_var)
"K"

julia> se_var = ClimaAnalysis.squared_error(sim_var, obs_var); # can also be plotted

julia> ClimaAnalysis.global_mse(sim, obs)
4.0

julia> ClimaAnalysis.global_rmse(sim, obs)
2.0

julia> ClimaAnalysis.units(se_var)
"K^2"
```


### 3D `OutputVar`s
For three-dimensional variables that have dimensions longitude, latitude, and z, pressure,
or time, the functions mentioned before will not work. To compute the bias or squared error,
one can use `ClimaAnalysis.slice` to slice across the z, pressure, or time dimension to get
a 2D variable defined on longitude and latitude. Then, any of the functions mentioned
earlier will work. See an example of this below, where the bias and global MSE are computed
between two `OutputVar`s, where the time dimension is sliced at one day from the start date.

```@julia
# Load in 3D temperature variable defined over longitude, latitude, and time
julia> obs_var = OutputVar("ta_1d_average.nc"); # load in observational data

# Load in 3D temperature variable defined over longitude, latitude, and time
julia> sim_var = get(simdir("simulation_output"), "ta"); # load in simulation data

# Slice to get variables defined over longitude and latitude
julia> obs_var = ClimaAnalysis.slice(obs_var, time =  86400)

julia> sim_var = ClimaAnalysis.slice(sim_var, time =  86400)

julia> ClimaAnalysis.bias(sim_var, obs_var);

julia> ClimaAnalysis.global_mse(sim, obs)
4.0
```

For 3D variables defined over longitude, latitude, and pressure, one can find the global
RMSE in pressure space using `ClimaAnalysis.global_rmse_pfull`. See an example of this
below, where global RMSE is computed between 3D variables in pressure space.

```@julia
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

## Masking
Bias and squared error can be computed only over the land or ocean through the `mask`
parameter. As of now, the mask parameter takes in `apply_oceanmask` or `apply_oceanmask`.
See the example below of this usage.

```julia
# Do not consider the ocean when computing the bias
ClimaAnalysis.bias(sim_var, obs_var, mask = apply_oceanmask)
ClimaAnalysis.global_bias(sim_var, obs_var, mask = apply_oceanmask)

# Do not consider the land when computing the squared error
ClimaAnalysis.squared_error(sim_var, obs_var, mask = apply_landmask)
ClimaAnalysis.global_mse(sim_var, obs_var, mask = apply_landmask)
ClimaAnalysis.global_rmse(sim_var, obs_var, mask = apply_landmask)
```

In other cases, you may want to generate a masking function using a `OutputVar`. For
instance, you are comparing against observational data over some parts of the land, and you
can't use an ocean mask, since not all the observational data is defined over the land. The
function `generate_lonlat_mask` allows you to generate a masking function. If the data is
already zeros and ones, then you can use `generate_lonlat_mask(var, zero_to, one_to)` which
generate a mask where all ones are replaced with `one_to` and zeros are replaced with
`zero_to`. If this is not the case, then you can use `replace` or `replace!` to make the
data all zeros and ones. You can also use the `threshold` keyword argument to adjust how the
values are rounded after resampling from the data of zeros and ones. Values greater than or
equal to `threshold` are rounded to one and all other values are rounded to zero.

```julia
# var is a OutputVar whose data is ones and NaN
var = replace(x -> isnan(x) ? 0.0 : 1.0, var)

# Any points that are NaNs should be zero in the mask
lonlat_mask = ClimaAnalysis.generate_lonlat_mask(
    var,
    NaN, # zero to NaN
    1.0; # one to one
    threshold = 0.99
)

# Apply mask to another OutputVar
another_masked_var = lonlat_mask(another_var)

# Compute squared error and global MSE with custom masking function
ClimaAnalysis.squared_error(sim_var, obs_var, mask = lonlat_mask)
ClimaAnalysis.global_mse(sim_var, obs_var, mask = lonlat_mask)
```
