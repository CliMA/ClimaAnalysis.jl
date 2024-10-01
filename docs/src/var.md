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
to transform between incompatible units). In this case, you an pass a function
that specify how to apply this transformation. For example, in the previous
case, we can assume that we are talking about water and transform units into a
mass flux:
```julia
new_var = ClimaAnalysis.convert_units(var, "kg m/s", conversion_function = (x) -> 1000x)
```

!!! note If you find some unparseable units, please open an issue. We can fix them!

If units do not exist or you want to change the name of the units, then one can uses the
`set_units` function.
```julia
new_var = ClimaAnalysis.set_units(var, "kg m s^-1")
```
!!! warning "Override existing units"
    If units already exist, this will override the units for data in `var`.

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
can be found in `var.attributes["start_date"]` and time is a dimension in the `OutputVar`.
The unit of time is expected to be second. The function `split_by_season(var)` returns a
vector of four `OutputVar`s with each `OutputVar` corresponding to a season. The months of
the seasons are March to May, June to August, September to November, and December to
February. The order of the vector is MAM, JJA, SON, and DJF. If there are no dates found for
a season, then the `OutputVar` for that season will be an empty `OutputVar`.

```@julia split_by_season
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
and some observations stored in "ta\_1d\_average.nc".

```@julia bias_and_mse
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
