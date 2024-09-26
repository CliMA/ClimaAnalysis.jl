# How-to guide and cookbook

## How do I make a line plot with variable on the y axis instead of the x one?

By default, the plotting functions in `MakieExt` place the variable on the
x axis. If you want it on the y axis instead (e.g., you are plotting the
vertical profile of a column), you can pass the `dim_on_y = true` argument to
the axis.

For instance,
```julia
plot!(var, more_kwargs = Dict(:axis => [:dim_on_y = true]))
```

`ClimaAnalysis.Utils` provides a convenience function `kwargs` to specify
arguments a little bit more easily without having to think about `Symbol`s too
much.
```julia
plot!(var, more_kwargs = Dict(:axis => kwargs(dim_on_y = true))
```

## How do I take an average of a a variable in a given window of time?

You can use the `window` function to select a portion of a given `var`. For
example, to select only the time from 10 seconds to 100 for `var`

```julia
reduced_var = window(var, "time", left = 10, right = 100)
```

Now, you can apply the usual average functions.

## How do I make the y axis logscale?

The plotting routines can pass additional arguments to `Makie` through the
`more_kwargs` keyword arguments. `more_kwargs` has to be a dictionary that maps
symbols to a list of pairs. The accepted symbols are `:axis`, `:plot`, and
`:cb`, and the pairs have to be pairs of symbols to values (the keyword
arguments you want to pass down). Given that the type structure is a little
complicated, `ClimaAnalysis` comes with a helper function for you to use. So, if
you want to set the logscale for the `y` axis, you would do something like
```julia
import ClimaAnalysis.Utils: kwargs as ca_kwargs

plot!(fig, var, more_kwargs = Dict(:axis => ca_kwargs(yscale = log)))
```
where inside `ca_kwargs` you pass the arguments you would pass to `Makie.Axis`.

## How do I center my longitude to 180 instead of 0?

You can use the `center_longitude!` function.

## How do I find the specific name of a dimension in a `OutputVar` given its conventional name?

Suppose you want to extract the `longitudes` of your data but you don't know if
the dimension was called `lon` or `long`. In this case, you can use the
`Var.longitude_name` function to find the name. This function scans the names
and compares them to a standard list in `ClimaAnalysis.Var.LONGITUDE_NAMES`.
You can also customize that variable if your name is not in that list.
Equivalent functions exist for others dimensions too.

## How do I move to pressure coordinates?

The `Atmos` module in `ClimaAnalysis` comes with a function,
`to_pressure_coordinates` that does precisely that. The function takes an input
`OutputVar` and a pressure `OutputVar`. If the two are compatible, a new
`OutputVar` is returned where the values are linearly interpolated on fixed
pressure levels.

## How do I reorder the dimensions in a `OutputVar` to match the dimensions in another `OutputVar`?

You can use the `reordered_as(src_var, dest_var)` function where `src_var` is a `OutputVar`
with the dimensions you want to reorder to match the dimensions in the OutputVar `dest_var`.

```@setup reordered_as
import ClimaAnalysis
import OrderedCollections: OrderedDict

src_long = 0.0:180.0 |> collect
src_lat = 0.0:90.0 |> collect
src_data = ones(length(src_long), length(src_lat))
src_dims = OrderedDict(["long" => src_long, "lat" => src_lat])
src_attribs = Dict("long_name" => "hi")
src_dim_attribs = OrderedDict([
    "long" => Dict("units" => "test_units1"),
    "lat" => Dict("units" => "test_units2"),
])
src_var =
    ClimaAnalysis.OutputVar(src_attribs, src_dims, src_dim_attribs, src_data)

dest_long = 20.0:180.0 |> collect
dest_lat = 30.0:90.0 |> collect
dest_data = zeros(length(dest_lat), length(dest_long))
dest_dims = OrderedDict(["lat" => dest_lat, "long" => dest_long])
dest_attribs = Dict("long_name" => "hi")
dest_dim_attribs = OrderedDict([
    "lat" => Dict("units" => "test_units4"),
    "long" => Dict("units" => "test_units3"),
])
dest_var = ClimaAnalysis.OutputVar(
    dest_attribs,
    dest_dims,
    dest_dim_attribs,
    dest_data,
)
```

```@repl reordered_as
src_var.dims |> keys |> collect
dest_var.dims |> keys |> collect
reordered_var = ClimaAnalysis.reordered_as(src_var, dest_var);
reordered_var.dims |> keys |> collect
```

## How do I resample the data in a `OutputVar` using the dimensions from another `OutputVar`?

You can use the `resampled_as(src_var, dest_var)` function where `src_var` is a
OutputVar with the data you want to resample using the dimensions in another
OutputVar `dest_var`. If resampling is possible, then a new `OutputVar` is
returned where the data in `src_var` is resampled using a linear interpolation
to fit the dimensions in `dest_var`. Resampling is not possible when the
dimensions in either `OutputVar`s are missing units, the dimensions between the
`OutputVar`s do not agree, or the data in `src_var` is not defined everywhere on
the dimensions in `dest_var`.

```@julia resampled_as
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
