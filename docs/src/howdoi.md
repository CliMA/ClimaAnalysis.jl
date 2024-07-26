# How-to guide and cookbook

## How do I make a line plot with variable on the y axis instead of the x one?

By default, the plotting functions in `CairoMakieExt` place the variable on the
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

