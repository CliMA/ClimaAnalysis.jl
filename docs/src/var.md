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
