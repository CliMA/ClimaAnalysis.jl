ClimaAnalysis.jl Release Notes
===============================

v0.5.8
------

## Features

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

## Bug fixes

- Increased the default value for `warp_string` to 72.

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
