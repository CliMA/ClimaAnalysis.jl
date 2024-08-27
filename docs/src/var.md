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
