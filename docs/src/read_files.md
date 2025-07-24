# Organizing and loading NetCDF files

`ClimaAnalysis` provides `SimDir` and `NCCatalog` for organizing and loading `OutputVar`s
from NetCDF files. For output from a `CliMA` simulation, use `SimDir`, and for external
NetCDF files, use `NCCatalog`.

## `SimDir`

For information about `SimDir`, see the [`SimDir`](index.md#SimDir) section in the
quickstart guide.

## `NCCatalog`

For NetCDF files not generated from a `CliMA` simulation, you can use
[`ClimaAnalysis.NCCatalog`](@ref) to organize and retrieve variables from the NetCDF files.

The first step is to instantiate a [`ClimaAnalysis.NCCatalog()`](@ref).

```julia
import ClimaAnalysis

catalog = ClimaAnalysis.NCCatalog()
```

To populate the catalog, you can use [`ClimaAnalysis.add_file!`](@ref).

You don't need to provide short names to `add_file!`, but if multiple files contain the same
variable name, `ClimaAnalysis` will load from the file that was added first. To ensure that
the variable is loaded from a specific file, provide the short name (e.g., `"rsdt"`) or an
alias (e.g., `"precip" => "pr"`). If an alias is provided, then the short name of the
`OutputVar` will be renamed accordingly.

```julia
ClimaAnalysis.add_file!(catalog, "gpp.nc")
ClimaAnalysis.add_file!(catalog, "precip.nc", "precip" => "pr")
ClimaAnalysis.add_file!(catalog, "radiation.nc", "rsdt", "rsut", "toa_lw_all_mon" => "rlut")
```

Finally, to retrieve a variable from the `NCCatalog`, you can call `get`. Additional keyword
arguments for constructing the [`ClimaAnalysis.OutputVar`](@ref) can be passed as
`var_kwargs`.

```julia
pr_var = get(catalog, "pr")
rlut_var = get(catalog, "rsut", var_kwargs = (shift_by = Dates.firstdayofmonth,))
```
