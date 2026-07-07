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

!!! note "Support for multi-file concatenation"
    In versions of ClimaAnalysis after v0.5.21, you can pass a vector of filepaths to
    `ClimaAnalysis.add_file!`. The NetCDF files are automatically concatenated along the
    time dimension.

    ```julia
    ClimaAnalysis.add_file!(catalog, ["era5_2008.nc", "era5_2009.nc"])
    ```

Finally, to retrieve a variable from the `NCCatalog`, you can call `get`. Additional keyword
arguments for constructing the [`ClimaAnalysis.OutputVar`](@ref) can be passed as
`var_kwargs`.

```julia
pr_var = get(catalog, "pr")
rlut_var = get(catalog, "rsut", var_kwargs = (shift_by = Dates.firstdayofmonth,))
```

If you have organized your NetCDF files into several catalogs, you can combine them into a
single catalog with `merge`. This is useful when different data sources are loaded into
separate catalogs that you later want to query through one catalog.

```julia
era5_catalog = ClimaAnalysis.NCCatalog()
ClimaAnalysis.add_file!(era5_catalog, "era5.nc", "rsdt")

ceres_catalog = ClimaAnalysis.NCCatalog()
ClimaAnalysis.add_file!(ceres_catalog, "ceres.nc", "toa_lw_all_mon" => "rlut")

combined_catalog = merge(era5_catalog, ceres_catalog)
```

The in-place variant `merge!` merges the other catalogs into the first one, modifying it.

!!! note "Duplicate names when merging"
    If the same short name or alias appears in more than one catalog, then the entry from the
    rightmost (last) catalog that owns the name is kept, following the same behavior as
    `merge` for dictionaries.

    ```julia
    # If both catalogs contain "rsdt", the entry from ceres_catalog is kept
    combined_catalog = merge(era5_catalog, ceres_catalog)
    ```
