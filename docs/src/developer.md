# Developer Documentation

## TemplateVar

!!! warning "Reading NetCDF files"
    For most users, it is sufficient to be able to create a `OutputVar` by using `SimDir`
    or by reading it directly from a NetCDF file (e.g. `read_var(nc_data_file)`).

!!! warning "Experimental feature"
    This feature is experimental, and the API may change.

The main struct in `ClimaAnalysis` is `OutputVar`. To manually construct a `OutputVar` from
scratch involves either making your own NetCDF file or manually defining the attributes,
dimensions, dimension attributes, and data. This process is prone to errors (e.g.,
mismatches between dimension size and data size or mismatch between variable names) and is
repetitive when initializing multiple `OutputVar`s for testing.

To solve this issue, we introduce `TemplateVar`. A `TemplateVar` can be thought of as an
uninitialized `OutputVar` which is not materialized until it is called with `initialize`.
See the example below of constructing a `TemplateVar` and initializing it into an
`OutputVar`.

```@example
using ClimaAnalysis.Template # hide
var =
    Template.TemplateVar() |>
    Template.add_attribs(long_name = "Test", short_name = "test") |>
    Template.add_time_dim(name = "t", dim = collect(0.0:2.0)) |>
    Template.add_lon_dim(units = "degrees") |>
    Template.add_lat_dim(axis = "Y") |>
    Template.one_to_n_data(collected = true) |>
    Template.initialize
```

The convenience functions for adding a dimension are `add_time_dim`, `add_lon_dim`,
`add_lat_dim`, `add_pfull_dim`, and `add_z_dim` which accept the keyword arguments
`dim_name`, `dim_array`, and `units`. Any other keyword arguments will be used as attributes
for the dimension. These functions are designed to be used with the pipe operator (`|>`).
Additionally, there are also functions of the form `add_x_dim!` which are designed to be
used with function composition. The `x` can be `time`, `lon`, `lat`, `pfull`, or `z`.

The order in which dimensions are added to the `TemplateVar` determines their order in the
`OutputVar`. If the same dimension is added again, then the existing dimension is
overwritten with the new dimension. As a result, when defining a dimension, one must add all
the dimension attributes, name of the dimension, and the dimension array at once.

When adding multiple attributes using `add_attribs` or `add_attribs!`, existing attributes
will be overwritten. There is currently no functionality to remove attributes.

If multiple functions are provided for specifying the data, only the last one will be used.

There is also the convenience function `make_template_var` which is shown below.

```@example
using ClimaAnalysis.Template # hide
# Order of arguments determine the order of the dimensions in OutputVar
# Names passed to `make_template_var` is not used for the names of the dimensions of the
# OutputVar
var =
    Template.make_template_var("t", "long", "lat"; abc = "def") |>
    Template.add_lon_dim(name = "lon", dim = [1.0, 2.0]) |>
    Template.add_time_dim(name = "t", dim = [3.0, 4.0]) |>
    Template.initialize
```

# API

```@meta
CurrentModule = ClimaAnalysis
```

## Template

```@docs
Template.TemplateVar()
Template.initialize
Template.add_attribs
Template.add_attribs!
Template.add_dim
Template.add_dim!
Template.ones_data
Template.ones_data!
Template.zeros_data
Template.zeros_data!
Template.one_to_n_data
Template.one_to_n_data!
Template.make_template_var
```
