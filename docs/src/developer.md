# Developer Documentation

## TemplateVar

!!! warning "Reading OutputVars"
    For most users, it is sufficient to be able to create an `OutputVar` by using `SimDir`
    or by reading it directly from a NetCDF file (e.g. `read_var(nc_data_file)`).

!!! warning "Experimental feature"
    This feature is experimental and the API may change.

The main struct of `ClimaAnalysis` is `OutputVar`. To manually construct an `OutputVar` from
scratch involves either making your own NetCDF file or manually defining the attributes,
dimensions, dimension attributes, and data. This process is prone to errors (e.g. mismatch
between dimension size and data size or mismatch between variable names) and is repetitive
when initializing multiple `OutputVar`s for testing.

To solve this issue, we introduce `TemplateVar`. A `TempleteVar` can be thought of as an
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
    Template.ones_to_n_data(lazy = true) |>
    Template.initialize
```

The conveience functions for adding a dimension are `add_time_dim`, `add_lon_dim`,
`add_lat_dim`, `add_pfull_dim`, and `add_z_dim` which accepts the keyword arguments
`dim_name`, `dim_array`, and `units`. Any other keyword arguments will be attributes for the
dimension. These functions are meant to be used with the pipe operator (`|>`). Additionally,
there are also functions of the form `add_x_dim!` which are meant to be used with function
composition. The `x` can be `time`, `lon`, `lat`, `pfull`, or `z`.

The order of the dimensions added to the `TemplateVar` determines the order of the
dimensions in the `OutputVar`. If the same dimension is added again, then the existing
dimension is overwritten with the new dimension. As a result, when defining a dimension, one
must add all the dimension attributes, name of the dimension, and the dimension array all at
once.

When adding multiple attributes using `add_attribs` or `add_attribs!`, existing attributes
will be overwritten. There is currently no functionaly remove attributes.

When adding multiple functions for specifying the data, only the last function will be used.

There is also the convenience function `make_template_var` which is shown below.

```@example
using ClimaAnalysis.Template # hide
# Order of arguments determine the order of the dimensions in OutputVar
# Order of a TemplateVar is fixed in a TemplateVar
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
Template.add_dim!
Template.ones_data
Template.ones_data!
Template.zeros_data
Template.zeros_data!
Template.ones_to_n_data
Template.ones_to_n_data!
Template.make_template_var
```
