# `FlatVar`

To help with data preparation and feeding data from `OutputVar`s into other pipelines,
ClimaAnalysis provides two functions [`flatten`](@ref) to flatten an `OutputVar` into a
`FlatVar` and [`unflatten`](@ref) to unflatten a `FlatVar` into an `OutputVar`. A `FlatVar`
consists of `data` which is a one-dimensional vector of the `OutputVar`'s data and
`metadata` which contain all the necessary information to reconstruct the original
`OutputVar`.

The problems that `FlatVar` aim to solve are

1. Flattening can be error prone if one does not keep track of the ordering of the
   dimensions,
2. Flattening a one-dimensional vector loses important information, such as what quantity is
   represented by the `OutputVar` and the ordering of the dimensions.

## Flatten

To solve the first problem, ClimaAnalysis enforces the ordering of the dimensions when
flattening to be `("longitude", "latitude", "pressure", "altitude", "time")` and omits the
dimensions that do not exist in the `OutputVar`.

```@setup flat
import ClimaAnalysis
import OrderedCollections: OrderedDict

time = [0.0, 1.0, 5.0]
lon = [-60.0, -30.0, 0.0, 30.0, 60.0]
lat = [-90.0, -30.0, 30.0, 90.0]
n_elts = length(lat) * length(lon) * length(time)
dims = OrderedDict("time" => time, "lon" => lon, "lat" => lat)
size_of_data = (length(time), length(lon), length(lat))
data = reshape(collect(1.0:n_elts), size_of_data...)
attribs = Dict("short_name" => "pr")
dim_attribs = OrderedDict(
    "time" => Dict("units" => "s"),
    "lon" => Dict("units" => "degrees_east"),
    "lat" => Dict("units" => "degrees_north"),
)
var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data)

# Construct an OutputVar with 3 NaNs
nan_data = copy(data)
nan_data[1,1,1] = NaN
nan_data[2,1,3] = NaN
nan_data[3,5,4] = NaN
nan_var = ClimaAnalysis.remake(var, data = nan_data)
```

The example below demonstrates how the order of the dimensions of a `OutputVar` does
not matter when flattening the data. The order of dimensions of `flat_var` is
`time`, `lon`, and `lat`, and the order of dimensions of `permuted_var` is `lon`, `lat`, and
`time`. When both `OutputVar`s are flattened, the flattened data of both `OutputVar`s are
the same. As such, using `flatten` removes the need to keep track of the order of the
dimensions of a `OutputVar`.

```@repl flat
var
flat_var = ClimaAnalysis.flatten(var);
permuted_var = permutedims(var, ("longitude", "latitude", "time"));
flat_permuted_var = ClimaAnalysis.flatten(permuted_var);
isequal(flat_permuted_var.data, flat_var.data)
```

The data can be extracted by `flat_var.data` and the metadata can be extracted by
`flatvar.metadata`. More information about the metadata will discussed in the section below.

Furthermore, if `ignore_nan = true`, then `NaNs` are excluded when flattening the data.

```@repl flat
count(isnan, nan_var.data) # nan_var is the same as var, but contains three NaNs
flat_nan_var = ClimaAnalysis.flatten(nan_var, ignore_nan = true); # default is true
length(flat_nan_var.data)
```

## Unflatten

To solve the second problem, there is the `metadata` field in `FlatVar` that stores the
necessary information to fully reconstruct the `OutputVar`. With `unflatten`, one can
call `unflatten` on both `FlatVar` or on `metadata` and `data` to reconstruct the
`OutputVar`. This decoupling of `metadata` and `data` also means that one needs to be
careful that the correct `metadata` and `data` are being used to reconstruct the
`Outputvar`.

See the example below of unflattening `flat_var`.

```@repl flat
unflatten_var = ClimaAnalysis.unflatten(flat_var)
isequal(unflatten_var.data, var.data)
```

One can also unflatten using the metadata and data of `flat_var`.

```@repl flat
unflatten_var = ClimaAnalysis.unflatten(flat_var.metadata, flat_var.data)
isequal(unflatten_var.data, var.data)
```
