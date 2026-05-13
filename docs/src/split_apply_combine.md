# Split Apply Combine

!!! note "Experimental"
    This is an initial implementation of the split-apply-combine pattern that
    supports a limited set of features. As a result, the public API will not be
    covered under semantic versioning and is subjected to changes, but we will
    try not to change it.

    If you need a feature that is not yet implemented, please let us know!

The paper ["The Split-Apply-Combine Strategy for Data Analysis"](https://www.jstatsoft.org/article/view/v040i01),
written by Hadley Wickham, illustrates that most data processing follows the
pattern of:
- splitting data into groups,
- applying a reduction, transformation, or filter on each group,
- combining the results together via concatenation.

ClimaAnalysis currently implements splitting operations via a single group via
[`GroupAll`](@ref) or seasonal splits via [`SplitSeason`](@ref) and reductions
via [`Reduce`](@ref).

## Tutorial

The split-apply-combine pattern is only defined for `OutputVar`s. Assuming you
already have `var`, a `OutputVar`, here's a complete example of computing a time
average in a functional style approach.

```@setup split-apply-combine
import Dates
import ClimaAnalysis
import ClimaAnalysis.Template:
    TemplateVar, add_dim, add_attribs, initialize

ref_date = Dates.DateTime(2010)
time =
    ClimaAnalysis.Utils.date_to_time.(
        Ref(ref_date),
        [
            Dates.DateTime(2010, 1),
            Dates.DateTime(2010, 2),
            Dates.DateTime(2010, 3),
            Dates.DateTime(2010, 4),
        ],
    )
lat = [-90.0, 0.0, 90.0]
lon = [-180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0]
var =
    TemplateVar() |>
    add_dim("time", time, units = "s") |>
    add_dim("lat", lat, units = "degrees") |>
    add_dim("lon", lon, units = "degrees") |>
    add_attribs(start_date = ref_date, short_name = "lwu") |>
    initialize
```

```@example split-apply-combine
var
```

```@example split-apply-combine
import ClimaAnalysis
import Statistics: mean

time_averaged_var =
    var |>
    ClimaAnalysis.GroupAll("time") |>
    ClimaAnalysis.Reduce(mean) |>
    ClimaAnalysis.combine
```

!!! note "What is the difference between this example and `average_time`?"
    There are no modifications to the attributes in this example and
    [`average_time`](@ref) does modifies the attributes to reflect that the
    operation happens. Also, `average_time` squeezes the time dimension and
    remove it while the time dimension is kept when using the
    split-apply-combine pattern.

First, a [`AbstractSplitOperation`](@ref) is called on a `OutputVar` which
produces a [`SplitApplyVar`](@ref). In this example, var is piped to
`GroupAll("time")`. The result is a `SplitApplyVar` which represents a lazy
evaluation of the split-apply-combine pattern on a `OutputVar`.

Second, a [`AbstractApplyOperation`](@ref) is called on the `SplitApplyVar`
which record the apply operation. In this example, we compute a time average
over every group which there is only one. For [`Reduce`](@ref), the reduction
function passed must accept an `Array` and a `dims` keyword argument, and must
return an array with the same number of dimensions where the size along `dims`
is 1 (i.e., the dimension is kept but collapsed to a single element). Functions
such as `Statistics.mean`, `sum`, `minimum`, and `maximum` satisfy this
requirement.

Third, [`combine`](@ref) is called on the `SplitApplyVar` which starts the
split, apply, and combine operations to produce the resulting `OutputVar`. The
final result is a time-averaged `OutputVar` that still have a time dimension.
The single value of the time dimension is the first date of the time dimension.
The split dimension is still present in the returned `OutputVar`, but its size
equals the number of groups. For [`GroupAll`](@ref), this is always 1. The
coordinate value for each group is taken from the first element of that group.

```@repl split-apply-combine
ClimaAnalysis.dates(var)
ClimaAnalysis.dates(time_averaged_var)
```
