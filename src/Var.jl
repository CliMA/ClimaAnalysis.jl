module Var

import Dates
import NCDatasets
import OrderedCollections: OrderedDict

import Interpolations as Intp
import Statistics
import Statistics: mean
import NaNStatistics: nanmean, nanvar

import ..Numerics
import ..Utils:
    nearest_index,
    seconds_to_prettystr,
    squeeze,
    split_by_season,
    split_by_season_across_time,
    split_by_month,
    time_to_date,
    date_to_time,
    _data_at_dim_vals,
    _isequispaced,
    find_season_and_year

export OutputVar,
    read_var,
    average_lat,
    weighted_average_lat,
    average_lon,
    average_x,
    average_y,
    average_xy,
    average_time,
    average_lonlat,
    weighted_average_lonlat,
    variance_lon,
    variance_lat,
    variance_time,
    is_z_1D,
    arecompatible,
    shift_longitude,
    center_longitude!,
    short_name,
    long_name,
    units,
    dim_units,
    range_dim,
    reordered_as,
    resampled_as,
    has_units,
    convert_units,
    convert_dim_units,
    integrate_lonlat,
    integrate_lon,
    integrate_lat,
    isempty,
    split_by_season,
    split_by_season_across_time,
    split_by_month,
    average_season_across_time,
    bias,
    global_bias,
    squared_error,
    global_mse,
    global_rmse,
    set_units,
    set_dim_units!,
    shift_to_start_of_previous_month,
    replace,
    replace!,
    reverse_dim,
    reverse_dim!,
    remake,
    permutedims

"""
    Representing an output variable
"""
struct OutputVar{T <: AbstractArray, A <: AbstractArray, B, C <: AbstractDict}

    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Array that contains all the data"
    data::A

    "Dictionary that maps dimension name to its array index"
    dim2index::Dict{String, Int}

    "Array that maps name array index to the dimension name"
    index2dim::Vector{String}
end

"""
    _make_interpolant(dims, data)

Make a linear interpolant from `dims`, a dictionary mapping dimension name to array and
`data`, an array containing data. Used in constructing a `OutputVar`.

If any element of the arrays in `dims` is a Dates.DateTime, then no interpolant is returned.
Interpolations.jl does not support interpolating on dates. If the longitudes span the entire
range and are equispaced, then a periodic boundary condition is added for the longitude
dimension. If the latitudes span the entire range and are equispaced, then a flat boundary
condition is added for the latitude dimension. In all other cases, an error is thrown when
extrapolating outside of `dim_array`.
"""
function _make_interpolant(dims, data)
    # If any element is DateTime, then return nothing for the interpolant because
    # Interpolations.jl do not support DateTimes
    for dim_array in values(dims)
        eltype(dim_array) <: Dates.DateTime && return nothing
    end

    # We can only create interpolants when we have 1D dimensions
    if isempty(dims) || any(d -> ndims(d) != 1 || length(d) == 1, values(dims))
        return nothing
    end

    # Dimensions are all 1D, check that the knots are in increasing order (as required by
    # Interpolations.jl)
    for (dim_name, dim_array) in dims
        if !issorted(dim_array)
            @warn "Dimension $dim_name is not in increasing order. An interpolant will not be created. See Var.reverse_dim and Var.reverse_dim! if the dimension is in decreasing order"
            return nothing
        end
    end

    # Find boundary conditions for extrapolation
    extp_bound_conds = (
        _find_extp_bound_cond(dim_name, dim_array) for
        (dim_name, dim_array) in dims
    )

    dims_tuple, data = _add_extra_lon_point(dims, data)

    extp_bound_conds_tuple = tuple(extp_bound_conds...)
    return Intp.extrapolate(
        Intp.interpolate(dims_tuple, data, Intp.Gridded(Intp.Linear())),
        extp_bound_conds_tuple,
    )
end

"""
    _add_extra_lon_point(dims, data)

Helper function for adding a extra longitude point when making an interpolant.

An extra longitude point is added if the points of the longitude dimension represent centers
instead of edges of the cells and the longitude dimension spans all 360 degrees.

For example, if the longitude dimension is [0.5, 1.5, ..., 359.5] or [0.0, 1.0, ..., 359.0],
then an extra longitude point will be added to the interpolant.
"""
function _add_extra_lon_point(dims, data)
    # Do not want to add an extra point to the dimensions of the OutputVar
    dims_tuple = tuple(deepcopy(values(dims))...)

    for (idx, (dim_name, dim_array)) in enumerate(dims)
        min_of_dim, max_of_dim = extrema(dim_array)
        dim_size = max_of_dim - min_of_dim
        dsize = dim_array[begin + 1] - dim_array[begin]
        if conventional_dim_name(dim_name) == "longitude" &&
           _isequispaced(dim_array) &&
           isapprox(dim_size + dsize, 360.0)
            lon = dims_tuple[idx]

            # Append extra lon point
            # For example, if the longitude dimension is [0.5, 1.5, ..., 359.5], then add
            # 360.5. If the value of 360 is evaluated, it should be the average of the
            # points at 0.5 and 359.5 which is what we compute by adding the extra longitude
            # point and using ongrid.
            push!(lon, lon[end] + dsize)

            # Add corresponding lon slice to the end of data along the index
            # corresponding to the longitude dimension
            first_lon_slice = selectdim(data, idx, 1)
            data = stack(
                (eachslice(data, dims = idx)..., first_lon_slice),
                dims = idx,
            )
        end
    end

    return (dims_tuple, data)
end

"""
    _find_extp_bound_cond(dim_name, dim_array)

Find the appropriate boundary condition for the `dim_name` dimension.
"""
function _find_extp_bound_cond(dim_name, dim_array)
    min_of_dim, max_of_dim = extrema(dim_array)
    dim_size = max_of_dim - min_of_dim
    dsize = dim_array[begin + 1] - dim_array[begin]

    # If the dimension array span the entire range and is equispaced, then add the
    # appropriate boundary condition
    # We do not handle the cases when the array is not equispaced
    (
        conventional_dim_name(dim_name) == "longitude" &&
        _isequispaced(dim_array) &&
        isapprox(dim_size + dsize, 360.0)
    ) && return Intp.Periodic()
    (
        conventional_dim_name(dim_name) == "longitude" &&
        (dim_array[end] - dim_array[begin]) ≈ 360.0
    ) && return Intp.Periodic()
    (
        conventional_dim_name(dim_name) == "latitude" &&
        _isequispaced(dim_array) &&
        isapprox(dim_size + dsize, 180.0)
    ) && return Intp.Flat()
    return Intp.Throw()
end

function OutputVar(attribs, dims, dim_attribs, data)
    index2dim = keys(dims) |> collect
    dim2index =
        Dict([dim_name => index for (index, dim_name) in enumerate(keys(dims))])

    # Check if the size of data matches with the size of dims
    if !(
        isempty(dims) ||
        any(d -> ndims(d) != 1 || length(d) == 1, values(dims))
    )
        size_data = size(data)
        for (dim_index, (dim_name, dim_array)) in enumerate(dims)
            dim_length = length(dim_array)
            data_length = size_data[dim_index]
            if dim_length != data_length
                error(
                    "Dimension $dim_name has inconsistent size with provided data ($dim_length != $data_length)",
                )
            end
        end
    end

    return OutputVar(
        attribs,
        OrderedDict(dims),
        OrderedDict(dim_attribs),
        data,
        dim2index,
        index2dim,
    )
end

function OutputVar(dims, data)
    return OutputVar(Dict{String, Any}(), dims, Dict{String, Dict}(), data)
end

"""
    OutputVar(path::Union{String, Vector{String}},
              short_name = nothing;
              new_start_date = nothing,
              shift_by = identity)

Read the NetCDF file(s) in `path` as a `OutputVar`.

If `short_name` is `nothing`, automatically find the name.

If `path` is a `Vector{String}`, then the files should contain the same variable, and the
times should be sorted. The non-temporal dimensions should be the same across all the files.

Dates in the time dimension are automatically converted to seconds with respect to the first
date in the time dimension array or the `new_start_date`. The parameter `new_start_date` can
be any string parseable by the [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/)
module or a `Dates.DateTime` object. The parameter `shift_by` is a function that takes in
Dates.DateTime elements and return Dates.DateTime elements. The start date is added to the
attributes of the `OutputVar`. The parameter `shift_by` is a function that takes in
`Dates.DateTime` elements and returns `Dates.DateTime` elements. This function is applied to
each element of the time array. Shifting the dates and converting to seconds is done in that
order.
"""
function OutputVar(
    path::Union{String, Vector{String}},
    short_name = nothing;
    new_start_date = nothing,
    shift_by = identity,
)
    var = read_var(path, short_name = short_name)
    # Check if it is possible to convert dates to seconds in the time dimension
    if (
        has_time(var) && (
            eltype(times(var)) <:
            Union{Dates.DateTime, NCDatasets.CFTime.AbstractCFDateTime}
        )
    )
        var = _dates_to_seconds(
            read_var(path; short_name),
            new_start_date = new_start_date,
            shift_by = shift_by,
        )
    end
    return var
end

"""
    read_var(path::String; short_name = nothing)

Read the `short_name` variable in the given NetCDF file.

When `short_name` is `nothing`, automatically identify the name of the variable. If multiple
variables are present, the last one in alphabetical order is chosen.

When `units` is among the attributes, try to parse it and convert it into an
[`Unitful`](https://painterqubits.github.io/Unitful.jl) object. `OutputVar`s with `Unitful`
support automatic unit conversions.

If you want to access `units` as a string, look at [`units`](@ref) function.

Example
=========

```julia
simdir = SimDir("my_output")
read_var(simdir.variable_paths["hu"]["inst"])

read_var("my_netcdf_file.nc", short_name = "ts")
```
"""
function read_var(path::String; short_name = nothing)
    NCDatasets.NCDataset(path) do nc
        # First, if short_name is nothing, we have to identify the name of the variable by
        # finding what is not a dimension
        unordered_dims = NCDatasets.dimnames(nc)
        isnothing(short_name) &&
            (short_name = pop!(setdiff(keys(nc), unordered_dims)))
        dims =
            map(NCDatasets.dimnames(nc[short_name])) do dim_name
                return dim_name => Array(nc[dim_name])
            end |> OrderedDict
        attribs = Dict(k => v for (k, v) in nc[short_name].attrib)

        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[short_name])

        return OutputVar(attribs, dims, dim_attribs, data)
    end
end

"""
    read_var(paths::Vector{String}; short_name = nothing)

Read the `short_name` variable in the vector of NetCDF files.

This function automatically aggregates the NetCDF files along the time dimension. If this
is not possible, then an error is thrown. The attributes of the OutputVar is the same as the
attributes of the first NetCDF file.

When `short_name` is `nothing`, automatically identify the name of the variable. If multiple
variables are present, the last one in alphabetical order is chosen.

When `units` is among the attributes, try to parse it and convert it into an
[`Unitful`](https://painterqubits.github.io/Unitful.jl) object. `OutputVar`s with `Unitful`
support automatic unit conversions.

If you want to access `units` as a string, look at [`units`](@ref) function.

Example
=========

```julia
simdir = SimDir("my_output")
read_var(simdir.variable_paths["hu"]["inst"])

read_var(["my_netcdf_file1.nc", "my_netcdf_file2.nc"], short_name = "ts")
```
"""
function read_var(paths::Vector{String}; short_name = nothing)
    # Assume that `paths` contain all the datasets we want to try to stitch together
    # Assume that `paths` is sorted (will rely on simdir to return it sorted)

    # If there is only one path, then use `read_var(path::String)`
    # We do this because SimDir always store file paths in a vector even if there is only one
    # file path
    length(paths) == 1 && return read_var(first(paths), short_name = short_name)

    # Helper function to find the name of the time dimension
    function find_time_name(dim_names)
        for dim_name in dim_names
            dim_name in TIME_NAMES && return dim_name
        end
        return nothing
    end
    ncs = (NCDatasets.NCDataset(path) for path in paths)
    time_names = map(nc -> find_time_name(NCDatasets.dimnames(nc)), ncs)

    # Check if each NCDataset got a time dimension
    if nothing in time_names
        error(
            "Time dimension does not exist in one of the NetCDF files. If you are loading NetCDF files from a CliMA simulation, provide the specific folder in SimDir",
        )
    end

    # Check if the name for the time dimension is the same across alll NetCDF files
    if length(unique(time_names)) != 1
        error(
            " Names of the time dimension are not the same across all NetCDF files. If you are loading NetCDF files from a CliMA simulation, provide the specific folder in SimDir",
        )
    end

    time_dim_name = first(time_names)

    # Check if the start dates are the same
    unordered_dims = (NCDatasets.dimnames(nc) for nc in ncs)
    short_names = (
        !isnothing(short_name) ? short_name :
        pop!(setdiff(keys(nc), unordered_dim)) for
        (unordered_dim, nc) in zip(unordered_dims, ncs)
    )
    start_dates = (
        get(nc[short_name].attrib, "start_date", "") for
        (nc, short_name) in zip(ncs, short_names)
    )
    if length(unique(start_dates)) != 1
        error(
            "Start dates are not the same across all the NetCDF files. If you are loading NetCDF files from a CliMA simulation, provide the specific folder in SimDir",
        )
    end

    # Check if it make sense to concat time dimension arrays
    time_dim_arrays = (Array(nc[time_dim_name]) for nc in ncs)
    time_dim_concat = vcat(time_dim_arrays...)
    if !issorted(time_dim_concat)
        error(
            "Time dimension is not in increasing order after aggregating datasets. If you are loading NetCDF files from a CliMA simulation, provide the specific folder in SimDir",
        )
    end

    NCDatasets.NCDataset(paths, aggdim = time_dim_name) do nc
        # First, if short_name is nothing, we have to identify the name of the variable by
        # finding what is not a dimension
        unordered_dims = NCDatasets.dimnames(nc)
        isnothing(short_name) &&
            (short_name = pop!(setdiff(keys(nc), unordered_dims)))

        dims =
            map(NCDatasets.dimnames(nc[short_name])) do dim_name
                return dim_name => Array(nc[dim_name])
            end |> OrderedDict
        attribs = Dict(k => v for (k, v) in nc[short_name].attrib)

        dim_attribs = OrderedDict(
            dim_name => Dict(nc[dim_name].attrib) for dim_name in keys(dims)
        )
        data = Array(nc[short_name])

        return OutputVar(attribs, dims, dim_attribs, data)
    end
end

"""
    short_name(var::OutputVar)

Return the `short_name` of the given `var`, if available.

If not available, return an empty string.
"""
function short_name(var::OutputVar)
    get(var.attributes, "short_name", "")
end

"""
    long_name(var::OutputVar)

Return the `long_name` of the given `var`, if available.

If not available, return an empty string.
"""
function long_name(var::OutputVar)
    get(var.attributes, "long_name", "")
end

"""
    units(var::OutputVar)

Return the `units` of the given `var`, if available.

If not available, return an empty string.
"""
function units(var::OutputVar)
    string(get(var.attributes, "units", ""))
end

"""
    has_units(var::OutputVar)

Return whether the given `var` has `units` or not.
"""
function has_units(var::OutputVar)
    return haskey(var.attributes, "units")
end

# Implemented in ClimaAnalysisUnitfulExt
function _maybe_convert_to_unitful end

"""
    Var.convert_units(var, new_units; conversion_function = nothing)

Return a new `OutputVar` with converted physical units of `var` to `new_units`, if possible.

Automatic conversion happens when the units for `var` and `new_units` are both parseable.
When `var` does not have units (see also [`Var.has_units`](@ref)) or has no parseable units,
a conversion function `conversion_function` is required.

`conversion_function` has to be a function that takes one data point and returns the
transformed value.

Being parseable means that `Unitful.uparse` can parse the expression. Please, refer to the
documentation for [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) for more
information.

Examples
=======

Let us set up a trivial 1D `OutputVar` with units of meters per second and automatically
convert it to centimeters per second.

```jldoctest example1
julia> values = 0:100.0 |> collect;

julia> data = copy(values);

julia> attribs = Dict("long_name" => "speed", "units" => "m/s");

julia> dim_attribs = Dict{String, Dict}();

julia> var = ClimaAnalysis.OutputVar(attribs, Dict("distance" => values), dim_attribs, data);

julia> ClimaAnalysis.has_units(var)
true

julia> var_cms = ClimaAnalysis.convert_units(var, "cm/s");

julia> extrema(var_cms.data)
(0.0, 10000.0)
```

Not all the units can be properly parsed, for example, assuming `bob=5lol`

```jldoctest example1
julia> attribs = Dict("long_name" => "speed", "units" => "bob/s");

julia> var_bob = ClimaAnalysis.OutputVar(attribs, Dict("distance" => values), dim_attribs, data);

julia> var_lols = ClimaAnalysis.convert_units(var, "lol/s", conversion_function = (x) -> 5x);

julia> extrema(var_lols.data)
(0.0, 500.0)
```

Failure to specify the `conversion_function` will produce an error.
"""
function convert_units end

"""
    Var.convert_dim_units(var, dim_name, new_units; conversion_function = nothing)

Return a new OutputVar with converted physical units of `dim_name` to `new_units` using
`conversion_function`

This function does not support Unitful, so the parameter `conversion_function` must be
supplied.
"""
function convert_dim_units(
    var,
    dim_name,
    new_units;
    conversion_function = nothing,
)
    dim_name = find_corresponding_dim_name_in_var(dim_name, var)
    (
        !haskey(var.dim_attributes, dim_name) ||
        !haskey(var.dim_attributes[dim_name], "units")
    ) && error(
        "Units do not exist for $dim_name. Use Var.set_dim_units! instead",
    )
    isnothing(conversion_function) && error(
        "Provide a conversion function. Converting units for dimensions using Unitful is not currently supported",
    )

    # Convert to new units
    ret_dims = deepcopy(var.dims)
    ret_dims[dim_name] = conversion_function.(ret_dims[dim_name])

    # Set units
    ret_dim_attribs = deepcopy(var.dim_attributes)
    ret_dim_attribs[dim_name]["units"] = new_units

    ret_attribs = deepcopy(var.attributes)
    ret_data = copy(var.data)
    return OutputVar(ret_attribs, ret_dims, ret_dim_attribs, ret_data)
end

"""
    set_units(var::OutputVar, units::AbstractString)

Set `units` for data in `var`.

!!! warning "Override existing units"
    If units already exist, this will override the units for data in `var`. To convert
    units, see [`Var.convert_units`](@ref)
"""
function set_units(var::OutputVar, units::AbstractString)
    converted_var = convert_units(var, units, conversion_function = identity)
    return converted_var
end

"""
    set_dim_units!(var::OutputVar, dim_name::AbstractString, units::AbstractString)

Set `units` for the `dim_name` dimension in `var`.

!!! warning "Override existing units"
    If units already exist for the dimension, this will override the units for the dimension
    in `var`.
"""
function set_dim_units!(
    var::OutputVar,
    dim_name::AbstractString,
    units::AbstractString,
)
    dim_name = find_corresponding_dim_name_in_var(dim_name, var)
    if haskey(var.dim_attributes, dim_name)
        push!(var.dim_attributes[dim_name], "units" => units)
    else
        var.dim_attributes[dim_name] = Dict("units" => units)
    end
    return nothing
end

"""
    remake(var; attributes, dims, dim_attributes, data)

Remake an `OutputVar`. If a keyword argument is not supplied, then it defaults
to `attributes`, `dims`, `dim_attributes`, or `data` of `var`.
"""
function remake(
    var;
    attributes = var.attributes |> deepcopy,
    dims = var.dims |> deepcopy,
    dim_attributes = var.dim_attributes |> deepcopy,
    data = var.data |> copy,
)
    return OutputVar(attributes, dims, dim_attributes, data)
end

"""
    is_z_1D(var::OutputVar)

Return whether the given `var`iable has an altitude dimension that is 1D.

When topography is present, the altitude dimension in the output variable is typically
multidimensional. The dimensions are (X, Y, Z), where (X, Y) are the horizontal dimensions.
In this case, `dims["z"]` is essentially a map that identifies the physical altitude of the
given point.
"""
function is_z_1D(var::OutputVar)
    has_altitude(var) || error("Variable does not have an altitude dimension")

    return length(size(altitudes(var))) == 1
end

"""
    isempty(var::OutputVar)

Determine whether an OutputVar is empty.
"""
function Base.isempty(var::OutputVar)
    # Do not include :interpolant because var.interpolant is Nothing if data is
    # zero dimensional and empty and isempty(Nothing) throws an error
    return map(
        fieldname -> isempty(getproperty(var, fieldname)),
        filter(x -> x != :interpolant, fieldnames(OutputVar)),
    ) |> all
end

function Base.copy(var::OutputVar)
    fields = fieldnames(OutputVar)
    # We have nested containers and we have to make sure we hand ownership off,
    # so we deepcopy
    return OutputVar([deepcopy(getfield(var, field)) for field in fields]...)
end

"""
    _reduce_over(reduction::F,
                 dims,
                 var::OutputVar,
                 args...;
                 kwargs...)

Apply the given reduction over multiple dimensions in `dims`.

`reduction` has to support the `dims` key. Additional arguments are passed to `reduction`.

The return type is an `OutputVar` with the same attributes, the new data, and the dimension
dropped.

Example
=========

Average over latitudes

```julia
import Statistics: mean
long = 0.:180. |> collect
lat = 0.:90. |> collect
data = reshape(1.:91*181., (181, 91))
dims = Dict(["lat" => lat, "long" => long])
var = OutputVar(dims, data)
_reduce_over(mean, "lat", var)
```
"""
function _reduce_over(
    reduction::F,
    dims,
    var::OutputVar,
    args...;
    kwargs...,
) where {F <: Function}
    dims = dims isa String ? (dims,) : dims

    dim_indices = Tuple(var.dim2index[dim_name] for dim_name in dims)

    # squeeze removes the unnecessary singleton dimension
    # TODO: We might substitute squeeze with dropdims
    data = squeeze(
        reduction(var.data, args...; dims = dim_indices, kwargs...),
        dims = dim_indices,
    )

    # If we reduce over a dimension, we have to remove it
    dims_dict = copy(var.dims)
    dim_attributes = copy(var.dim_attributes)
    for dim in dims
        pop!(dims_dict, dim)
        haskey(var.dim_attributes, dim) && pop!(dim_attributes, dim)
    end
    return OutputVar(
        copy(var.attributes),
        dims_dict,
        dim_attributes,
        copy(data),
    )
end

"""
    average_lat(var::OutputVar; ignore_nan = true, weighted = false)

Return a new OutputVar where the values on the latitudes are averaged arithmetically.

When `weighted` is `true`, weight the average over `cos(lat)`.
"""
function average_lat(var; ignore_nan = true, weighted = false)
    if weighted
        var = copy(var)
        lats = latitudes(var)
        abs(maximum(lats)) >= 0.5π ||
            @warn "Detected latitudes are small. If units are radians, results will be wrong"

        weights_1d = cosd.(lats)
        lat_index = var.dim2index[latitude_name(var)]
        weights = ones(size(var.data))
        # Create a bitmask for the NaN's, we use this to remove weights in the normalization (with nanmean)
        nan_mask = ifelse.(isnan.(var.data), NaN, 1.0)
        for index in CartesianIndices(weights)
            index_tuple =
                ntuple(d -> d == lat_index ? Colon() : index[d], ndims(weights))

            weights[index_tuple...] .= weights_1d
            weights[index_tuple...] ./=
                nanmean(nan_mask[index_tuple...] .* weights_1d)
        end
        var.data .*= weights
    end

    reduced_var =
        _reduce_over(ignore_nan ? nanmean : mean, latitude_name(var), var)
    weighted &&
        haskey(var.attributes, "long_name") &&
        (reduced_var.attributes["long_name"] *= " weighted")
    _update_long_name_generic!(reduced_var, var, latitude_name(var), "averaged")


    return reduced_var
end

"""
    weighted_average_lat(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values on the latitudes are averaged arithmetically
with weights of `cos(lat)`.
"""
weighted_average_lat(var; ignore_nan = true) =
    average_lat(var; ignore_nan, weighted = true)

"""

    average_lon(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values on the longitudes are averaged arithmetically.
"""
function average_lon(var; ignore_nan = true)
    reduced_var =
        _reduce_over(ignore_nan ? nanmean : mean, longitude_name(var), var)
    _update_long_name_generic!(
        reduced_var,
        var,
        longitude_name(var),
        "averaged",
    )
    return reduced_var
end

"""
    average_x(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `x` dimension are averaged arithmetically.
"""
function average_x(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "x", var)
    _update_long_name_generic!(reduced_var, var, "x", "averaged")
    return reduced_var
end

"""
    average_y(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along the `y` dimension are averaged arithmetically.
"""
function average_y(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, "y", var)
    _update_long_name_generic!(reduced_var, var, "y", "averaged")
    return reduced_var
end

"""
    average_xy(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values along both horizontal dimensions `x` and `y`
are averaged arithmetically.
"""
function average_xy(var; ignore_nan = true)
    reduced_var = _average_dims(
        var,
        ("x", "y"),
        ignore_nan = ignore_nan,
        update_long_name = false,
    )

    first_x, last_x = range_dim(var, "x")
    first_y, last_y = range_dim(var, "y")
    units_x = dim_units(var, "x")
    units_y = dim_units(var, "y")
    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " averaged horizontally over x ($first_x to $last_x$units_x) and y ($first_y to $last_y$units_y)"
    end

    return reduced_var
end

"""
    weighted_average_lonlat(var; ignore_nan = true)

Return a new `OutputVar` where the values along the longitude and latitude dimensions
are averaged arithmetically with weights of `cos(lat)` along the latitude dimension.

!!! note "Difference from average_lon and weighted_average_lat"
    The computation `average_lon(weighted_average_lat(var))` computes an average of
    averages. This function computes the global latitude-weighted average over all the
    values across both the longitude and latitude dimensions. In particular, the results
    differ when there are `NaN`s.
"""
function weighted_average_lonlat(var; ignore_nan = true)
    return average_lonlat(var; ignore_nan = ignore_nan, weighted = true)
end

"""
    average_lonlat(var; ignore_nan = true)

Return a new `OutputVar` where the values along the longitude and latitude dimensions
are averaged arithmetically.

!!! note "Difference from average_lon and average_lat"
    The computation `average_lon(average_lat(var))` computes an average of averages. This
    function computes the global average over all the values across both the longitude and
    latitude dimensions. In particular, the results will differ when there are `NaN`s.
"""
function average_lonlat(var; ignore_nan = true, weighted = false)
    lat_name = latitude_name(var)
    lon_name = longitude_name(var)
    !weighted &&
        return _average_dims(var, (lat_name, lon_name), ignore_nan = ignore_nan)

    # Treat the weighted case separately
    abs(maximum(latitudes(var))) >= 0.5π ||
        @warn "Detected latitudes are small. If units are radians, results will be wrong"
    function weighted_lat_avg(data, lats, lat_idx; dims, ignore_nan)
        lat_weights = cosd.(lats)
        # Reshape to broadcast correctly
        size_to_reshape =
            (i == lat_idx ? length(lat_weights) : 1 for i in 1:ndims(data))
        lat_weights = reshape(lat_weights, size_to_reshape...)

        # NaNStatistics provides a way to do weighted mean with NaNs
        ignore_nan && return nanmean(data, lat_weights; dims = dims)

        # Otherwise, we do it the weighted mean by hand if ignore_nan is false
        weighted_avg_data = data .* lat_weights
        # Compute normalization term
        dims_to_drop = Tuple(filter!(i -> i ∉ dims, collect(1:ndims(data))))
        lat_weights = dropdims(lat_weights; dims = dims_to_drop)
        normalization = mapslices(data, dims = dims) do lonlat_slice
            mask = ifelse.(isnan.(lonlat_slice), NaN, 1.0)
            mask .*= lat_weights
            sum(mask)
        end
        weighted_avg_data ./= normalization
        weighted_avg_data = sum(weighted_avg_data, dims = dims)
        return weighted_avg_data
    end

    lat_idx = var.dim2index[latitude_name(var)]
    reduced_var = _reduce_over(
        weighted_lat_avg,
        (longitude_name(var), latitude_name(var)),
        var,
        latitudes(var),
        lat_idx,
        ignore_nan = ignore_nan,
    )

    # Update long name in attributes
    haskey(reduced_var.attributes, "long_name") &&
        (reduced_var.attributes["long_name"] *= " weighted")
    _update_long_name_generic!(
        reduced_var,
        var,
        (lat_name, lon_name),
        "averaged",
    )
    return reduced_var
end

"""
    _average_dims(var,
                  dims;
                  ignore_nan = true,
                  update_long_name = true)

Return a new `OutputVar` where the values along the dimensions in `dims` are averaged
arithmetically.

If `update_long_name` is `true`, then the long name is updated by using
`_update_long_name_generic!`.
"""
function _average_dims(var, dims; ignore_nan = true, update_long_name = true)
    function reduction(data; dims, ignore_nan)
        return ignore_nan ? nanmean(data, dims = dims) : mean(data, dims = dims)
    end

    reduced_var = _reduce_over(reduction, dims, var, ignore_nan = ignore_nan)

    update_long_name &&
        _update_long_name_generic!(reduced_var, var, dims, "averaged")
    return reduced_var
end

"""
    average_time(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the values are averaged arithmetically in time.
"""
function average_time(var; ignore_nan = true)
    reduced_var = _reduce_over(ignore_nan ? nanmean : mean, time_name(var), var)
    _update_long_name_generic!(reduced_var, var, time_name(var), "averaged")
    return reduced_var
end

"""
    variance_time(var; ignore_nan = true)

Return a new OutputVar where the values are the variances along the time dimension.

If `corrected` is `true`, then the variance is computed by dividing the sample mean by `n -
1`, whereas if `corrected` is `false`, then the variance is computed by dividing the sample
mean by `n`, where `n` is the number of elements that the variance is being computed over.
"""
function variance_time(var; ignore_nan = true, corrected = true)
    reduced_var = _reduce_over(
        ignore_nan ? nanvar : Statistics.var,
        time_name(var),
        var,
        corrected = corrected,
    )
    _update_long_name_generic!(reduced_var, var, time_name(var), "variance")
    return reduced_var
end

"""
    variance_lon(var; ignore_nan = true)

Return a new OutputVar where the values are the variances along the longitude dimension.

If `corrected` is `true`, then the variance is computed by dividing the sample mean by `n -
1`, whereas if `corrected` is `false`, then the variance is computed by dividing the sample
mean by `n`, where `n` is the number of elements that the variance is being computed over.
"""
function variance_lon(var; ignore_nan = true, corrected = true)
    reduced_var = _reduce_over(
        ignore_nan ? nanvar : Statistics.var,
        longitude_name(var),
        var,
        corrected = corrected,
    )
    _update_long_name_generic!(
        reduced_var,
        var,
        longitude_name(var),
        "variance",
    )
    return reduced_var

end

"""
    variance_lat(var; ignore_nan = true, corrected = true)

Return a new OutputVar where the values are the variances along the latitude dimension.

If `corrected` is `true`, then the variance is computed by dividing the sample mean by `n -
1`, whereas if `corrected` is `false`, then the variance is computed by dividing the sample
mean by `n`, where `n` is the number of elements that the variance is being computed over.
"""
function variance_lat(var; ignore_nan = true, corrected = true)
    reduced_var = _reduce_over(
        ignore_nan ? nanvar : Statistics.var,
        latitude_name(var),
        var,
        corrected = corrected,
    )
    _update_long_name_generic!(reduced_var, var, latitude_name(var), "variance")
    return reduced_var
end

"""
    dim_units(var::OutputVar, dim_name)

Return the `units` of the given `dim_name` in `var`, if available.

If not available, return an empty string.
"""
function dim_units(var::OutputVar, dim_name)
    dim_name = find_corresponding_dim_name_in_var(dim_name, var)
    # Double get because var.dim_attributes is a dictionary whose values are dictionaries
    string(get(get(var.dim_attributes, dim_name, Dict()), "units", ""))
end

"""
    range_dim(var::OutputVar, dim_name)

Return the range of the dimension `dim_name` in `var`.

Range here is a tuple with the minimum and maximum of `dim_name`.
"""
function range_dim(var::OutputVar, dim_name)
    dim_name = find_corresponding_dim_name_in_var(dim_name, var)
    first_elt = first(var.dims[dim_name])
    last_elt = last(var.dims[dim_name])
    return first_elt, last_elt
end

"""
    _update_long_name_generic!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_name::S,
    operation_name) where {S <: String}

Used by reductions (e.g., average) to update the long name of `reduced_var` by describing
the operation being used to reduce the data and the associated units.
"""
function _update_long_name_generic!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_name::S,
    operation_name,
) where {S <: String}
    dim_of_units = dim_units(var, dim_name)
    first_elt, last_elt = range_dim(var, dim_name)

    if haskey(var.attributes, "long_name")
        reduced_var.attributes["long_name"] *= " $operation_name over $dim_name ($first_elt to $last_elt$dim_of_units)"
    end
    return nothing
end

"""
    _update_long_name_generic!(reduced_var::OutputVar,
                               var::OutputVar,
                               dim_names,
                               operation_name,

Used by reductions (e.g., average) to update the long name of `reduced_var` by describing
the operation being used to reduce the data and the associated units.
"""
function _update_long_name_generic!(
    reduced_var::OutputVar,
    var::OutputVar,
    dim_names,
    operation_name,
)
    N = length(dim_names)
    !haskey(reduced_var.attributes, "long_name") && return nothing
    reduced_var.attributes["long_name"] *= " $operation_name over"
    for (i, dim_name) in enumerate(dim_names)
        dim_of_units = dim_units(var, dim_name)
        first_elt, last_elt = range_dim(var, dim_name)
        # For the last item, append "and"
        i == N && (reduced_var.attributes["long_name"] *= " and")
        reduced_var.attributes["long_name"] *= " $dim_name ($first_elt to $last_elt$dim_of_units)"
        # If there are more than 2 dimensions being reduced and it is not the last
        # dimension, append a comma
        i != N && N > 2 && (reduced_var.attributes["long_name"] *= ",")
    end
    return nothing
end

"""
    center_longitude!(var::OutputVar, lon::Real)

Shift the longitudes in `var` so that `lon` is the center one.

This is useful to center the global projection to the 180 meridian instead of the 0.

!!! warn "Deprecated"
    This function is deprecated and users are encouraged to use
    [`shift_longitude`](@ref) instead.
"""
function center_longitude!(var, lon)
    Base.depwarn(
        "This function is deprecated and may be incorrect. Users are encouraged to use
        `center_longitude(var, lower_lon, upper_lon; shift_by = 0.0)` instead.",
        :center_longitude!,
    )
    lon_name = longitude_name(var)

    old_center_lon_index = nearest_index(var.dims[lon_name], lon)
    half_index = div(length(var.dims[lon_name]), 2)

    # half_index = old_index + shift => shift = half_index - old_index
    shift = half_index - old_center_lon_index

    # We do not use circshift! because it can lead to unpredictable problems when mutating
    shifted_lon = circshift(var.dims[lon_name], shift)
    var.dims[lon_name] = shifted_lon

    lon_dim_index = var.dim2index[lon_name]

    # Prepare the shift tuple for the data array: do not shift, except for the dimension
    # corresponding to the longitude
    shift_tuple = zeros(length(var.dims))
    shift_tuple[lon_dim_index] = shift

    shifted_data = circshift(var.data, shift)
    var.data .= shifted_data
end

"""
    shift_longitude(var, lower_lon, upper_lon; shift_by = 0.0)

Shift the longitudes in `var` to `lower_lon` to `upper_lon` degrees.

We assume that the units of the longitude dimension is degrees. If this is not the case,
then one can use [`convert_dim_units`](@ref) to convert the units to degrees.

This function assumes the prime meridian (0th degree) in `var` is the same before and after
shifting the longitudes.

If two points share the same longitude (e.g. -180 degrees and 180 degrees on longitudes
spanning from -180 degrees to 180 degrees), then we remove the last longitude before
shifting the longitudes. This is necessary to prevent duplicated longitudes after shifting
longitudes.

To shift from -180 to 180 degrees to 0 to 360 degrees, use
`shift_longitude(var, 0.0, 360.0)` and to shift from 0 to 360 degrees to -180 to 180
degrees, use `shift_longitude(var, -180.0, 180.0)`.
"""
function shift_longitude(var, lower_lon, upper_lon)
    width = upper_lon - lower_lon
    width ≈ 2π &&
        @warn "Result may be incorrect if radians are used instead of degrees"
    width ≈ 360.0 || error(
        "The range of longitudes ($lower_lon to $upper_lon degrees) is not 360 degrees",
    )

    # It is possible for the first and last longitude to represent the same latitude on the
    # globe. If this is the case, we remove the last point
    lon_name = longitude_name(var)
    lon = copy(var.dims[lon_name])
    if (lon[end] - lon[begin] ≈ 360.0)
        pop!(lon)
    end

    # Center longitudes
    # If the longitude is already in [lower_lon, upper_lon), return it. Otherwise, we want
    # to map the longitude to [lower_lon, upper_lon]. It is simpler to map it to
    # `[0, upper_lon - lower_lon]` first, so we compute
    # `(lon - lower_lon) % (upper_lon - lower_lon)` to map from `lon` to `lon - lower_lon`
    # to a number in `[0, upper_lon - lower_lon]`. Then, we add `lower_lon` to shift the
    # input to be in range `[lower_lon, upper_lon]`. Note that % is mod1 here.
    wrap_longitude(lon, lower_lon, upper_lon) =
        (lon >= upper_lon) || (lon < lower_lon) ?
        mod1(lon - lower_lon, upper_lon - lower_lon) + lower_lon : lon
    lon .= wrap_longitude.(lon, lower_lon, upper_lon)
    sort_indices = sortperm(lon)
    lon .= lon[sort_indices]

    # Rearrange data
    lon_idx = var.dim2index[lon_name]
    index_tuple =
        (i == lon_idx ? sort_indices : Colon() for i in 1:ndims(var.data))
    data = var.data[index_tuple...]

    # Copy lon to remake the OutputVar
    dims = deepcopy(var.dims)
    dims[lon_name] = lon
    return remake(var, dims = dims, data = data)
end

"""
    date_to_time(var::OutputVar, date::Dates.DateTime)

Convert the given calendar date to a time (in seconds) where t=0 is
`var.attributes["start_date"]`.

This function throws an error if this conversion is not possible with the information in
`var`.
"""
function date_to_time(var::OutputVar, date::Dates.DateTime)
    # NOTE: this a method in Utils

    var_has_start_date = haskey(var.attributes, "start_date")

    if var_has_start_date
        return date_to_time(Dates.DateTime(var.attributes["start_date"]), date)
    else
        error(
            "$date is a Date but `var` does not contain `start_date` in the attributes",
        )
    end
end

"""
    (x::OutputVar)(target_coord)

Interpolate variable `x` onto the given `target_coord` coordinate using
multilinear interpolation.

Extrapolation is now allowed and will throw a `BoundsError` in most cases.

If any element of the arrays of the dimensions is a Dates.DateTime, then interpolation is
not possible. Interpolations.jl do not support making interpolations for dates. If the
longitudes span the entire range and are equispaced, then a periodic boundary condition is
added for the longitude dimension. If the latitudes span the entire range and are
equispaced, then a flat boundary condition is added for the latitude dimension. In all other
cases, an error is thrown when extrapolating outside of the array of the dimension.

Example
=======
```jldoctest
julia> import ClimaAnalysis;

julia> time = 100.0:110.0 |> collect;

julia> z = 0.0:20.0 |> collect;

julia> data = reshape(1.0:(11 * 21), (11, 21));

julia> var2d = ClimaAnalysis.OutputVar(Dict("time" => time, "z" => z), data); var2d.([[105., 10.], [105.5, 10.5]])
2-element Vector{Float64}:
 116.0
 122.0
```
"""
function (x::OutputVar)(target_coord)
    itp = _make_interpolant(x.dims, x.data)
    return itp(target_coord...)
end

"""
    arecompatible(x::OutputVar, y::OutputVar)

Return whether two `OutputVar` are defined on the same physical space

This is accomplished by comparing `dims` and `dim_attributes` (the latter because they might contain information about the units).

We assume that:
- `dim2index` and `index2dim` where correctly created and they reflect `dims`
- `data` is also consistent with `dims`,

We also *do not* check units for `data`.
"""
function arecompatible(x::OutputVar, y::OutputVar)
    x_dims = collect(keys(x.dims))
    y_dims = collect(keys(y.dims))
    x_dims_conventional_names = [conventional_dim_name(name) for name in x_dims]
    y_dims_conventional_names = [conventional_dim_name(name) for name in x_dims]
    x_dim_arrays = collect(values(x.dims))
    y_dim_arrays = collect(values(y.dims))
    x_units = (dim_units(x, dim_name) for dim_name in x_dims)
    y_units = (dim_units(y, dim_name) for dim_name in y_dims)


    for (x_dim, x_unit, y_dim, y_unit) in zip(x_dims, x_units, y_dims, y_units)
        x_unit == "" && @warn "Missing units for dimension $x_dim in x"
        y_unit == "" && @warn "Missing units for dimension $y_dim in y"
        x_unit != y_unit && return false
    end
    x_dims_conventional_names != y_dims_conventional_names && return false
    return x_dim_arrays == y_dim_arrays
end

"""
    _check_dims_consistent(x::OutputVar, y::OutputVar, dim_names = nothing)

Check if the number, name, and unit of dimensions in `x` and `y` are consistent.
Order of the dimensions is not checked.

If the unit for a dimension is missing, then the unit is not consistent for that dimension.

If the iterable or string `dim_names` is supplied, then the dimensions in `dim_names` are
the dimensions that are checked.
"""
function _check_dims_consistent(x::OutputVar, y::OutputVar; dim_names = nothing)
    x_dims = x.dims
    y_dims = y.dims
    if dim_names isa AbstractString
        dim_names = [dim_names]
    end
    if !isnothing(dim_names)
        # Determine if dim_names are in both x and y
        x_dim_names = conventional_dim_name.(keys(x.dims))
        y_dim_names = conventional_dim_name.(keys(y.dims))
        dim_names = conventional_dim_name.(collect(dim_names))
        for dim_name in dim_names
            dim_name in x_dim_names || error(
                "Cannot find $dim_name in the dimension names of x ($x_dim_names)",
            )
            dim_name in y_dim_names || error(
                "Cannot find $dim_name in the dimension names of y ($y_dim_names)",
            )
        end
        # Keep only the dimensions we care about
        keep_dim_name(dim_name) = conventional_dim_name(dim_name) in dim_names
        x_dims =
            filter(dim_name_arr -> keep_dim_name(first(dim_name_arr)), x_dims)
        y_dims =
            filter(dim_name_arr -> keep_dim_name(first(dim_name_arr)), y_dims)
    end

    # Check if the number of dimensions is the same
    x_num_dims = length(x_dims)
    y_num_dims = length(y_dims)
    x_num_dims != y_num_dims && error(
        "Number of dimensions do not match between x ($x_num_dims) and y ($y_num_dims)",
    )

    # Check if the dimensions agree with each other (order does not matter)
    conventional_dim_name_x = Set(conventional_dim_name.(keys(x_dims)))
    conventional_dim_name_y = Set(conventional_dim_name.(keys(y_dims)))
    mismatch_conventional_dim_name =
        conventional_dim_name_x != conventional_dim_name_y
    mismatch_conventional_dim_name && error(
        "Dimensions do not agree between x ($conventional_dim_name_x) and y ($conventional_dim_name_y)",
    )

    # Reorder dimensions to check for units
    # We only care about the name of the dimensions
    x_dim_names_reordered = collect(keys(x_dims))
    y_dim_names_reordered = empty(collect(keys(y_dims)))
    for dim_name in keys(x_dims)
        push!(
            y_dim_names_reordered,
            find_corresponding_dim_name(dim_name, keys(y_dims)),
        )
    end

    x_units = [dim_units(x, dim_name) for dim_name in x_dim_names_reordered]
    y_units = [dim_units(y, dim_name) for dim_name in y_dim_names_reordered]

    # Check for any missing units (missing units are represented with an empty string)
    missing_x = (x_units .== "")
    missing_y = (y_units .== "")
    (any(missing_x) && any(missing_y)) && error(
        "Units for dimensions $(x_dim_names_reordered[missing_x]) are missing in x and units for dimensions $(y_dim_names_reordered[missing_y]) are missing in y",
    )
    any(missing_x) && error(
        "Units for dimensions $(x_dim_names_reordered[missing_x]) are missing in x",
    )
    any(missing_y) && error(
        "Units for dimensions $(y_dim_names_reordered[missing_y]) are missing in y",
    )

    # Check if units match between dimensions
    not_consistent_units = (x_units .!= y_units)
    any(not_consistent_units) && error(
        "Units for dimensions $(x_dim_names_reordered[not_consistent_units]) in x is not consistent with units for dimensions $(y_dim_names_reordered[not_consistent_units]) in y",
    )
    return nothing
end

"""
    permutedims(var::OutputVar, perm)

Permute the dimensions of `var` according to `perm`, an iterable of dimension names
specifying the permutation.

The dimension names in `perm` does not need to be the same as the dimensions in `var`. For
example, dimensions with names `lon` and `long` are identified as the longitude dimension.
"""
function Base.permutedims(var::OutputVar, perm)
    # Get the conventional dim names for var and perm
    conventional_dim_name_var = conventional_dim_name.(keys(var.dims))
    conventional_dim_name_perm = conventional_dim_name.(collect(perm))

    # Check if the dimensions are the same (order does not matter)
    Set(conventional_dim_name_var) == Set(conventional_dim_name_perm) || error(
        "Dimensions are not the same between var ($conventional_dim_name_var) and perm ($conventional_dim_name_perm)",
    )

    # Find permutation indices to reorder dims
    reorder_indices =
        indexin(conventional_dim_name_perm, conventional_dim_name_var)

    # Reorder dims, dim_attribs, and data, but not attribs
    ret_dims = deepcopy(var.dims)
    ret_dims = OrderedDict(collect(ret_dims)[reorder_indices])

    # Cannot assume that every dimension is present in dim_attribs so we loop to reorder the
    # best we can and merge with var.dim_attributes to add any remaining pairs to
    # ret_dim_attribs
    ret_dim_attribs = empty(var.dim_attributes)
    var_dim_attribs = var.dim_attributes |> deepcopy
    var_dim_names = collect(keys(var.dims))
    for idx in reorder_indices
        dim_name = var_dim_names[idx]
        haskey(var_dim_attribs, dim_name) &&
            (ret_dim_attribs[dim_name] = var_dim_attribs[dim_name])
    end
    merge!(ret_dim_attribs, var_dim_attribs)

    ret_data = copy(var.data)
    ret_data = permutedims(ret_data, reorder_indices)
    return remake(
        var,
        dims = ret_dims,
        data = ret_data,
        dim_attributes = ret_dim_attribs,
    )
end

"""
    reordered_as(src_var::OutputVar, dest_var::OutputVar)

Reorder the dimensions in `src_var` to match the ordering of dimensions in `dest_var`.
"""
function reordered_as(src_var::OutputVar, dest_var::OutputVar)
    # Get the conventional dim names for both src_var and dest_var
    conventional_dim_name_src = conventional_dim_name.(keys(src_var.dims))
    conventional_dim_name_dest = conventional_dim_name.(keys(dest_var.dims))

    # Check if the dimensions are the same (order does not matter)
    Set(conventional_dim_name_src) == Set(conventional_dim_name_dest) || error(
        "Dimensions are not the same between src ($conventional_dim_name_src) and dest ($conventional_dim_name_dest)",
    )

    return permutedims(src_var, keys(dest_var.dims))
end

"""
    resampled_as(src_var::OutputVar, dest_var::OutputVar, dim_names = nothing)

Resample `data` in `src_var` to `dims` in `dest_var`.

The resampling performed here is a 1st-order linear resampling.

If the string or iterable `dim_names` is `nothing`, then resampling is done over all
dimensions. Otherwise, resampling is done over the dimensions in `dim_names`.

!!! note "Automatic reordering"
    If resampling is done over all dimensions, then reordering the dimensions of the
    resulting `OutputVar` is automatically done.

The correct dimension to resample on is identified by using `conventional_dim_name`.
For example, if `src_var` has a longitude dimension named `lon`, `dest_var` has a longitude
dimension named `long`, and `dim_names` is `longitude`, then resampling is done on the
longitude dimension because `conventional_dim_name` maps `lon`, `long`, and `longitude` to
`longitude`.

!!! compat "`dim_names` keyword argument"
    The keyword argument `dim_names` is introduced in ClimaAnalysis v0.5.14.

"""
function resampled_as(
    src_var::OutputVar,
    dest_var::OutputVar;
    dim_names = nothing,
)
    if dim_names isa AbstractString
        dim_names = [dim_names]
    end
    # Check for units and if the dimensions exist in src_var and dest_var
    _check_dims_consistent(src_var, dest_var, dim_names = dim_names)

    # If dim_names is nothing, then resample over all dimensions
    if isnothing(dim_names)
        return _resampled_as_all(src_var, dest_var)
    end

    # If the dimensions are the same between both OutputVars and dim_names are the same as
    # well, then resample over all dimensions
    src_var_dim_names = Set(conventional_dim_name.(keys(src_var.dims)))
    dest_var_dim_names = Set(conventional_dim_name.(keys(dest_var.dims)))
    conventional_dim_names = Set(conventional_dim_name.(dim_names))
    if (src_var_dim_names == dest_var_dim_names) &&
       (src_var_dim_names == conventional_dim_names)
        return _resampled_as_all(src_var, dest_var)
    end

    return _resampled_as_partial(src_var, dest_var, dim_names)
end

"""
    _resampled_as_all(src_var::OutputVar, dest_var::OutputVar)

Resample `data` in `src_var` to `dims` in `dest_var` over all dimensions.

Reordering is automatically done.
"""
function _resampled_as_all(src_var::OutputVar, dest_var::OutputVar)
    conventional_names_src = collect(conventional_dim_name.(keys(src_var.dims)))
    conventional_names_dest =
        collect(conventional_dim_name.(keys(dest_var.dims)))
    if conventional_names_src != conventional_names_dest
        src_var = reordered_as(src_var, dest_var)
    end

    itp = _make_interpolant(src_var.dims, src_var.data)
    src_resampled_data =
        [itp(pt...) for pt in Base.product(values(dest_var.dims)...)]

    # Make new dimensions for OutputVar
    src_var_ret_dims = empty(src_var.dims)

    # Loop because names could be different in src_var compared to dest_var
    # (e.g., `long` in one and `lon` in the other)
    for (dim_name, dim_data) in zip(keys(src_var.dims), values(dest_var.dims))
        src_var_ret_dims[dim_name] = copy(dim_data)
    end
    return remake(src_var, dims = src_var_ret_dims, data = src_resampled_data)
end

"""
    _resampled_as_partial(src_var::OutputVar, dest_var::OutputVar, dim_names...)

Resample `data` in `src_var` to `dim_names` in `dest_var`.

!!! note "No reordering"
    Dimensions are not reordered in `src_var` to match the order of the dimensions in
    `dest_var` because the dimensions in `dest_var` and `src_var` respectively are not
    necessarily the same.
"""
function _resampled_as_partial(
    src_var::OutputVar,
    dest_var::OutputVar,
    dim_names,
)
    dim_names = conventional_dim_name.(collect(dim_names))

    # Build grid to resample over
    src_var_ret_dims = empty(src_var.dims)
    for (dim_name, dim) in src_var.dims
        if conventional_dim_name(dim_name) in dim_names
            corresponding_dim_name =
                find_corresponding_dim_name(dim_name, keys(dest_var.dims))
            src_var_ret_dims[dim_name] =
                copy(dest_var.dims[corresponding_dim_name])
        else
            src_var_ret_dims[dim_name] = copy(dim)
        end
    end

    itp = _make_interpolant(src_var.dims, src_var.data)
    src_resampled_data =
        [itp(pt...) for pt in Base.product(values(src_var_ret_dims)...)]

    return remake(src_var, dims = src_var_ret_dims, data = src_resampled_data)
end

"""
    resampled_as(src_var::OutputVar, kwargs...)

Resample `data` in `src_var` to dimensions as defined by the keyword arguments.

The resampling performed here is a 1st-order linear resampling.

For example, to resample on the longitude dimension of `[0.0, 1.0, 2.0]`, one can do
`resampled_as(src_var, lon = [0.0, 1.0, 2.0])`.

If the dimensions in `dims` and `src_var.dims` match (ignoring order), the resulting
`OutputVar` will have its dimensions reordered to match the order in `dims`. Otherwise, the
dimensions of the resulting `OutputVar` will remain unchanged. If the dimensions of `dims`
is a strict subset of the dimensions in `src_var`, then partial resampling is done instead.

Note that there is no checking for units of the dimensions.
"""
function resampled_as(src_var::OutputVar; kwargs...)
    # Construct tuple
    dims = (String(dim_name) => dim_array for (dim_name, dim_array) in kwargs)

    # Find corresponding dimension names in src_var
    src_var_dim_names = collect(keys(src_var.dims))
    corresponding_dest_dim_names = collect(
        find_corresponding_dim_name_in_var(dim_name, src_var) for
        (dim_name, _) in dims
    )

    # Construct dest_var to resample on
    reorder_indices = indexin(corresponding_dest_dim_names, src_var_dim_names)
    dest_dims = OrderedDict(dims)

    # Copy attributes from src_var
    dest_dim_attribs = OrderedDict(
        dim_name => src_var.dim_attributes[src_var_dim_names[i]] for
        ((dim_name, _), i) in zip(dims, reorder_indices)
    )
    size_tuple = Tuple(length.(values(dest_dims)))
    total_size = prod(size_tuple)
    data = reshape(1.0:total_size, size_tuple)
    dest_var = remake(
        src_var;
        data = data,
        dims = dest_dims,
        dim_attributes = dest_dim_attribs,
    )

    # Do not pass to resampled_as as we do not want to check units
    return length(dest_dims) == length(src_var.dims) ?
           _resampled_as_all(src_var, dest_var) :
           _resampled_as_partial(src_var, dest_var, keys(dest_dims))
end

"""
    integrate_lonlat(var::OutputVar)

Integrate `data` in `var` on longitude and latitude with a first-order scheme. `data` has to
be discretized on longitude and latitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints for integrating
longitude and latitude. The units for longitude and latitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lonlat(var::OutputVar)
    var_integrate_lon = var |> integrate_lon
    # Update long name so that we get "...integrated over lat... and integrated over lon..."
    # instead of "...integrated over lat... integrated over lon..."
    if haskey(var_integrate_lon.attributes, "long_name")
        var_integrate_lon.attributes["long_name"] *= " and"
    end
    return var_integrate_lon |> integrate_lat
end

"""
    integrate_lon(var::OutputVar)

Integrate `data` in `var` on longitude with a first-order scheme. `data` has to be
discretized on longitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
longitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lon(var::OutputVar)
    has_longitude(var) || error("var does not has longitude as a dimension")
    lon_name = longitude_name(var)
    return _integrate_over_angle(var, Numerics._integrate_lon, lon_name)
end

"""
    integrate_lat(var::OutputVar)

Integrate `data` in `var` on latitude with a first-order scheme. `data` has to be
discretized on latitude.

If the points are equispaced, it is assumed that each point correspond to the midpoint of a
cell which results in rectangular integration using the midpoint rule. Otherwise, the
integration being done is rectangular integration using the left endpoints. The unit for
latitude should be degrees.

All `NaN`s in `var.data` are treated as zeros when integrating.
"""
function integrate_lat(var::OutputVar)
    has_latitude(var) || error("var does not has latitude as a dimension")
    lat_name = latitude_name(var)
    return _integrate_over_angle(var, Numerics._integrate_lat, lat_name)
end

"""
    _integrate_over_angle(var::OutputVar, integrate_on, angle_dim_name)

Integrate `data` in `var` on `angle_dim_name` with a first-order scheme. `data` has to be
discretized on `angle_dim_name`.

`angle_dim_name` is the name of the angle that is being integrated over. `integrate_on` is a
function that compute the integration for data over `angle_dim_name`.
"""
function _integrate_over_angle(var::OutputVar, integrate_on, angle_dim_name)
    # Enforce constraint that unit is degree because we compute integration weights assuming
    # degrees (see functions _integration_weights_lon_left and
    # _integration_weights_lon_equispaced for examples in Numerics.jl)
    deg_unit_names = [
        "degrees",
        "degree",
        "deg",
        "degs",
        "°",
        "degrees_north",
        "degrees_east",
    ]
    angle_dim_unit = dim_units(var, angle_dim_name)
    lowercase(angle_dim_unit) in deg_unit_names ||
        error("The unit for $angle_dim_name is missing or is not degree")

    integrated_var = _reduce_over(
        integrate_on,
        angle_dim_name,
        var,
        var.dims[angle_dim_name],
    )

    _update_long_name_generic!(
        integrated_var,
        var,
        angle_dim_name,
        "integrated",
    )
    return integrated_var
end


"""
    split_by_season(var::OutputVar; seasons = ("MAM", "JJA", "SON", "DJF"))

Return a vector of four `OutputVar`s split by season.

The months of the seasons are March to May (MAM), June to August (JJA), September to
November (SON), and December to February (JDF). The order of the vector is MAM, JJA, SON,
and DJF. If there are no dates found for a season, then the `OutputVar` for that season will
be an empty `OutputVar`. For non-empty OutputVars, the season can be found by
`var.attributes["season"]`.

The function will use the start date in `var.attributes["start_date"]`. The unit of time is
expected to be second.

!!! note "Interpolating between seasons"
    Interpolations will be inaccurate in time intervals outside of their respective season
    for the returned `OutputVar`s. For example, if an `OutputVar` has the dates 2010-2-1,
    2010-3-1, 2010-4-1, and 2011-2-1 after splitting by seasons, then any interpolation in
    time between the dates 2010-4-1 and 2011-2-1 will be inaccurate.

This function differs from `split_by_season_across_time` as `split_by_season_across_time`
splits dates by season for each year.
"""
function split_by_season(var::OutputVar; seasons = ("MAM", "JJA", "SON", "DJF"))
    _check_time_dim(var)
    start_date = Dates.DateTime(var.attributes["start_date"])

    season_dates = split_by_season(
        time_to_date.(start_date, times(var)),
        seasons = seasons,
    )
    season_times =
        (date_to_time.(start_date, season) for season in season_dates)

    season_vars = _split_along_dim(var, time_name(var), season_times)
    for (season, season_var) in zip(seasons, season_vars)
        isempty(season_var) || (season_var.attributes["season"] = season)
    end
    return season_vars
end

"""
    split_by_season_across_time(var::OutputVar)

Split `var` into `OutputVar`s representing seasons, sorted in chronological order. Each
`OutputVar` corresponds to a single season, and the ordering of the `OutputVar`s is
determined by the dates of the season. The return type is a vector of `OutputVar`s.

The months of the seasons are March to May (MAM), June to August (JJA), September to
November (SON), and December to February (DJF). If there are no dates found for a season,
then the `OutputVar` for that season will be an empty `OutputVar`. The first `OutputVar` is
guaranteed to not be empty. For non-empty `OutputVar`s, the season can be found by
`var.attributes["season"]`.

Also, for non-empty `OutputVar`s, the year can be found by `var.attributes["year"]`. The
convention used is that the second month of the season determines the year. For example, the
year of DJF is the same year as Janauary.

The function will use the start date in `var.attributes["start_date"]`. The unit of time is
expected to be second.

This function differs from `split_by_season` as `split_by_season` splits dates by
season and ignores that seasons can come from different years.
"""
function split_by_season_across_time(var::OutputVar)
    _check_time_dim(var)
    start_date = Dates.DateTime(var.attributes["start_date"])

    seasons_across_year_dates =
        split_by_season_across_time(time_to_date.(start_date, times(var)))
    seasons_across_year_times = (
        date_to_time.(start_date, season) for
        season in seasons_across_year_dates
    )

    split_by_season_vars =
        _split_along_dim(var, time_name(var), seasons_across_year_times)

    for var in split_by_season_vars
        if !isempty(var)
            season, year = find_season_and_year(
                first(time_to_date.(start_date, times(var))),
            )
            # This override the attributes if the user sets something for season and year
            # already
            var.attributes["season"] = season
            var.attributes["year"] = string(year)
        end
    end
    return split_by_season_vars
end

"""
    split_by_month(var::OutputVar)

Split `var` into `OutputVar`s representing months, sorted in chronological order. Each
`OutputVar` corresponds to a single month, and the ordering of the `OutputVar`s is
determined by the dates of the month. The return type is a vector of `OutputVar`s.

If there are no dates found for a month, then the `OutputVar` for that season will be an
empty `OutputVar`. For non-empty `OutputVar`s, the month can be found by
`var.attributes["month"]`.

The function will use the start date in `var.attributes["start_date"]`. The unit of time is
expected to be second.
"""
function split_by_month(var::OutputVar)
    _check_time_dim(var)
    start_date = Dates.DateTime(var.attributes["start_date"])

    monthly_dates = split_by_month(time_to_date.(start_date, times(var)))
    monthly_times =
        (date_to_time.(start_date, month) for month in monthly_dates)

    monthly_vars = _split_along_dim(var, time_name(var), monthly_times)

    for (month, var) in enumerate(monthly_vars)
        if !isempty(var)
            var.attributes["month"] = Dates.monthname(month)
        end
    end
    return monthly_vars
end

"""
    check_time_dim(var::OutputVar)

Check time dimension exists, unit for the time dimension is second, and a
start date is present.
"""
function _check_time_dim(var::OutputVar)
    has_time(var) || error("Time is not a dimension in var")
    dim_units(var, time_name(var)) == "s" ||
        error("Unit for time is not second")
    haskey(var.attributes, "start_date") ||
        error("Start date is not found in var")
    return nothing
end

"""
    average_season_across_time(var::OutputVar; ignore_nan = true)

Return a new OutputVar where the seasons are averaged arithmetically in time.

The values of time dimension of `var` correspond to the first date for the season.

The months of the seasons are March to May (MAM), June to August (JJA), September to
November (SON), and December to February (DJF). If there are no dates found for a season,
then the season is skipped. The season can be found by `var.attributes["season"]`, which
returns a vector of season names.

The year can be found by `var.attributes["year"]`, which returns a vector of years as
strings. The convention used is that the second month of the season determines the year. For
example, the year of DJF is the same year as Janauary.
"""
function average_season_across_time(var; ignore_nan = true)
    season_vars = split_by_season_across_time(var)
    nonempty_season_vars = filter(!isempty, season_vars)
    season_names =
        [season_var.attributes["season"] for season_var in nonempty_season_vars]
    year_names =
        [season_var.attributes["year"] for season_var in nonempty_season_vars]
    season_times =
        [first(times(season_var)) for season_var in nonempty_season_vars]

    # Cannot use cat because it is not possible to concatenate along a dimension that does
    # not exist
    # Cannot use average_time because average_time squeeze the singleton dimension
    avg = ignore_nan ? nanmean : mean
    time_idx = var.dim2index[time_name(var)]
    ret_data = cat(
        (avg(var.data, dims = time_idx) for var in nonempty_season_vars)...,
        dims = time_idx,
    )

    ret_dims = deepcopy(var.dims)
    ret_dims[time_name(var)] = season_times

    # Might need to promote the the type of the keys and values
    attribs_tuple =
        tuple(var.attributes..., "season" => season_names, "year" => year_names)
    ret_attribs = Dict(attribs_tuple...)
    reduced_var =
        remake(var, data = ret_data, dims = ret_dims, attributes = ret_attribs)
    _update_long_name_generic!(
        reduced_var,
        var,
        time_name(var),
        "season averaged",
    )
    return reduced_var
end

"""
    _split_along_dim(var::OutputVar, dim_name, split_vectors)

Given `dim_name` in `var`, split the `OutputVar` by the values in `split_vectors`
and return a vector of `OutputVar`s.

For example, if `dim_name = "time" and `split_vectors = [[0.0, 3.0], [2.0,
4.0]]`, the result is a vector of two `OutputVar`s, where the first OutputVar
has a time dimension of `[0.0, 3.0]` and the second OutputVar has a time
dimension of `[2.0, 4.0]`.

If the vector in `split_vectors` is empty, then an empty OutputVar is returned.
Additonally, there is no checks that are performed in the values in the vectors
in `split_vectors` as the nearest values in `var.dims[dim_name]` are used for
splitting.
"""
function _split_along_dim(var::OutputVar, dim_name, split_vectors)
    # Split data by vectors in split_vectors
    split_data = (
        collect(
            _data_at_dim_vals(
                var.data,
                var.dims[dim_name],
                var.dim2index[dim_name],
                split,
            ),
        ) for split in split_vectors
    )

    # Construct an OutputVar for each season
    return map(split_vectors, split_data) do split, data
        if isempty(split)
            dims = empty(var.dims)
            data = similar(var.data, 0)
            return OutputVar(dims, data)
        end
        ret_dims = deepcopy(var.dims)
        ret_dims[dim_name] = split
        remake(var, dims = ret_dims, data = data)
    end
end

"""
    _check_sim_obs_units_consistent(sim::OutputVar, obs::OutputVar, num_dim)

Check if the number of dimensions are two, the `data` in `sim` and `obs` is missing units or
not, and if the units of data are the same in `sim` and `obs`.

This function does not check if the dimensions are longitude and latitude in `sim` and `obs`
because `integrate_lonlat` (in `bias` and `squared_error`) handles that. The function also
does not check if the units of dimensions in `sim` and `obs` are the same because
`resampled_as` (in `bias` and `squared_error`) handles that.
"""
function _check_sim_obs_units_consistent(
    sim::OutputVar,
    obs::OutputVar,
    num_dim,
)
    # Check number of dimensions
    sim_num_dims = length(sim.dims)
    obs_num_dims = length(obs.dims)
    ((sim_num_dims != num_dim) || (obs_num_dims != num_dim)) && error(
        "There are not only $num_dim dimensions in sim ($sim_num_dims) or obs ($obs_num_dims).",
    )

    # Check units for data is not missing
    sim_data_units = units(sim)
    obs_data_units = units(obs)
    sim_data_units == "" && error("Unit is missing in data for sim")
    obs_data_units == "" && error("Unit is missing in data for obs")

    # Check if units of data match between sim and obs
    sim_data_units == obs_data_units || error(
        "Units do not match between the data in sim ($sim_data_units) and obs ($obs_data_units)",
    )
    return nothing
end

"""
    bias(sim::OutputVar, obs::OutputVar; mask = nothing)

Return a `OutputVar` whose data is the bias (`sim.data - obs.data`) and compute the global
bias of `data` in `sim` and `obs` over longitude and latitude. The result is stored in
`var.attributes["global_bias"]`.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`. Attributes in `sim` and `obs` will be thrown away. The long name and short
name of the returned `OutputVar` will be updated to reflect that a bias is computed.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`global_bias`](@ref), [`squared_error`](@ref), [`global_mse`](@ref),
[`global_rmse`](@ref).
"""
function bias(sim::OutputVar, obs::OutputVar; mask = nothing)
    _check_sim_obs_units_consistent(sim, obs, 2)

    # Resample obs on sim to ensure the size of data in sim and obs are the same and the
    # dims are the same
    obs_resampled = resampled_as(obs, sim)

    # Compute bias
    bias = sim - obs_resampled

    # Do this because we do not want to store global bias as a string and unit could be Unitful
    ret_attributes = Dict{keytype(bias.attributes), Any}(bias.attributes)

    # Add units back for bias
    ret_attributes["units"] = units(sim)

    # Add short and long name
    ret_attributes["short_name"] = "sim-obs"
    ret_attributes["long_name"] = "SIM - OBS"
    if !isempty(short_name(sim))
        ret_attributes["short_name"] *= "_" * short_name(sim)
        ret_attributes["long_name"] *= " " * short_name(sim)
    end

    # Compute global bias and store it as an attribute
    !isnothing(mask) && (bias = mask(bias))
    integrated_bias = integrate_lonlat(bias).data
    ones_var = OutputVar(
        bias.attributes,
        bias.dims,
        bias.dim_attributes,
        ones(size(bias.data)),
    )
    !isnothing(mask) && (ones_var = mask(ones_var))
    normalization = integrate_lonlat(ones_var).data
    # Do ./ instead of / because we are dividing between zero dimensional arrays
    global_bias = integrated_bias ./ normalization
    ret_attributes["global_bias"] = global_bias
    return OutputVar(ret_attributes, bias.dims, bias.dim_attributes, bias.data)
end

"""
    global_bias(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global bias of `data` in `sim` and `obs` over longitude and latitude.

This function is currently only implemented for `OutputVar`s with only the dimensions
longitude and latitude. Units must be supplied for data and dimensions in `sim` and `obs`.
The units for longitude and latitude should be degrees. Resampling is done automatically by
resampling `obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`bias`](@ref), [`squared_error`](@ref), [`global_mse`](@ref),
[`global_rmse`](@ref).
"""
function global_bias(sim::OutputVar, obs::OutputVar; mask = nothing)
    bias_var = bias(sim, obs, mask = mask)
    return bias_var.attributes["global_bias"]
end

"""
    squared_error(sim::OutputVar, obs::OutputVar; mask = nothing)

Return a `OutputVar` whose data is the squared error (`(sim.data - obs.data)^2`) and compute
the global mean squared error (MSE) and global root mean squared error (RMSE) of `data` in
`sim` and `obs` over longitude and latitude. The result is stored in `var.attributes["mse"]`
and `var.attributes["rmse"]`.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`. Attributes in `sim` and `obs` will be thrown away. The long name and short
name of the returned `OutputVar` will be updated to reflect that a squared error is computed.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`global_mse`](@ref), [`global_rmse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function squared_error(sim::OutputVar, obs::OutputVar; mask = nothing)
    _check_sim_obs_units_consistent(sim, obs, 2)

    # Resample obs on sim to ensure the size of data in sim and obs are the same and the dims are the same
    obs_resampled = resampled_as(obs, sim)

    # Compute squared error
    # Do not use ^ since ^ is not defined between a OutputVar and Real
    squared_error = (sim - obs_resampled) * (sim - obs_resampled)

    # Do this because we do not want to store global mse and rmse as strings
    ret_attributes = Dict{String, Any}(squared_error.attributes)

    # Add units back for bias
    # Always treat as a string since the string representation of something type Unitful is
    # not always parseable as a Unitful.Unit (see:
    # https://github.com/PainterQubits/Unitful.jl/issues/412)
    ret_attributes["units"] = "($(units(sim)))^2"

    # Add short and long name
    ret_attributes["short_name"] = "(sim-obs)^2"
    ret_attributes["long_name"] = "(SIM - OBS)^2"
    if !isempty(short_name(sim))
        ret_attributes["short_name"] *= "_" * short_name(sim)
        ret_attributes["long_name"] *= " " * short_name(sim)
    end

    # Compute global mse and global rmse and store it as an attribute
    !isnothing(mask) && (squared_error = mask(squared_error))
    integrated_squared_error = integrate_lonlat(squared_error).data
    ones_var = OutputVar(
        squared_error.attributes,
        squared_error.dims,
        squared_error.dim_attributes,
        ones(size(squared_error.data)),
    )
    !isnothing(mask) && (ones_var = mask(ones_var))
    normalization = integrate_lonlat(ones_var).data
    # Do ./ instead of / because we are dividing between zero dimensional arrays
    mse = integrated_squared_error ./ normalization
    ret_attributes["global_mse"] = mse
    ret_attributes["global_rmse"] = sqrt(mse)
    return OutputVar(
        ret_attributes,
        squared_error.dims,
        squared_error.dim_attributes,
        squared_error.data,
    )
end

"""
    global_mse(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global mean squared error (MSE) of `data` in `sim` and `obs` over longitude and
latitude.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`squared_error`](@ref), [`global_rmse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function global_mse(sim::OutputVar, obs::OutputVar; mask = nothing)
    squared_error_var = squared_error(sim, obs, mask = mask)
    return squared_error_var.attributes["global_mse"]
end

"""
    global_rmse(sim::OutputVar, obs::OutputVar; mask = nothing)

Return the global root mean squared error (RMSE) of `data` in `sim` and `obs` over longitude
and latitude.

This function is currently implemented for `OutputVar`s with only the dimensions longitude
and latitude. Units must be supplied for data and dimensions in `sim` and `obs`. The units
for longitude and latitude should be degrees. Resampling is done automatically by resampling
`obs` on `sim`.

The parameter `mask` is a function that masks a `OutputVar`. See [`apply_landmask`](@ref)
and [`apply_oceanmask`](@ref).

See also [`squared_error`](@ref), [`global_mse`](@ref), [`bias`](@ref), [`global_bias`](@ref).
"""
function global_rmse(sim::OutputVar, obs::OutputVar; mask = nothing)
    squared_error_var = squared_error(sim, obs, mask = mask)
    return squared_error_var.attributes["global_rmse"]
end

"""
    _dates_to_seconds(var::OutputVar;
                      new_start_date = nothing,
                      shift_by = identity)

Convert dates in time dimension to seconds with respect to the first date in the time
dimension array or the `new_start_date`.

Dates in the time dimension are automatically converted to seconds with respect to the first
date in the time dimension array or the `new_start_date`. The parameter `new_start_date` can
be any string parseable by the [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/)
module or a `Dates.DateTime` object. The parameter `shift_by` is a function that takes in
Dates.DateTime elements and return Dates.DateTime elements. The start date is added to the
attributes of the `OutputVar`. The parameter `shift_by` is a function that takes in
`Dates.DateTime` elements and returns `Dates.DateTime` elements. This function is applied to
each element of the time array. Shifting the dates and converting to seconds is done in that
order.

Note that this function only works for the time dimension and will not work for the date
dimension.
"""
function _dates_to_seconds(
    var::OutputVar;
    new_start_date = nothing,
    shift_by = identity,
)
    has_time(var) || error(
        "Converting from dates to seconds is only supported for the time dimension",
    )
    eltype(times(var)) <:
    Union{Dates.DateTime, NCDatasets.CFTime.AbstractCFDateTime} ||
        error("Type of time dimension is not dates")

    # Reinterpret everything as DateTime.Dates
    time_arr = copy(times(var))
    if eltype(time_arr) <: NCDatasets.CFTime.AbstractCFDateTime
        time_arr = map(t -> reinterpret(Dates.DateTime, t), time_arr)
    end

    # Preprocess time_arr by shifting dates
    if !isnothing(shift_by)
        time_arr .= shift_by.(time_arr)
    end

    # Convert from dates to seconds using the first date in the time dimension array as the
    # start date or the new_start_date
    start_date = isnothing(new_start_date) ? time_arr[begin] : new_start_date

    # Handle the case if start_date is a DateTime or string; if it is the latter, then try
    # to parse it as a DateTime
    start_date isa AbstractString && (start_date = Dates.DateTime(start_date))
    time_arr = map(date -> date_to_time(start_date, date), time_arr)

    # Remake OutputVar
    ret_attribs = deepcopy(var.attributes)
    ret_attribs["start_date"] = string(start_date) # add start_date as an attribute
    ret_dim_attribs = deepcopy(var.dim_attributes)
    ret_dim_attribs[time_name(var)]["units"] = "s" # add unit
    var_dims = deepcopy(var.dims)
    ret_dims_generator = (
        conventional_dim_name(dim_name) == "time" ? dim_name => time_arr :
        dim_name => dim_data for (dim_name, dim_data) in var_dims
    )
    ret_dims = OrderedDict(ret_dims_generator...)
    return remake(
        var,
        attributes = ret_attribs,
        dims = ret_dims,
        dim_attributes = ret_dim_attribs,
    )
end

"""
    shift_to_start_of_previous_month(var::OutputVar)

Shift the times in the time dimension to the start of the previous month.

After applying this function, the start date in the attributes correspond to the first
element in the time array.

This function is helpful in ensuring consistency in dates between simulation and
observational data. One example of this is when adjusting monthly averaged data. For
instance, suppose that data on 2010-02-01 in the `OutputVar` corresponds to the monthly
average for January. This function shifts the times so that 2010-01-01 will correspond to
the monthly average for January.

Note that this function only works for the time dimension and will not work for the date
dimension.
"""
function shift_to_start_of_previous_month(var::OutputVar)
    # Check if time dimension exists, floats are in the array, and unit of data is
    # second
    has_time(var) || error("Time is not a dimension of var")
    eltype(times(var)) <: Dates.DateTime && ("Dates found in time array")
    dim_units(var, time_name(var)) != "s" &&
        error("Unit of data is not in second")

    # Convert to seconds to dates
    date_arr =
        Dates.Second.(times(var)) .+
        Dates.DateTime.(var.attributes["start_date"])

    # Apply transformations (find first day of month and subtract one month)
    date_arr .=
        date_arr .|> Dates.firstdayofmonth .|> date -> date - Dates.Month(1)

    # Check for duplicate dates
    unique_dates = unique(date_arr)
    unique_dates != date_arr && error(
        "Dates are not unique after applying shift_to_start_of_previous_month",
    )

    # Convert from dates to seconds
    start_date = date_arr[begin]
    time_arr = map(date -> date_to_time(start_date, date), date_arr)

    ret_attribs = deepcopy(var.attributes)
    ret_attribs["start_date"] = string(start_date)
    ret_dims = deepcopy(var.dims)
    ret_dims["time"] = time_arr
    return remake(var, attributes = ret_attribs, dims = ret_dims)
end

"""
    replace(var::OutputVar, old_new::Pair...)

Return a `OutputVar` where, for each pair `old  => new`, all occurences of `old` are
replaced by `new` in `var.data`

This function is useful if there are `NaN`s or `missing` values in the data. For instance,
you want to use the ocean mask, but there are `NaN`s in the ocean. You can replace all the
`NaN` and `missing` values with 0.0 and apply the ocean mask afterward.
"""
function Base.replace(
    var::OutputVar,
    old_new::Pair...;
    count::Integer = typemax(Int),
)
    replaced_data = replace(var.data, old_new..., count = count)
    return remake(var, data = replaced_data)
end

"""
    replace!(var::OutputVar, old_new::Pair...; [count::Integer])

For each pair `old  => new`, all occurences of `old` are replaced by `new` in
`var.data`. See [`replace`](@ref).

This function is useful if there are `NaN`s or `missing` values in the data. For instance,
you want to use the ocean mask, but there are `NaN`s in the ocean. You can replace all the
`NaN` and `missing` values with 0.0 and apply the ocean mask afterward.
"""
function Base.replace!(
    var::OutputVar,
    old_new::Pair...;
    count::Integer = typemax(Int),
)
    replace!(var.data, old_new..., count = count)
    return nothing
end

"""
    replace!(new::Union{Function, Type}, var::OutputVar; [count::Integer])

Return a new `OutputVar where each value of `var.data` is replaced by `new(x)`. If `count`
is specified, then replace at most `count` values in total.
"""
function Base.replace(
    new::Union{Function, Type},
    var::OutputVar;
    count::Integer = typemax(Int),
)
    replaced_data = replace(new, var.data, count = count)
    return remake(var, data = replaced_data)
end

"""
    replace!(new::Union{Function, Type}, var::OutputVar; [count::Integer])

Return each value of `var.data` by `new(x)`. If `count` is specified, then replace at most
`count` values in total.
"""
function Base.replace!(
    new::Union{Function, Type},
    var::OutputVar;
    count::Integer = typemax(Int),
)
    replace!(new, var.data, count = count)
    return nothing
end

"""
    reverse_dim(var::OutputVar, dim_name)

Reverse a dimension by name.

This function is helpful if the order of a dimension need to be reversed, so that an
interpolant can be made.

See also in-place [`reverse_dim!`](@ref).
"""
function reverse_dim(var::OutputVar, dim_name)
    # Check if dim_name exists
    !haskey(var.dims, dim_name) &&
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    # Check if array is 1D
    ndims(var.dims[dim_name]) != 1 &&
        error("Can only reverse 1D array for dimensions")
    var = remake(var)
    reverse_dim!(var, dim_name)
    return var
end

"""
    reverse_dim!(var::OutputVar, dim_name)

Like [`reverse_dim`](@ref), but operates in-place in `var`.

This function is helpful if the order of a dimension need to be reversed, so that an
interpolant can be made.
"""
function reverse_dim!(var::OutputVar, dim_name)
    # Check if dim_name exists
    !haskey(var.dims, dim_name) &&
        error("Var does not have dimension $dim_name, found $(keys(var.dims))")
    # Check if array is 1D
    ndims(var.dims[dim_name]) != 1 &&
        error("Can only reverse 1D array for dimensions")
    reverse!(var.dims[dim_name])
    dim_idx = var.dim2index[dim_name]
    reverse!(var.data, dims = dim_idx)
    return nothing
end

"""
    Base.cat(vars::OutputVar...; dims::String)

Concatenate `vars` along the dimension `dims`.

This function is helpful if you need to concatenate `OutputVar`s after applying
`split_by_season_across_time` or `window`.

This function does not support concatenating `OutputVar`s of different quantities as
determined by the `short_name`.

Attributes that are not `start_date`, `short_name`, or "units" are discarded in the process.

Note that concatenation along multiple dimensions is not possible.
"""
function Base.cat(vars::OutputVar...; dims::String)
    last_var = last(vars)
    vars_dims = (var.dims for var in vars)

    # Find the index of `dim`
    dim_idx =
        last_var.dim2index[find_corresponding_dim_name_in_var(dims, last_var)]

    # Check the number of dimensions for all vars
    all(var -> length(var.dims) == length(last_var.dims), vars) ||
        error("The number of dimensions are not the same across all OutputVars")

    # Check order of dimensions are the same
    vars_dim_names = collect.(keys.(vars_dims))
    all(dim_name -> dim_name == first(vars_dim_names), vars_dim_names) ||
        return error(
            "Not all OutputVars has the same order of dimensions; use reordered_as or permutedims to reorder the dimensions of the OutputVars",
        )

    # Check units of dimensions
    for dim_name in keys(last_var.dim_attributes)
        all(
            var -> dim_units(var, dim_name) == dim_units(last_var, dim_name),
            vars,
        ) || error(
            "Not all OutputVars have the same units for the dimension $dim_name",
        )
    end

    # Check values of dimensions except for dim
    dim_values = values.(vars_dims)
    for (idx, dim_arrays) in enumerate(zip(dim_values...))
        if idx != dim_idx
            all(dim_array -> dim_array ≈ first(dim_arrays), dim_arrays) ||
                return error(
                    "Not all OutputVars have the same values for the $(collect(keys(last_var.dims))[idx]) dimension",
                )
        end
    end

    # Check units of the data
    all(var -> units(var) == units(first(vars)), vars) ||
        error("Not all OutputVars have the same units")

    # Check start date and short name
    keep_attribs = ("start_date", "short_name")
    for attrib in keep_attribs
        all(
            var ->
                get(var.attributes, attrib, nothing) ==
                get(last_var.attributes, attrib, nothing),
            vars,
        ) || error("Not all OutputVars have the same $attrib")
    end

    # Keep attributes of the last OutputVar. There are multiple approaches for choosing
    # which attributes to keep, such as merging the attributes recursively, taking an
    # intersection of the attributes, or keeping the attributes of the last OutputVar. We
    # take a conservative approach that keep only the short name, units, and start date if
    # those exist.
    keep_attribs = (keep_attribs..., "units")
    ret_attribs = empty(last_var.attributes)
    for attrib in keep_attribs
        attrib in keys(last_var.attributes) &&
            (ret_attribs[attrib] = last_var.attributes[attrib])
    end

    # Add long name
    long_name = join(
        (
            var.attributes["long_name"] for
            var in vars if haskey(var.attributes, "long_name")
        ),
        " concatenated ",
    )
    ret_attribs["long_name"] = long_name

    # Concat all the dim of vars
    ret_dims = deepcopy(last_var.dims)
    cat_dim = vcat(
        (
            var.dims[find_corresponding_dim_name_in_var(dims, var)] for
            var in vars
        )...,
    )
    ret_dims[find_corresponding_dim_name_in_var(dims, last_var)] = cat_dim

    # Concat all the data along dim
    ret_data = cat((var.data for var in vars)..., dims = dim_idx)

    # Keep dim attribs of the last OutputVar
    return remake(
        last_var;
        attributes = ret_attribs,
        dims = ret_dims,
        data = ret_data,
    )
end

"""
    Base.show(io::IO, var::OutputVar)

Pretty print the contents of an `OutputVar`.

Print the attributes, dimension attributes, and dimensions that the data is defined
over.
"""
function Base.show(io::IO, var::OutputVar)
    # Print the key value pairs of attributes
    printstyled(io, "Attributes:\n", bold = true, color = :green)
    # Find spacing to pad out key of var.attributes
    max_length_attribs =
        maximum(length(x) for (x, _) in var.attributes; init = 0)
    for (key, val) in var.attributes
        print(io, "  " * rpad(key, max_length_attribs) * " => ")
        printstyled(io, "$(val)\n", color = :light_cyan)
    end

    # Print the key value pairs of dimension attributes recursively
    printstyled(io, "Dimension attributes:\n", bold = true, color = :green)
    # Find spacing to pad out key of var.attributes
    for (dim, dim_attrib) in var.dim_attributes
        printstyled(io, "  $dim:\n", color = :light_green)
        max_length_dim_attribs =
            maximum(length(x) for (x, _) in dim_attrib; init = 0)
        for (key, val) in dim_attrib
            print(io, "    " * rpad(key, max_length_dim_attribs) * " => ")
            printstyled(io, "$(val)\n", color = :light_cyan)
        end
    end

    # Print the dimensions that the data is defined over
    printstyled(io, "Data defined over:", bold = true, color = :green)
    # Do not add a new line if there is nothing to print in var.dims
    !isempty(var.dims) && print(io, "\n")
    max_length_dims = maximum(length(x) for (x, _) in var.dims; init = 0)
    for (i, (dim, array)) in enumerate(var.dims)
        printstyled(io, "  " * rpad(dim, max_length_dims), color = :light_green)
        print(io, " with ")
        printstyled(io, "$(length(array)) ", color = :light_cyan)

        # Print contents depending on size of the dimension
        if length(array) >= 2
            print(io, "elements")
            sorted_at_all = issorted(array) || issorted(array, rev = true)
            if sorted_at_all
                print(io, " (")
                printstyled(io, "$(array[begin])", color = :light_cyan)
                print(io, " to ")
                printstyled(io, "$(array[end])", color = :light_cyan)
                print(io, ")")
            else
                print(io, " (")
                printstyled(io, "not sorted", color = :red)
                print(io, ")")
            end
        elseif length(array) == 1
            print(io, "element")
            print(io, " (")
            printstyled(io, "$(array[begin])", color = :light_cyan)
            print(io, ")")
        elseif length(array) == 0
            print(io, "element")
        end
        # Do not add a new line on the last dimension
        if i != length(var.dims)
            print(io, "\n")
        end
    end
end

"""
    Base.identity(var::OutputVar)

Return the given `var`.
"""
function Base.identity(var::OutputVar)
    return var
end

include("outvar_operators.jl")
include("outvar_dimensions.jl")
include("outvar_selectors.jl")
include("flat.jl")
include("masks.jl")

end
