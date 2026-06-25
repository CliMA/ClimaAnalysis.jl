export flatten, unflatten, flatten_dim_order, flattened_length

"""
    Representing the metadata of an `OutputVar` and contain all the necessary
    information to reconstruct an `OutputVar`.
"""
struct Metadata{T <: AbstractArray, B, C, FT, N}
    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Order of dimensions when flattening the data of the OutputVar"
    ordered_dims::NTuple{N, String}

    "Bit mask of the dropped values when flattening"
    drop_mask::BitVector

    "Values that are dropped when flattening"
    dropped_values::Vector{FT}
end

"""
    Representing a flat output variable
"""
struct FlatVar{MD <: Metadata, A <: AbstractVector}
    "Metadata associated with the OutputVar"
    metadata::MD

    "Vector that contains all the data"
    data::A
end

"""
    flatten(var::OutputVar;
            dims = ("longitude", "latitude", "pressure_level", "z", "time"),
            ignore_nan = true,
            mask = nothing)

Flatten `var` into a `FlatVar` according to the ordering `dims`.

The default order of dimensions for flattening is
`("longitude", "latitude", "pressure_level", "z", "time")`. Dimensions not present in `var`
are excluded from this ordering.

If `ignore_nan = true`, then `NaNs` are excluded from the flattened data. Otherwise, `NaN`s
are included.

If `mask` is `nothing`, then values are not removed by masking. If `mask` is a
`ClimaAnalysis.LonLatMask`, then values are removed according to the mask. Any value of
`var` whose coordinates correspond to zeros on the mask are excluded from the flattened
data.

!!! note "Mask aware flatten"
    The `mask` keyword argument is only available for versions of ClimaAnalysis after
    v0.5.19.
"""
function flatten(
    var::OutputVar;
    dims = ("longitude", "latitude", "pressure_level", "z", "time"),
    ignore_nan = true,
    mask = nothing,
    drop_mask = nothing,
)
    # The keyword argument `drop_mask` is not exposed to the user and is needed
    # to implement flatten(var::OutputVar, flat_var::FlatVar) and
    # flatten(var::OutputVar, metadata::Metadata). The `drop_mask` is a
    # `BitVector` with length equal to the product of all dimension lengths. If
    # `true`, the corresponding value is dropped. If `false`, the value is kept.
    # When `drop_mask` is provided, both `ignore_nan` and `mask` are ignored.

    # Filter unnecessary dimension names
    dims = conventional_dim_name.(dims)
    var_dims = conventional_dim_name.(collect(dim_names(var)))
    dims = Tuple(
        find_corresponding_dim_name_in_var(dim, var) for
        dim in dims if dim in var_dims
    )

    # Check that all dimensions are present
    length(dims) == length(var.dims) || error(
        "All the dimensions in var ($(dim_names(var))) are not present in dims ($dims)",
    )

    # To flatten data, we need to permute, vectorize, and mask the data
    perm = Tuple(indexin(dims, collect(dim_names(var))))
    permute_data = PermutedDimsArray(var.data, perm)
    vec_data = vec(permute_data)

    if isnothing(drop_mask)
        drop_mask = _drop_mask(var, ignore_nan, mask)
        permute_drop_mask = PermutedDimsArray(drop_mask, perm)
        drop_mask = vec(permute_drop_mask)
    end

    dropped_values = vec_data[drop_mask]
    vec_data = view(vec_data, .!drop_mask)

    # Make metadata
    metadata = Metadata(
        deepcopy(var.attributes),
        deepcopy(var.dims),
        deepcopy(var.dim_attributes),
        dims::NTuple{ndims(var.data), String},
        copy(drop_mask),
        dropped_values,
    )
    return FlatVar(metadata, copy(vec_data))
end

"""
    flatten(var::OutputVar, flat_var::FlatVar)

Flatten `var` into a `FlatVar` according to the `metadata` of `flat_var`.

!!! note "Dropped values"
    This function flattens `var` such that the coordinates of the dropped values
    in `var` are the same as the coordinates of the dropped values in
    `metadata`.
"""
function flatten(var::OutputVar, flat_var::FlatVar)
    return flatten(var, flat_var.metadata)
end

"""
    flatten(var::OutputVar, metadata::Metadata)

Flatten `var` into a `FlatVar` according to `metadata`.

This function will not succeed if:
- the dimension names in the metadata and var are not the same,
- the coordinate values of corresponding dimensions are not approximately equal,
- the units of corresponding dimensions are not the same.

Ordering of the dimensions does not matter.

!!! note "Dropped values"
    This function flattens `var` such that the coordinates of the kept values
    in `var` are the same as the coordinates of the kept values in
    `metadata`.
"""
function flatten(var::OutputVar, metadata::Metadata)
    # Check dimensions are the same in metadata and var
    Set(conventional_dim_name.(dim_names(var))) ==
    Set(conventional_dim_name.(keys(metadata.dims))) || error(
        "Dimensions in var ($(dim_names(var))) is not the same as the dimensions in the metadata ($(keys(metadata.dims)))",
    )

    # Check the dimension arrays are the same
    for var_dim_name in dim_names(var)
        md_dim_name = find_corresponding_dim_name_in_var(var_dim_name, metadata)
        if conventional_dim_name(var_dim_name) != "time"
            all(isapprox(var.dims[var_dim_name], metadata.dims[md_dim_name])) ||
                error(
                    "Dimensions in var ($var_dim_name) and metadata ($md_dim_name) are not the same",
                )
        else
            _dates_or_times(var) == _dates_or_times(metadata) || error(
                "The dates or times between var and metadata are not the same",
            )
        end
    end

    # Check the units of the dimensions are the same
    for var_dim_name in dim_names(var)
        md_dim_name = find_corresponding_dim_name_in_var(var_dim_name, metadata)
        var_dim_units = dim_units(var, var_dim_name)
        md_dim_units = dim_units(metadata, md_dim_name)
        var_dim_units == md_dim_units || error(
            "Units of $var_dim_name in var ($var_dim_units) is not the same as the units of $md_dim_name in metadata",
        )
        if var_dim_units == "" || md_dim_units == ""
            @warn(
                "Units for $(conventional_dim_name(var_dim_name)) is missing in var or metadata"
            )
        end
    end

    return flatten(
        var,
        dims = flatten_dim_order(metadata),
        drop_mask = metadata.drop_mask,
    )
end

"""
    _drop_mask(var::OutputVar, ignore_nan, mask)

Create a bitmask for values to exclude during flattening.

Returns a `BitArray` where `true` indicates a value should be dropped and `false` indicates
a value should be kept in the flattened data.
"""
function _drop_mask(var::OutputVar, ignore_nan, mask)
    bitmask = falses(size(var.data))
    ignore_nan && (bitmask .|= isnan.(var.data))
    if !isnothing(mask)
        bitmask .|= .!(_generate_binary_mask(mask, var))
    end
    return bitmask
end

"""
    unflatten(var::FlatVar)

Unflatten `var` and reconstruct the `OutputVar`.
"""
function unflatten(var::FlatVar)
    return unflatten(var.metadata, var.data)
end

"""
    unflatten(metadata::Metadata, data::AbstractVector)

Unflatten `data` according to the `metadata` and reconstruct the `OutputVar`.

This function assumes that order of `data` before flattened is the same as
`flatten_dim_order(metadata)`.
"""
function unflatten(metadata::Metadata, data::AbstractVector)
    # Check length of data match the length of the dimensions in metadata
    flattened_length_without_dropped =
        flattened_length(metadata, ignore_dropped = true)
    length(data) == flattened_length_without_dropped || error(
        "Flattened data should be of length $flattened_length_without_dropped; got length $(length(data))",
    )

    # To unflatten data, we need to undo the mask, vectorize, and permute. Note that this is
    # the operations of flattening data, but in reverse.
    flattened_length_with_dropped =
        values(metadata.dims) |> collect .|> length |> prod
    flat_data = fill(zero(eltype(data)), flattened_length_with_dropped)
    # Cast to the type of flat_data since NaN can be NaN32 or NaN64
    flat_data[metadata.drop_mask] .= eltype(flat_data).(metadata.dropped_values)
    flat_data[.!metadata.drop_mask] = data

    # Reshape data
    reshape_dims = (length(metadata.dims[dim]) for dim in metadata.ordered_dims)
    unflattened_data = reshape(flat_data, reshape_dims...)

    # Permute dimensions of data
    dim2index = Dict(
        dim_name => index for
        (index, dim_name) in enumerate(keys(metadata.dims))
    )
    perm = invperm(
        collect(dim2index[dim_name] for dim_name in metadata.ordered_dims),
    )
    data = permutedims(unflattened_data, perm)

    return OutputVar(
        deepcopy(metadata.attributes),
        deepcopy(metadata.dims),
        deepcopy(metadata.dim_attributes),
        data,
    )
end

"""
    arecompatible(x::FlatVar, y::FlatVar; ignore_dims = ())

Return whether two `FlatVar` are defined on the same physical space, with values dropped at
the same coordinates.

This means that `x.data[i]` and `y.data[i]` correspond to the same physical coordinates for
every index `i`.

See the documentation for
[`arecompatible(x::ClimaAnalysis.Var.Metadata, y::ClimaAnalysis.Var.Metadata)`](@ref) for
more information.
"""
function arecompatible(x::FlatVar, y::FlatVar; ignore_dims = ())
    return arecompatible(x.metadata, y.metadata; ignore_dims = ignore_dims)
end

"""
    arecompatible(x::Metadata, y::Metadata; ignore_dims = ())

Return whether two `Metadata` are compatible.

Two `Metadata` are compatible when the `FlatVar`s built from them are defined on
the same physical space, with values dropped at the same coordinates.

This is accomplished by comparing `dims` and `dim_attributes` (the latter
because they might contain information about the units), the coordinates at
which the values are dropped, and the order that the dimensions are flattened.
For the time dimension, dates are compared if they are available and times are
compared otherwise.

The keyword argument `ignore_dims` is a dimension name or collection of dimension names
whose coordinate values and dimension units are not compared between the `Metadata`s. The
lengths of these dimensions are still compared, because there is no element-wise
correspondence between `FlatVar`s whose dimensions differ in length.
"""
function arecompatible(x::Metadata, y::Metadata; ignore_dims = ())
    ignore_dims = ignore_dims isa AbstractString ? (ignore_dims,) : ignore_dims
    ignored_names = Set(conventional_dim_name.(ignore_dims))

    # Check the order of the dimensions when flattening are the same and the
    # type of dimensions
    conventional_dim_name.(flatten_dim_order(x)) ==
    conventional_dim_name.(flatten_dim_order(y)) || return false

    for x_dim_name in dim_names(x)
        y_dim_name = find_corresponding_dim_name_in_var(x_dim_name, y)

        # Check the lengths of the corresponding dimensions are the same
        length(x.dims[x_dim_name]) == length(y.dims[y_dim_name]) || return false

        conventional_dim_name(x_dim_name) in ignored_names && continue

        # Check coordinate values of the corresponding dimensions are the same
        if conventional_dim_name(x_dim_name) != "time"
            isapprox(x.dims[x_dim_name], y.dims[y_dim_name]) || return false
        else
            _dates_or_times(x) == _dates_or_times(y) || return false
        end

        # Check the units of the corresponding dimensions are the same
        x_dim_units = dim_units(x, x_dim_name)
        y_dim_units = dim_units(y, y_dim_name)
        x_dim_units == y_dim_units || return false

        # Warn if dimension units are empty
        isempty(x_dim_units) &&
            @warn "Missing units for dimension $x_dim_name in FlatVar/Metadata with short name $(short_name(x))"
        isempty(y_dim_units) &&
            @warn "Missing units for dimension $y_dim_name in FlatVar/Metadata with short name $(short_name(y))"
    end

    # Check the drop_mask are the same, since we care about which coordinates
    # are kept and dropped
    return x.drop_mask == y.drop_mask
end

"""
    _dates_or_times(var_or_metadata)

Helper function to return the coordinate values of the time dimension of `var_or_metadata`,
preferring dates over times.

This function does not check if `var_or_metadata` has a time dimension.
"""
function _dates_or_times(var_or_metadata)
    temporal_dim = try
        dates(var_or_metadata)
    catch
        times(var_or_metadata)
    end
    return temporal_dim
end

"""
    flatten_dim_order(var::FlatVar)

Return the order of the dimensions when flattening the `data` of `var`.
"""
function flatten_dim_order(var::FlatVar)
    return flatten_dim_order(var.metadata)
end

"""
    flatten_dim_order(metadata::Metadata)

Return the order of the dimensions when flattening the `data` of `var` as determined from
the `metadata`.
"""
function flatten_dim_order(metadata::Metadata)
    return metadata.ordered_dims
end

"""
    flattened_length(var::FlatVar; ignore_dropped = true)

Get the length of the flattened `data` according to the `var`.

If `ignore_dropped = true`, then the length is computed excluding the dropped values. If
`ignore_dropped = false`, then the length is computed including the dropped values.
"""
function flattened_length(var::FlatVar; ignore_dropped = true)
    return flattened_length(var.metadata; ignore_dropped = ignore_dropped)
end

"""
    flattened_length(metadata::Metadata; ignore_dropped = true)

Get the length of the flattened `data` according to the `metadata`.

If `ignore_dropped = true`, then the length is computed excluding the dropped values. If
`ignore_dropped = false`, then the length is computed including the dropped values.
"""
function flattened_length(metadata::Metadata; ignore_dropped = true)
    data_size_with_nans = prod(length.(values(metadata.dims)))
    return ignore_dropped ? data_size_with_nans - sum(metadata.drop_mask) :
           data_size_with_nans
end

"""
    Base.getproperty(var::FlatVar, s::Symbol)

Add `attributes`, `dims`, and `dim_attributes` as properties for `var`.

This is necessary to make functions compatible with both `OutputVar` and `FlatVar`.
"""
function Base.getproperty(var::FlatVar, s::Symbol)
    if s === :attributes
        return Base.getfield(Base.getfield(var, :metadata), :attributes)
    elseif s === :dims
        return Base.getfield(Base.getfield(var, :metadata), :dims)
    elseif s === :dim_attributes
        return Base.getfield(Base.getfield(var, :metadata), :dim_attributes)
    else
        return Base.getfield(var, s)
    end
end

"""
    Base.:(==)(metadata1::Metadata, metadata2::Metadata)

Return whether `metadata1` and `metadata2` are equal, using `==` on each field.

This is similar to `isequal`, except `==` treats floating-point `NaN` values as not equal to
each other, `-0.0` as equal to `0.0`, and `missing` as not equal to `missing`.
"""
function Base.:(==)(metadata1::Metadata, metadata2::Metadata)
    # Use & instead of &&, since we want to propagate missing
    return (metadata1.attributes == metadata2.attributes) &
           (metadata1.dims == metadata2.dims) &
           (metadata1.dim_attributes == metadata2.dim_attributes) &
           (metadata1.ordered_dims == metadata2.ordered_dims) &
           (metadata1.drop_mask == metadata2.drop_mask) &
           (metadata1.dropped_values == metadata2.dropped_values)
end

"""
    Base.isequal(metadata1::Metadata, metadata2::Metadata)

Return whether `metadata1` and `metadata2` are equal, using `isequal` on each field.

This is similar to `==`, except `isequal` treats floating-point `NaN` values as equal to
each other, `-0.0` as unequal to `0.0`, and `missing` as equal to `missing`.
"""
function Base.isequal(metadata1::Metadata, metadata2::Metadata)
    return isequal(metadata1.attributes, metadata2.attributes) &&
           isequal(metadata1.dims, metadata2.dims) &&
           isequal(metadata1.dim_attributes, metadata2.dim_attributes) &&
           isequal(metadata1.ordered_dims, metadata2.ordered_dims) &&
           isequal(metadata1.drop_mask, metadata2.drop_mask) &&
           isequal(metadata1.dropped_values, metadata2.dropped_values)
end

"""
    Base.hash(metadata::Metadata, h::UInt)

Compute an integer hash code of `metadata` starting with the hash code `h`.
"""
function Base.hash(metadata::Metadata, h::UInt)
    h = hash(metadata.attributes, h)
    h = hash(metadata.dims, h)
    h = hash(metadata.dim_attributes, h)
    h = hash(metadata.ordered_dims, h)
    h = hash(metadata.drop_mask, h)
    h = hash(metadata.dropped_values, h)
    return h
end

"""
    Base.:(==)(var1::FlatVar, var2::FlatVar)

Return whether `var1` and `var2` are equal, using `==` on each field.

This is similar to `isequal`, except `==` treats floating-point `NaN` values as not equal to
each other, `-0.0` as equal to `0.0`, and `missing` as not equal to `missing`.
"""
function Base.:(==)(var1::FlatVar, var2::FlatVar)
    # Use & instead of &&, since we want to propagate missing
    return (var1.metadata == var2.metadata) & (var1.data == var2.data)
end

"""
    Base.isequal(var1::FlatVar, var2::FlatVar)

Return whether `var1` and `var2` are equal, using `isequal` on each field.

This is similar to `==`, except `isequal` treats floating-point `NaN` values as equal to
each other, `-0.0` as unequal to `0.0`, and `missing` as equal to `missing`.
"""
function Base.isequal(var1::FlatVar, var2::FlatVar)
    return isequal(var1.metadata, var2.metadata) &&
           isequal(var1.data, var2.data)
end

"""
    Base.hash(var::FlatVar, h::UInt)

Compute an integer hash code of `var` starting with the hash code `h`.
"""
function Base.hash(var::FlatVar, h::UInt)
    h = hash(var.metadata, h)
    h = hash(var.data, h)
    return h
end
