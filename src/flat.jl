export flatten, unflatten, flatten_dim_order

"""
    Representing the metadata of an `OutputVar` and contain all the necessary
    information to reconstruct an `OutputVar`.
"""
struct Metadata{T <: AbstractArray, B, C, FT}
    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Order of dimensions when flattening the data of the OutputVar"
    ordered_dims::Tuple{Vararg{String}}

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
)

    # Filter unnecessary dimension names
    dims = conventional_dim_name.(dims)
    var_dims = conventional_dim_name.(collect(keys(var.dims)))
    dims = Tuple(
        find_corresponding_dim_name_in_var(dim, var) for
        dim in dims if dim in var_dims
    )

    # Check that all dimensions are present
    length(dims) == length(var.dims) || error(
        "All the dimensions in var ($(keys(var.dims))) are not present in dims ($dims)",
    )

    # To flatten data, we need to permute, vectorize, and mask the data
    perm = Tuple(indexin(dims, collect(keys(var.dims))))
    permute_data = PermutedDimsArray(var.data, perm)
    vec_data = vec(permute_data)

    drop_mask = _drop_mask(var, ignore_nan, mask)
    permute_drop_mask = PermutedDimsArray(drop_mask, perm)
    drop_mask = vec(permute_drop_mask)

    dropped_values = vec_data[drop_mask]
    vec_data = view(vec_data, .!drop_mask)

    # Make metadata
    metadata = Metadata(
        deepcopy(var.attributes),
        deepcopy(var.dims),
        deepcopy(var.dim_attributes),
        dims,
        copy(drop_mask),
        dropped_values,
    )
    return FlatVar(metadata, copy(vec_data))
end

"""
    _drop_mask(var, mask, ignore_nan)

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

Unflatten `data` according the `metadata` and reconstruct the `OutputVar`.

This function assumes that order of `data` before flattened is the same as
`flatten_dim_order(metadata)`.
"""
function unflatten(metadata::Metadata, data::AbstractVector)
    # Check length of data match the length of the dimensions in metadata
    data_size = _data_length(metadata, ignore_dropped = true)
    length(data) == data_size || error(
        "Flattened data should be of length $data_size; got length $(length(data))",
    )

    # To unflatten data, we need to undo the mask, vectorize, and permute. Note that this is
    # the operations of flattening data, but in reverse.
    data_length = values(metadata.dims) |> collect .|> length |> prod
    flat_data = fill(zero(eltype(data)), data_length)
    # Cast to the type of flat_data since NaN can be NaN32 or NaN64
    flat_data[metadata.drop_mask] .= eltype(flat_data).(metadata.dropped_values)
    flat_data[.!metadata.drop_mask] = data

    # Reshape data
    reshape_dims = (length(metadata.dims[dim]) for dim in metadata.ordered_dims)
    unflattened_data = reshape(flat_data, reshape_dims...)

    # Permute dimensions of data
    dim2index = Dict([
        dim_name => index for
        (index, dim_name) in enumerate(keys(metadata.dims))
    ])
    perm = invperm(
        collect(dim2index[dim_name] for dim_name in metadata.ordered_dims),
    )
    data = permutedims(unflattened_data, perm)

    return OutputVar(
        deepcopy(metadata.attributes),
        deepcopy(metadata.dims),
        deepcopy(metadata.dim_attributes),
        copy(data),
    )
end

"""
    flatten_dim_order(var::FlatVar)

Return the order of the dimensions before flattening `data`.
"""
function flatten_dim_order(var::FlatVar)
    return flatten_dim_order(var.metadata)
end

"""
    flatten_dim_order(metadata::Metadata)

Return the order of the dimensions before flattening `data`.
"""
function flatten_dim_order(metadata::Metadata)
    return metadata.ordered_dims
end

"""
    _data_length(var::FlatVar; ignore_dropped = true)

Get the length of `var.data`.

If `ignore_dropped = true`, then the length is computed excluding the dropped values. If
`ignore_dropped = false`, then the length is computed including the dropped values.
"""
function _data_length(var::FlatVar; ignore_dropped = true)
    return _data_length(var.metadata; ignore_dropped = ignore_dropped)
end

"""
    _data_length(metadata::Metadata; ignore_dropped = true)

Get the length of the flattened `data` according to the `metadata`.

If `ignore_dropped = true`, then the length is computed excluding the dropped values. If
`ignore_dropped = false`, then the length is computed including the dropped values.
"""
function _data_length(metadata::Metadata; ignore_dropped = true)
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
