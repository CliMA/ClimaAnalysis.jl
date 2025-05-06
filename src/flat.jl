export flatten, unflatten, flatten_dim_order

"""
    Representing the metadata of an `OutputVar` and contain all the necessary
    information to reconstruct an `OutputVar`.
"""
struct Metadata{T <: AbstractArray, B, C, I <: Integer}
    "Attributes associated to this variable, such as short/long name"
    attributes::Dict{String, B}

    "Dimensions over which the variable is defined"
    dims::OrderedDict{String, T}

    "Attributes associated to the dimensions"
    dim_attributes::OrderedDict{String, C}

    "Order of dimensions when flattening the data of the OutputVar"
    ordered_dims::Tuple{Vararg{String}}

    "Indices of NaNs in data"
    nan_indices::Vector{I}
end

"""
    Representing a flat output variable
"""
struct FlatVar{A <: AbstractVector}
    "Metadata associated with the OutputVar"
    metadata::Metadata

    "Vector that contains all the data"
    data::A
end

"""
    flatten(var::OutputVar;
            dims = ("longitude", "latitude", "pressure_level", "z", "time"),
            ignore_nan = true)

Flatten `var` into a `FlatVar` according to the ordering `dims`.

The default order of dimensions for flattening is 
`("longitude", "latitude", "pressure_level", "z", "time")`. Dimensions not present in `var`
are excluded from this ordering.

If `ignore_nan = true`, then `NaNs` are excluded from the flattened data. Otherwise, `NaN`s
are included.
"""
function flatten(
    var::OutputVar;
    dims = ("longitude", "latitude", "pressure_level", "z", "time"),
    ignore_nan = true,
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

    # To flatten data, we need to permute, vectorize, and remove NaNs
    # Permute and vectorize data
    perm = Tuple(indexin(dims, collect(keys(var.dims))))
    permute_data = PermutedDimsArray(var.data, perm)
    vec_data = vec(permute_data)

    # Remove NaNs if necessary
    if ignore_nan
        nan_indices = findall(isnan, vec_data)
        if !isempty(nan_indices)
            valid_indices = setdiff(collect(eachindex(vec_data)), nan_indices)
            vec_data = @view vec_data[valid_indices]
        end
    end

    # Make metadata
    metadata = Metadata(
        deepcopy(var.attributes),
        deepcopy(var.dims),
        deepcopy(var.dim_attributes),
        dims,
        nan_indices,
    )
    return FlatVar(metadata, copy(vec_data))
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
    data_size = _data_size(metadata, ignore_nan = true)
    length(data) == data_size || error(
        "Flattened data should be of length $data_size; got length $(length(data))",
    )

    # To unflatten data, we need to remove NaNs, vectorize, and permute. Note that this is
    # the operations of flattening data, but in reverse.
    # Include NaNs
    if !isempty(metadata.nan_indices)
        data_length = values(metadata.dims) |> collect .|> length |> prod
        nan_data = fill(zero(eltype(data)), data_length)
        # NaN can be NaN32 or NaN64
        nan_data[metadata.nan_indices] .= eltype(data)(NaN)
        valid_indices = findall(!isnan, nan_data)
        nan_data[valid_indices] = data
        data = nan_data
    end

    # Reshape data
    reshape_dims = (length(metadata.dims[dim]) for dim in metadata.ordered_dims)
    data = reshape(data, reshape_dims...)

    # Permute dimensions of data
    dim2index = Dict([
        dim_name => index for
        (index, dim_name) in enumerate(keys(metadata.dims))
    ])
    perm = invperm(
        collect(dim2index[dim_name] for dim_name in metadata.ordered_dims),
    )
    data = permutedims(data, perm)

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
    _data_size(var::FlatVar; ignore_nan = true)

Get the length of `var.data`.

If `ignore_nan = true`, then the length is computed excluding `NaN`s. If
`ignore_nan = false`, then the length is computed including `NaN`s.
"""
function _data_size(var::FlatVar; ignore_nan = true)
    return _data_size(var.metadata; ignore_nan = ignore_nan)
end

"""
    _data_size(metadata::Metadata; ignore_nan = true)

Get the length of the flattened `data` according to the `metadata`.

If `ignore_nan = true`, then the length is computed excluding `NaN`s. If
`ignore_nan = false`, then the length is computed including `NaN`s.
"""
function _data_size(metadata::Metadata; ignore_nan = true)
    data_size_with_nans = prod(length.(values(metadata.dims)))
    return ignore_nan ? data_size_with_nans - length(metadata.nan_indices) :
           data_size_excluding_nans
end
