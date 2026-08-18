"""
    group_and_reduce_by(var::OutputVar, group_by, reduce_by)

Group the dimension `dim_name` using `group_by` and apply the reduction
`reduce_by` along the dimension.

The first element in each group is used as the dimension of the resulting `OutputVar`.

Only the short name, long name, units, and start date are kept. All other attributes are
discarded in the process.

The function `group_by` takes in a vector of values for a dimension and returns a vector of
vectors of values of the dimension. The function `reduce_by` must be of the form `f(A;
dims)`, where `A` is a slice of `var.data` and `dims` is the index of `A` to reduce against
and return the reduction over the `dims` dimension.

Note `group_by` does not need to partition the values of the dimension. For example, the
`group_by` function can be `group_by(dim) = [[first(dim)]]`, which return the first slice of
the `OutputVar`.
"""
function group_and_reduce_by(var::OutputVar, dim_name, group_by, reduce_by)
    dim_name_in_var =
        ClimaAnalysis.Var.find_corresponding_dim_name_in_var(dim_name, var)
    dim_idx = var.dim2index[dim_name_in_var]

    # Group by
    dim_vals_groups = group_by(var.dims[dim_name_in_var])
    dim_indices_groups =
        indexin.(dim_vals_groups, Ref(var.dims[dim_name_in_var]))
    data_groups = map(dim_indices_groups) do dim_indices
        index_tuple = ntuple(
            idx -> idx == dim_idx ? dim_indices : Colon(),
            ndims(var.data),
        )
        view(var.data, index_tuple...)
    end

    # Reduce by and concat
    ret_data = cat(reduce_by.(data_groups, dims = dim_idx)..., dims = dim_idx)

    # Get the elements for constructing the new dimension
    ret_dim_indices = [
        first(dim_indices) for
        dim_indices in dim_indices_groups if !isempty(dim_indices)
    ]

    # New dimension to return
    dim = var.dims[dim_name_in_var][ret_dim_indices]

    # Make new dimensions for OutputVar
    ret_dims = deepcopy(var.dims)
    ret_dims[dim_name_in_var] = dim

    # Keep short name, long name, units, and start_date and discard the rest
    keep_attribs = ("long_name", "short_name", "units", "start_date")
    ret_attribs = Dict(
        attrib => var.attributes[attrib] for
        attrib in keep_attribs if attrib in keys(var.attributes)
    )
    return ClimaAnalysis.remake(
        var,
        attributes = ret_attribs,
        data = ret_data,
        dims = ret_dims,
    )
end
