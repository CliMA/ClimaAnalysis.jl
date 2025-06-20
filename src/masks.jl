export apply_landmask, apply_oceanmask, make_lonlat_mask, generate_lonlat_mask

import NCDatasets
import OrderedCollections: OrderedDict
import Artifacts

"""
The variable `LANDSEA_MASK` is an `OutputVar` whose data consists of zeros and ones. The
ones indicate land and zeros indicate ocean.
"""
const LANDSEA_MASK = let
    path = joinpath(
        Artifacts.artifact"landsea_mask_30arcseconds",
        "landsea_mask.nc",
    )
    lon, lat, data = NCDatasets.NCDataset(path) do ds

        # Only get every third element because ClimaAnalysis.Utils._isequispaced(lat) failed
        # for ds["lat"][:] and ds["lat"][begin:2:end]
        lon = ds["lon"][begin:3:end]
        lat = ds["lat"][begin:3:end]

        # Cast to a BitMatrix to further reduce memory footprint. Note that the values of
        # ds["landsea"][:,:] are stored as Int8. There are 14400 elements for the longitude
        # dimension and 7200 elements in the latitude dimensions if every third element is
        # loaded. Hence, the memory footprint of this is about 98.88 MB before casting to a
        # BitMatrix.
        data = Bool.(ds["landsea"][begin:3:end, begin:3:end])
        lon, lat, data
    end

    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "Land sea mask")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = OutputVar(attribs, dims, dim_attribs, data)
end

"""
    generate_lonlat_mask(var::OutputVar, zero_to, one_to; threshold = 0.5)

Returns a masking function that takes an `OutputVar` and masks its data using `var.data`,
which contains only zeros and ones.

The data in `var` is resampled to match the dimensions of the `OutputVar`. After resampling,
values greater than or equal to `threshold` are rounded to one and all other values are
rounded to zero. Zeros are replaced with `zero_to`, and ones are replaced with `one_to`. The
masking is applied as element wise multiplication.
"""
function generate_lonlat_mask(
    mask_var::OutputVar,
    zero_to,
    one_to;
    threshold = 0.5,
)
    # Check mask is binary
    all(x -> iszero(x) || isone(x), mask_var.data) || error(
        "var is not a OutputVar whose data contains only zeros and ones. Use replace or replace! to to ensure that all values in the data are zeros and ones",
    )

    # Check if lon and lat are the only dimensions
    has_longitude(mask_var) ||
        error("Longitude does not exist as a dimension in var")
    has_latitude(mask_var) ||
        error("Latitude does not exist as a dimension in var")
    length(mask_var.dims) == 2 ||
        error("Number of dimensions ($(length(mask_var.dims))) is not two")

    function apply_lonlat_mask(var)
        # Check if longitude and latitude exist in var
        has_longitude(var) || error("var does not has a longitude dimension")
        has_latitude(var) || error("var does not has a latitude dimension")
        # Use _resampled_as_partial as we do not want to do units checking as
        # it would be too restrictive
        resampled_mask_var =
            _resampled_as_partial(mask_var, var, ("lon", "lat"))
        # Apply threshold and zero_to and one_to in a single step
        mask = replace(resampled_mask_var.data) do val
            val >= threshold ? one_to : zero_to
        end

        # Reshape data for broadcasting
        lon_idx = var.dim2index[longitude_name(var)]
        lat_idx = var.dim2index[latitude_name(var)]
        lon_length = var.dims[longitude_name(var)] |> length
        lat_length = var.dims[latitude_name(var)] |> length
        if lon_idx > lat_idx
            mask = transpose(mask)
        end
        size_to_reshape = (
            if i == lon_idx
                lon_length
            elseif i == lat_idx
                lat_length
            else
                1
            end for i in 1:ndims(var.data)
        )
        mask = reshape(mask, size_to_reshape...)

        # Apply mask
        data = copy(var.data)
        data .*= mask

        return remake(var, data = data)
    end
    return apply_lonlat_mask
end

"""
    apply_landmask(var::OutputVar)

Apply a land mask to `var` by NaNing any data whose coordinates are located on land.

!!! note "Threshold keyword argument"
    The `threshold` keyword argument is available in ClimaAnalysis v0.5.18 and beyond.
"""
function apply_landmask(var::OutputVar; threshold = 0.5)
    mask_fn =
        generate_lonlat_mask(LANDSEA_MASK, 1.0, NaN, threshold = threshold)
    return mask_fn(var)
end

"""
    apply_oceanmask(var::OutputVar)

Apply an ocean mask to `var` by NaNing any data whose coordinates are in the ocean.

!!! note "Threshold keyword argument"
    The `threshold` keyword argument is available in ClimaAnalysis v0.5.18 and beyond.
"""
function apply_oceanmask(var::OutputVar, threshold = 0.5)
    mask_fn =
        generate_lonlat_mask(LANDSEA_MASK, NaN, 1.0, threshold = threshold)
    return mask_fn(var)
end

"""
    make_lonlat_mask(var;
                     set_to_val = nothing,
                     true_val = NaN,
                     false_val = 1.0)

Return a masking function that takes in a `OutputVar` and mask the data using `var.data`.

The parameter `set_to_val` is a function that takes in an element of `var.data` and return
`true` or `false`. If `set_to_val` returns `true`, then the value will be `true_val` in the
mask and if `set_to_nan` returns `false`, then the value will be `false_val` in the mask.

If `set_to_val` is `nothing`, then no transformation is done and `var.data` is used. This is
helpful if `var.data` is already an array of NaNs and ones or an array of zeros and ones.

!!! warning "Deprecated"
    This function is deprecated and users are encouraged to use
    [`generate_lonlat_mask`](@ref) instead.
"""
function make_lonlat_mask(
    var;
    set_to_val = nothing,
    true_val = NaN,
    false_val = 1.0,
)
    Base.depwarn(
        "This function is deprecated. Users are encouraged to use `generate_lonlat_mask` instead.",
        :make_lonlat_mask!,
    )
    # Check if lon and lat are the only dimensions
    has_longitude(var) ||
        error("Longitude does not exist as a dimension in var")
    has_latitude(var) || error("Latitude does not exist as a dimension in var")
    length(var.dims) == 2 ||
        error("Number of dimensions ($(length(var.dims))) is not two")

    # Preprocess data for the mask if needed
    if !isnothing(set_to_val)
        true_indices = findall(set_to_val, var.data)
        false_indices = findall(!set_to_val, var.data)
        cleaned_up_data = var.data |> copy
        cleaned_up_data[true_indices] .= true_val
        cleaned_up_data[false_indices] .= false_val
        var = OutputVar(
            var.attributes |> deepcopy,
            var.dims |> deepcopy,
            var.dim_attributes |> deepcopy,
            cleaned_up_data,
        )
    end

    mask_var = permutedims(var, ("lon", "lat"))

    return apply_lonlat_mask(input_var) = begin
        # Check if lon and lat exist
        has_longitude(input_var) ||
        error("Longitude does not exist as a dimension in var")
        has_latitude(input_var) ||
        error("Latitude does not exist as a dimension in var")

        # Resample so that the mask match up with the grid of var
        # Round because linear resampling is done and we want the mask to be only ones and zeros
        intp = _make_interpolant(mask_var.dims, mask_var.data)
        mask_arr =
            [
                intp(pt...) for pt in Base.product(
                    input_var.dims[longitude_name(input_var)],
                    input_var.dims[latitude_name(input_var)],
                )
            ] .|> round

        # Reshape data for broadcasting
        lon_idx = input_var.dim2index[longitude_name(input_var)]
        lat_idx = input_var.dim2index[latitude_name(input_var)]
        lon_length = input_var.dims[longitude_name(input_var)] |> length
        lat_length = input_var.dims[latitude_name(input_var)] |> length
        if lon_idx > lat_idx
            mask_arr = transpose(mask_arr)
        end
        size_to_reshape = (
            if i == lon_idx
                lon_length
            elseif i == lat_idx
                lat_length
            else
                1
            end for i in 1:ndims(input_var.data)
        )

        # Mask data
        mask_arr = reshape(mask_arr, size_to_reshape...)
        masked_data = input_var.data .* mask_arr

        return remake(input_var, data = masked_data)
    end
end
