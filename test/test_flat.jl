using Test
import ClimaAnalysis

import ClimaAnalysis.Template:
    TemplateVar,
    make_template_var,
    add_attribs,
    add_dim,
    add_time_dim,
    add_lon_dim,
    add_lat_dim,
    add_data,
    ones_data,
    zeros_data,
    one_to_n_data,
    initialize

@testset "FlatVar" begin
    lat = [-90.0, -30.0, 30.0, 90.0]
    lon = [-60.0, -30.0, 0.0, 30.0, 60.0]
    time = [0.0, 1.0, 5.0]
    n_elts = length(lat) * length(lon) * length(time)
    size_of_data = (length(lat), length(lon), length(time))
    data = reshape(collect(1.0:n_elts), size_of_data...)
    var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", time, units = "seconds") |>
        add_attribs(long_name = "hi") |>
        add_data(data = data) |>
        initialize

    flat_var = ClimaAnalysis.flatten(var, dims = ("latitude", "longitude", "t"))
    @test ClimaAnalysis.flatten_dim_order(flat_var) == ("lat", "lon", "time")
    @test ClimaAnalysis.flattened_length(flat_var) == n_elts
    @test flat_var.data == vec(data)
    @test flat_var.metadata.attributes == var.attributes
    @test flat_var.metadata.dims == var.dims
    @test flat_var.metadata.dim_attributes == var.dim_attributes
    @test flat_var.metadata.ordered_dims == ("lat", "lon", "time")

    unflatten_var = ClimaAnalysis.unflatten(flat_var)
    @test var.data == unflatten_var.data
    @test var.attributes == unflatten_var.attributes
    @test var.dim_attributes == unflatten_var.dim_attributes
    @test var.dims == unflatten_var.dims

    flat_var = ClimaAnalysis.Var.flatten(var)
    @test ClimaAnalysis.Var.flatten_dim_order(flat_var) ==
          ("lon", "lat", "time")
    @test flat_var.data == vec(permutedims(data, (2, 1, 3)))
    @test flat_var.metadata.attributes == var.attributes
    @test flat_var.metadata.dims == var.dims
    @test flat_var.metadata.dim_attributes == var.dim_attributes
    @test flat_var.metadata.ordered_dims == ("lon", "lat", "time")

    unflatten_var = ClimaAnalysis.Var.unflatten(flat_var)
    @test var.data == unflatten_var.data
    @test var.attributes == unflatten_var.attributes
    @test var.dim_attributes == unflatten_var.dim_attributes
    @test var.dims == unflatten_var.dims

    # Test unflattening with different datas
    flat_var = ClimaAnalysis.Var.flatten(var)
    var2 = ClimaAnalysis.Var.remake(
        var;
        data = reshape(collect((-n_elts):-1), size_of_data...),
    )
    var2 = permutedims(var2, ("lon", "time", "lat"))
    flat_var2 = ClimaAnalysis.Var.flatten(var2)
    reconstructed_var =
        ClimaAnalysis.Var.unflatten(flat_var.metadata, flat_var2.data)
    reconstructed_var = permutedims(reconstructed_var, ("lon", "time", "lat"))
    @test reconstructed_var.data == var2.data
    @test reconstructed_var.attributes == var.attributes
    @test reconstructed_var.dim_attributes == var.dim_attributes
    @test reconstructed_var.dims == var.dims

    # Test with ignore_nan = false
    flat_var = ClimaAnalysis.Var.flatten(var2, ignore_nan = false)
    reconstructed_var =
        ClimaAnalysis.unflatten(flat_var.metadata, flat_var.data)
    @test isequal(reconstructed_var.data, var2.data)
    @test reconstructed_var.attributes == var2.attributes
    @test reconstructed_var.dim_attributes == var2.dim_attributes
    @test reconstructed_var.dims == var2.dims

    # Test with NaNs and ignore_nan = true
    nan_data = reshape(collect(1.0:n_elts), size_of_data...)
    nan_data[1, 1, 1] = NaN
    nan_data[3, 2, 3] = NaN
    nan_var = ClimaAnalysis.remake(var, data = nan_data)
    flat_var = ClimaAnalysis.Var.flatten(nan_var)
    @test ClimaAnalysis.Var.flatten_dim_order(flat_var) ==
          ("lon", "lat", "time")
    @test ClimaAnalysis.flattened_length(flat_var) == n_elts - 2
    @test flat_var.data == filter(!isnan, vec(permutedims(nan_data, (2, 1, 3))))
    @test flat_var.metadata.attributes == nan_var.attributes
    @test flat_var.metadata.dims == nan_var.dims
    @test flat_var.metadata.dim_attributes == nan_var.dim_attributes
    @test flat_var.metadata.ordered_dims == ("lon", "lat", "time")

    unflatten_var = ClimaAnalysis.Var.unflatten(flat_var)
    @test isequal(nan_var.data, unflatten_var.data)
    @test nan_var.attributes == unflatten_var.attributes
    @test nan_var.dim_attributes == unflatten_var.dim_attributes
    @test nan_var.dims == unflatten_var.dims

    # Test unflattening with different NaN datas
    nan_data2 = reshape(collect(Float64.((-n_elts):-1)), size_of_data...)
    nan_data2[1, 1, 1] = NaN
    nan_data2[3, 2, 3] = NaN
    nan_var2 = ClimaAnalysis.remake(var, data = nan_data2)
    nan_var2 = permutedims(nan_var2, ("lon", "time", "lat"))
    flat_var2 = ClimaAnalysis.flatten(nan_var2)
    # Reconstruct with different data
    reconstructed_var =
        ClimaAnalysis.Var.unflatten(flat_var.metadata, flat_var2.data)
    reconstructed_var = permutedims(reconstructed_var, ("lon", "time", "lat"))
    @test isequal(reconstructed_var.data, nan_var2.data)
    @test reconstructed_var.attributes == var.attributes
    @test reconstructed_var.dim_attributes == var.dim_attributes
    @test reconstructed_var.dims == var.dims

    # Test with masking and ignore_nan = false
    mask_var = ClimaAnalysis.slice(var, t = 0.0)
    replace!(x -> 1.0, mask_var)
    mask_var = permutedims(mask_var, ("lon", "lat"))
    mask_var.data[:, 1] .= 0.0
    mask_var =
        ClimaAnalysis.generate_lonlat_mask(mask_var, NaN, 1.0, threshold = 0.5)

    flat_var = ClimaAnalysis.flatten(var, ignore_nan = false, mask = mask_var)
    @test ClimaAnalysis.flattened_length(flat_var) == 45

    permuted_var = permutedims(var, ("lon", "lat", "time"))
    data = copy(permuted_var.data)
    data[:, 1, :] .= NaN
    @test flat_var.data == filter(!isnan, vec(data))

    reconstructed_var = ClimaAnalysis.unflatten(flat_var)
    @test isequal(var.data, reconstructed_var.data)
    @test var.dims == reconstructed_var.dims
    @test var.dim_attributes == reconstructed_var.dim_attributes
    @test var.attributes == reconstructed_var.attributes

    # Test comparing apply_oceanmask and ignore_nan
    ocean_mask = ClimaAnalysis.generate_ocean_mask(NaN, 1.0)
    masked_var = ocean_mask(var)

    flat_var_mask = ClimaAnalysis.flatten(var, mask = ocean_mask)
    flat_var_ignore = ClimaAnalysis.flatten(masked_var, ignore_nan = true)

    @test flat_var_mask.data == flat_var_ignore.data
    @test flat_var_mask.metadata.drop_mask == flat_var_ignore.metadata.drop_mask
    @test flat_var.metadata.dims == flat_var_ignore.metadata.dims
    @test flat_var.metadata.dim_attributes ==
          flat_var_ignore.metadata.dim_attributes
    @test flat_var.metadata.attributes == flat_var_ignore.metadata.attributes

    reconstructed_var_mask = ClimaAnalysis.unflatten(flat_var_mask)
    reconstructed_var_ignore = ClimaAnalysis.unflatten(flat_var_ignore)
    @test reconstructed_var_mask.attributes ==
          reconstructed_var_ignore.attributes
    @test reconstructed_var_mask.dims == reconstructed_var_ignore.dims
    @test reconstructed_var_mask.dim_attributes ==
          reconstructed_var_ignore.dim_attributes

    # Test with ignore_nan = true and masking
    nan_data = copy(var.data)
    nan_data[:, 1, 1] .= NaN
    nan_var3 = ClimaAnalysis.remake(var, data = nan_data)
    flat_var =
        ClimaAnalysis.flatten(nan_var3, ignore_nan = true, mask = mask_var)
    @test ClimaAnalysis.flattened_length(flat_var) == 42

    permuted_var = permutedims(var, ("lon", "lat", "time"))
    data = copy(permuted_var.data)
    data[:, 1, :] .= NaN
    data[1, :, 1] .= NaN
    @test flat_var.data == filter(!isnan, vec(data))

    reconstructed_var = ClimaAnalysis.unflatten(flat_var)
    @test isequal(nan_var3.data, reconstructed_var.data)
    @test nan_var3.dims == reconstructed_var.dims
    @test nan_var3.dim_attributes == reconstructed_var.dim_attributes
    @test nan_var3.attributes == reconstructed_var.attributes

    # Test with ignore_nan = false and masking
    mask_var = ClimaAnalysis.slice(var, t = 0.0)
    replace!(x -> 1.0, mask_var)
    mask_var = permutedims(mask_var, ("lon", "lat"))
    mask_var.data[:, 1] .= 0.0
    mask_var =
        ClimaAnalysis.generate_lonlat_mask(mask_var, NaN, 1.0, threshold = 0.5)

    flat_var = ClimaAnalysis.flatten(var, ignore_nan = false, mask = mask_var)
    @test ClimaAnalysis.flattened_length(flat_var) == 45

    permuted_var = permutedims(var, ("lon", "lat", "time"))
    data = copy(permuted_var.data)
    data[:, 1, :] .= NaN
    @test flat_var.data == filter(!isnan, vec(data))

    reconstructed_var = ClimaAnalysis.unflatten(flat_var)
    @test isequal(var.data, reconstructed_var.data)
    @test var.dims == reconstructed_var.dims
    @test var.dim_attributes == reconstructed_var.dim_attributes
    @test var.attributes == reconstructed_var.attributes

    # Test with ignore_nan = true and masking
    nan_data = copy(var.data)
    nan_data[:, 1, 1] .= NaN
    nan_var3 = ClimaAnalysis.remake(var, data = nan_data)
    flat_var =
        ClimaAnalysis.flatten(nan_var3, ignore_nan = true, mask = mask_var)
    @test ClimaAnalysis.flattened_length(flat_var) == 42

    permuted_var = permutedims(var, ("lon", "lat", "time"))
    data = copy(permuted_var.data)
    data[:, 1, :] .= NaN
    data[1, :, 1] .= NaN
    @test flat_var.data == filter(!isnan, vec(data))

    reconstructed_var = ClimaAnalysis.unflatten(flat_var)
    @test isequal(nan_var3.data, reconstructed_var.data)
    @test nan_var3.dims == reconstructed_var.dims
    @test nan_var3.dim_attributes == reconstructed_var.dim_attributes
    @test nan_var3.attributes == reconstructed_var.attributes

    # Missing dimensions
    @test_throws ErrorException ClimaAnalysis.Var.flatten(
        var,
        dims = ("lat", "lon"),
    )

    # Length of data is not right
    @test_throws ErrorException ClimaAnalysis.Var.unflatten(
        flat_var.metadata,
        [1.0, 2.0, 3.0],
    )
end

@testset "Flatten with metadata" begin
    lat = [-90.0, -30.0, 30.0, 90.0]
    lon = [-60.0, -30.0, 0.0, 30.0, 60.0]
    time = [0.0, 1.0, 5.0]
    n_elts = length(lat) * length(lon) * length(time)
    size_of_data = (length(lat), length(lon), length(time))
    data = reshape(collect(1.0:n_elts), size_of_data...)
    start_date_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", time, units = "seconds") |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = data) |>
        initialize

    flat_var = ClimaAnalysis.flatten(start_date_var)

    # Number of dimensions are not the same
    less_dim_var =
        TemplateVar() |>
        add_dim("time", time, units = "seconds") |>
        add_attribs(long_name = "hi") |>
        initialize
    @test_throws ErrorException ClimaAnalysis.flatten(less_dim_var, flat_var)

    # Start date present in one and not the other
    no_start_date_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", time, units = "seconds") |>
        add_attribs(long_name = "hello") |>
        initialize
    @test_throws ErrorException ClimaAnalysis.flatten(
        no_start_date_var,
        flat_var,
    )

    # Absolute times not the same with the same start dates
    diff_start_date_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", [0.0, 1.0, 2.0], units = "seconds") |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = data) |>
        initialize

    @test_throws ErrorException ClimaAnalysis.flatten(
        diff_start_date_var,
        flat_var,
    )

    # Absolute times not the same with no start dates
    no_start_date_var2 =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", [0.0, 1.0, 2.0], units = "seconds") |>
        add_attribs(long_name = "hi") |>
        add_data(data = data) |>
        initialize
    no_start_date_metadata2 = ClimaAnalysis.flatten(no_start_date_var2).metadata
    @test_throws ErrorException ClimaAnalysis.flatten(
        no_start_date_var,
        no_start_date_metadata2,
    )

    # Units of dimensions do not match
    radian_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "radians") |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = lat) |>
        initialize
    degrees_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = lat) |>
        initialize
    @test_throws ErrorException ClimaAnalysis.flatten(
        radian_var,
        ClimaAnalysis.flatten(degrees_var),
    )

    # Units are missing for both
    missing_units_var =
        TemplateVar() |>
        add_dim("lat", lat) |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = lat) |>
        initialize
    @test_logs (:warn, r"missing in") ClimaAnalysis.flatten(
        missing_units_var,
        ClimaAnalysis.flatten(missing_units_var),
    )

    # Good cases
    function test_flatten(src_var, target_var, perms)
        for perm in perms
            flat_var = ClimaAnalysis.flatten(permutedims(src_var, perm))
            for perm in perms
                flat_target_var = ClimaAnalysis.flatten(
                    permutedims(target_var, perm),
                    flat_var,
                )
                @test flat_target_var.data == flat_var.data
                @test ClimaAnalysis.flatten_dim_order(flat_target_var) ==
                      ClimaAnalysis.flatten_dim_order(flat_var)
            end
        end
    end
    perms_3d = (
        ("lon", "lat", "time"),
        ("lat", "lon", "time"),
        ("lat", "lon", "time"),
        ("lat", "time", "lon"),
        ("time", "lon", "lat"),
        ("time", "lat", "lon"),
    )

    # Flatten with the OutputVar and metadata created from flatten for a OutputVar with
    # a start date and no start date
    test_flatten(start_date_var, start_date_var, perms_3d)
    test_flatten(no_start_date_var, no_start_date_var, perms_3d)

    # Absolute times are different with different start dates
    different_start_date_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_dim("time", [-5.0, -4.0, 0.0], units = "seconds") |>
        add_attribs(long_name = "hi", start_date = "2010-12-01T00:00:05") |>
        add_data(data = data) |>
        initialize
    test_flatten(start_date_var, different_start_date_var, perms_3d)
    test_flatten(different_start_date_var, start_date_var, perms_3d)

    # No time dimension but one OutputVar has a start date and the other one does not
    start_date_no_time_var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees") |>
        add_dim("lon", lon, units = "degrees") |>
        add_attribs(long_name = "hi", start_date = "2010-12-1") |>
        add_data(data = data[:, :, 1]) |>
        initialize

    attribs = Dict("long_name" => "hi")
    no_start_date_no_time_var =
        ClimaAnalysis.remake(start_date_no_time_var, attributes = attribs)
    perms_2d = (("lon", "lat"), ("lat", "lon"))
    test_flatten(start_date_no_time_var, no_start_date_no_time_var, perms_2d)
    test_flatten(no_start_date_no_time_var, start_date_no_time_var, perms_2d)

    # Flatten with dropping NaNs
    nan_data = reshape(collect(1.0:n_elts), size_of_data...)
    nan_data[1, 1, 1] = NaN
    nan_data[1, 2, 3] = NaN
    nan_data[2, 3, 2] = NaN
    nan_data[4, 5, 1] = NaN
    nan_var = ClimaAnalysis.remake(start_date_var, data = nan_data)
    test_flatten(nan_var, nan_var, perms_3d)
end

@testset "Extracting dimensions, units, and names for FlatVar" begin
    lat = collect(range(-89.5, 89.5, 3))
    lon = collect(range(-179.5, 179.5, 4))
    time = [0.0, 1.0, 2.0]
    pfull = [0.0, 2.0, 4.0, 6.0]
    z = [0.0, 5.0]
    var =
        TemplateVar() |>
        add_dim("lat", lat, units = "degrees_north") |>
        add_dim("lon", lon, units = "degrees_east") |>
        add_dim("time", time, units = "s") |>
        add_dim("pfull", pfull, units = "Pa") |>
        add_dim("z", z, units = "m") |>
        add_attribs(
            short_name = "hi",
            long_name = "hello",
            units = "idk",
            start_date = "1979-1-1",
        ) |>
        initialize

    flat_var = ClimaAnalysis.flatten(var)

    @test ClimaAnalysis.has_latitude(flat_var)
    @test ClimaAnalysis.has_longitude(flat_var)
    @test ClimaAnalysis.has_time(flat_var)
    @test ClimaAnalysis.has_pressure(flat_var)
    @test ClimaAnalysis.has_altitude(flat_var)
    @test !ClimaAnalysis.has_date(flat_var)
    @test_throws ErrorException ClimaAnalysis.date_name(flat_var)
    @test ClimaAnalysis.latitudes(flat_var) == lat
    @test ClimaAnalysis.longitudes(flat_var) == lon
    @test ClimaAnalysis.times(flat_var) == time
    @test ClimaAnalysis.pressures(flat_var) == pfull
    @test ClimaAnalysis.altitudes(flat_var) == z
    @test ClimaAnalysis.dates(flat_var) == ClimaAnalysis.dates(var)
    @test ClimaAnalysis.has_units(flat_var)
    @test ClimaAnalysis.units(flat_var) == "idk"
    @test ClimaAnalysis.dim_units(flat_var, "lon") == "degrees_east"
    @test ClimaAnalysis.short_name(flat_var) == "hi"
    @test ClimaAnalysis.long_name(flat_var) == "hello"

    metadata = flat_var.metadata
    @test ClimaAnalysis.has_latitude(metadata)
    @test ClimaAnalysis.has_longitude(metadata)
    @test ClimaAnalysis.has_time(metadata)
    @test ClimaAnalysis.has_pressure(metadata)
    @test ClimaAnalysis.has_altitude(metadata)
    @test !ClimaAnalysis.has_date(metadata)
    @test_throws ErrorException ClimaAnalysis.date_name(metadata)
    @test ClimaAnalysis.latitudes(metadata) == lat
    @test ClimaAnalysis.longitudes(metadata) == lon
    @test ClimaAnalysis.times(metadata) == time
    @test ClimaAnalysis.pressures(metadata) == pfull
    @test ClimaAnalysis.altitudes(metadata) == z
    @test ClimaAnalysis.dates(metadata) == ClimaAnalysis.dates(var)
    @test ClimaAnalysis.has_units(metadata)
    @test ClimaAnalysis.units(metadata) == "idk"
    @test ClimaAnalysis.dim_units(metadata, "lon") == "degrees_east"
    @test ClimaAnalysis.short_name(metadata) == "hi"
    @test ClimaAnalysis.long_name(metadata) == "hello"
end
