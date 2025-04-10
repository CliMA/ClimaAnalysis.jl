import NCDatasets
import OrderedCollections: OrderedDict
import Artifacts

"""
The variable `LAND_MASK` is an `OutputVar` whose data consists of ones and NaNs. The ones
indicate ocean and NaNs indicate the land.
"""
const LAND_MASK = let
    path = joinpath(
        Artifacts.artifact"landsea_mask_30arcseconds",
        "landsea_mask.nc",
    )
    ds = NCDatasets.NCDataset(path)
    # Taking every 16th element corresponds to roughly 0.1 degree resolution
    # This is roughly 0.05 GB in size
    # Every 8th element -> ~0.25 GB
    # Every 4th element -> ~1 GB
    lon = ds["lon"][begin:16:end]
    lat = ds["lat"][begin:16:end]
    data = ds["landsea"][begin:16:end, begin:16:end]
    close(ds)

    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "Land sea mask")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = OutputVar(attribs, dims, dim_attribs, data)
    var = replace(var, 0.0 => 1.0, 1.0 => NaN)
end

"""
The variable `OCEAN_MASK` is an `OutputVar` whose data consists of ones and NaNs. The ones
indicate land and NaNs indicate the ocean.
"""
const OCEAN_MASK = let
    path = joinpath(
        Artifacts.artifact"landsea_mask_30arcseconds",
        "landsea_mask.nc",
    )
    ds = NCDatasets.NCDataset(path)
    # See comment in LAND_MASK for why 16 is chosen
    lon = ds["lon"][begin:16:end]
    lat = ds["lat"][begin:16:end]
    data = ds["landsea"][begin:16:end, begin:16:end]
    close(ds)

    dims = OrderedDict(["lon" => lon, "lat" => lat])
    attribs = Dict("long_name" => "Land sea mask")
    dim_attribs = OrderedDict([
        "lat" => Dict("units" => "deg"),
        "lon" => Dict("units" => "deg"),
    ])
    var = OutputVar(attribs, dims, dim_attribs, data)
    var = replace(var, 0.0 => NaN)
end
