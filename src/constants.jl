"""
The variable LAND_MASK is a filepath to a NCDataset whose dimensions are latitude and
longitude and data are ones and zeros. The zeros indicate land and ones indicate everything
else.

The land mask is generated from the ETOPO2022 dataset at 60 arc-second resolutions
(https://www.ncei.noaa.gov/products/etopo-global-relief-model). The dimensions latitude and
longitude are sampled every 16th element to thin out the dimensions. Any value for z greater
than 0 is assigned 0.0 and all other values for z is assigned 1.0.
"""
LAND_MASK = joinpath(@__DIR__, "..", "masks", "land_mask16.nc")

"""
The variable OCEAN_MASK is a filepath to a NCDataset whose dimensions are latitude and
longitude and data are 1s and 0s. The zeros indicate ocean and ones indicate everything
else.

The land mask is generated from the ETOPO2022 dataset at 60 arc-second resolutions
(https://www.ncei.noaa.gov/products/etopo-global-relief-model). The dimensions latitude and
longitude are sampled every 16th element to thin out the dimensions. Any value for z less
than or equal to 0 is assigned 0.0 and all other values for z is assigned 1.0.
"""
OCEAN_MASK = joinpath(@__DIR__, "..", "masks", "ocean_mask16.nc")
