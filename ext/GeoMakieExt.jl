module GeoMakieExt

import GeoMakie
import GeoMakie: Makie
import ClimaAnalysis
import ClimaAnalysis: Visualize

MakiePlace = Union{Makie.Figure, Makie.GridLayout}

const LONGITUDE_NAMES = Set(["lon", "long"])
const LATITUDE_NAMES = Set(["lat"])

"""
    heatmap2D_on_globe!(fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar;
                        p_loc = (1,1),
                        plot_coastline = true,
                        more_kwargs)
    heatmap2D_on_globe!(grid_layout::CairoMakie.GridLayout,
                        var::ClimaAnalysis.OutputVar;
                        p_loc = (1,1),
                        plot_coastline = true,
                        more_kwargs)


Plot a heatmap of the given 2D `var`iable on a projected geoid.

The plot comes with labels, units, and a colorbar.

This function assumes that the following attributes are available:
- long_name
- short_name
- units

The dimensions have to be longitude and latitude.

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)
- the coastline (`:coast`)

The coastline is plotted from `GeoMakie.coastline` using the `lines!` plotting function.

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.

"""
function Visualize.heatmap2D_on_globe!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    plot_coastline = true,
    more_kwargs = Dict(
        :plot => Dict(),
        :cb => Dict(),
        :axis => Dict(),
        :coast => Dict(:color => :black),
    ),
)
    length(var.dims) == 2 || error("Can only plot 2D variables")

    lon_name, lat_name = var.index2dim
    lon_name in LONGITUDE_NAMES || error("$lon_name is not a longitude")
    lat_name in LATITUDE_NAMES || error("$lat_name is not a latitude")
    lon = var.dims[lon_name]
    lat = var.dims[lat_name]

    units = var.attributes["units"]
    short_name = var.attributes["short_name"]
    colorbar_label = "$short_name [$units]"

    axis_kwargs = get(more_kwargs, :axis, Dict())
    plot_kwargs = get(more_kwargs, :plot, Dict())
    cb_kwargs = get(more_kwargs, :cb, Dict())
    coast_kwargs = get(more_kwargs, :coast, Dict(:color => :black))

    title = get(axis_kwargs, :title, var.attributes["long_name"])

    GeoMakie.GeoAxis(place[p_loc...]; title, axis_kwargs...)

    plot = Makie.surface!(lon, lat, var.data; plot_kwargs...)
    plot_coastline && Makie.lines!(GeoMakie.coastlines(); coast_kwargs...)

    p_loc_cb = Tuple([p_loc[1], p_loc[2] + 1])
    Makie.Colorbar(
        place[p_loc_cb...],
        plot,
        label = colorbar_label;
        cb_kwargs...,
    )

end

end
