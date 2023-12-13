module CairoMakieExt

import CairoMakie
import ClimaAnalysis
import ClimaAnalysis: Visualize

"""
    contour_plot2D!(
                    fig::CairoMakie.Figure,
                    var::ClimaAnalysis.OutputVar;
                    p_loc = (1,1)
                    )

Plot a contour plot of the given 2D `var`iable in the given `fig`ure and location.

The plot comes with labels, units, and a colorbar.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)
"""
function Visualize.contour_plot2D!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
)
    length(var.dims) == 2 || error("Can only plot 2D variables")

    dim1_name, dim2_name = var.index2dim
    dim1 = var.dims[dim1_name]
    dim2 = var.dims[dim2_name]

    units = var.attributes["units"]
    short_name = var.attributes["short_name"]
    colorbar_label = "$short_name [$units]"
    dim1_units = var.dim_attributes[dim1_name]["units"]
    dim2_units = var.dim_attributes[dim2_name]["units"]

    title = var.attributes["long_name"]
    xlabel = "$dim1_name [$dim1_units]"
    ylabel = "$dim2_name [$dim2_units]"

    CairoMakie.Axis(fig[p_loc...]; title, xlabel, ylabel)

    # custom_levels is a workaround for plotting constant fields with CairoMakie
    Z = var.data
    custom_levels =
        minimum(Z) ≈ maximum(Z) ? (minimum(Z):0.1:(minimum(Z) + 0.2)) : 25
    plot = CairoMakie.contourf!(dim1, dim2, Z)

    p_loc_cb = Tuple([p_loc[1], p_loc[2] + 1])
    CairoMakie.Colorbar(fig[p_loc_cb...], plot, label = colorbar_label)
end

"""
Private function to define `sliced_` functions.

It slices a given variable and applies `func` to it.
"""
function _sliced_plot_generic(func, fig, var, cut; p_loc)
    isnothing(cut) && (cut = Dict())

    var_sliced = var

    for (dim_name, val) in cut
        var_sliced = ClimaAnalysis.slice_general(var_sliced, val, dim_name)
    end

    func(fig, var_sliced; p_loc)
end

"""
Private function to define `plot!` functions.

It composes a `cut` from given `kwargs`. Used with `sliced` functions.
"""
function _plot_generic_kwargs(func, fig, var; p_loc, kwargs...)
    cut = Dict("$k" => v for (k, v) in kwargs)
    length(cut) == 0 && (cut = nothing)
    return func(fig, var, cut; p_loc)
end

"""
    sliced_contour_plot!(
                        fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar,
                        cut::Union{Nothing, AbstractDict{String, <: Real}};
                        p_loc = (1,1)
                        )

Take a `var`iable, slice as directed, and plot a 2D contour plot in the given `fig`ure and
location.

The plot comes with labels, units, and a colorbar.

Arguments
=========

If the variable is not 2D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` contour plot at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)
"""
function Visualize.sliced_contour_plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
)
    return _sliced_plot_generic(Visualize.contour_plot2D!, fig, var, cut; p_loc)
end

"""
    contour_plot!(
                fig::CairoMakie.Figure,
                var::ClimaAnalysis.OutputVar;
                p_loc = (1,1),
                kwargs...
                )

Syntactic sugar for `sliced_contour_plot` with `kwargs` instead of `cut`.

Example
=======

`contour_plot!(fig, var, time = 100, lat = 70)` plots a contour plot by slicing `var` along
the time nearest to 100 and latitude nearest 70.

"""
function Visualize.contour_plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_contour_plot!,
        fig,
        var;
        p_loc,
        kwargs...,
    )
end

"""
    contour_plot2D!(
                    fig::CairoMakie.Figure,
                    var::ClimaAnalysis.OutputVar;
                    p_loc = (1,1)
                    )

Plot a contour plot of the given 2D `var`iable in the given `fig`ure and location.

The plot comes with labels, units, and a colorbar.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)
"""
function Visualize.line_plot1D!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
)
    length(var.dims) == 1 || error("Can only plot 1D variables")

    dim_name = var.index2dim[]
    dim = var.dims[dim_name]

    units = var.attributes["units"]
    short_name = var.attributes["short_name"]

    dim_units = var.dim_attributes[dim_name]["units"]

    title = var.attributes["long_name"]
    xlabel = "$dim_name [$dim_units]"
    ylabel = "$short_name [$units]"

    CairoMakie.Axis(fig[p_loc...]; title, xlabel, ylabel)
    CairoMakie.lines!(dim, var.data; title)
end

"""
    sliced_line_plot!(
                        fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar,
                        cut::Union{Nothing, AbstractDict{String, <: Real}};
                        p_loc = (1,1)
                        )

Take a `var`iable, slice as directed, and plot a 1D line plot in the given `fig`ure and
location.

The plot comes with labels, and units.

Arguments
=========

If the variable is not 1D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` contour plot at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)
"""
function Visualize.sliced_line_plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
)
    return _sliced_plot_generic(Visualize.line_plot1D!, fig, var, cut; p_loc)
end

"""
    line_plot!(
                fig::CairoMakie.Figure,
                var::ClimaAnalysis.OutputVar;
                p_loc = (1,1),
                kwargs...
                )

Syntactic sugar for `sliced_line_plot` with `kwargs` instead of `cut`.

Example
=======

`line_plot!(fig, var, time = 100, lat = 70)` plots a line plot by slicing `var` along
the time nearest to 100 and latitude nearest 70.

"""
function Visualize.line_plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_line_plot!,
        fig,
        var;
        p_loc,
        kwargs...,
    )
end

"""
    sliced_plot!(
                        fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar,
                        cut::Union{Nothing, AbstractDict{String, <: Real}};
                        p_loc = (1,1)
                        )

Take a `var`iable, slice as directed, and plot a 1D line plot or 2D contour plot in the given `fig`ure and
location.

The plot comes with labels, and units (and possibly a colorbar).

Arguments
=========

If the variable is not 1D/2D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` contour plot at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)
"""
function Visualize.sliced_plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
)
    initial_dim = length(var.dims)
    removed_dims = isnothing(cut) ? 0 : length(cut)
    final_dim = initial_dim - removed_dims

    if final_dim == 1
        fun = Visualize.line_plot1D!
    elseif final_dim == 2
        fun = Visualize.contour_plot2D!
    else
        error("Sliced variable has $dim_num dimensions (needed 1 or 2)")
    end

    return _sliced_plot_generic(fun, fig, var, cut; p_loc)
end


"""
    plot!(
          fig::CairoMakie.Figure,
          var::ClimaAnalysis.OutputVar;
          p_loc = (1,1),
          kwargs...
          )

Syntactic sugar for `sliced_plot` with `kwargs` instead of `cut`.

Example
=======

`line_plot!(fig, var, time = 100, lat = 70)` plots a line plot or a contour plot by slicing
`var` along the time nearest to 100 and latitude nearest 70.

"""
function Visualize.plot!(
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    kwargs...,
)
    _plot_generic_kwargs(Visualize.sliced_plot!, fig, var; p_loc, kwargs...)
end

end
