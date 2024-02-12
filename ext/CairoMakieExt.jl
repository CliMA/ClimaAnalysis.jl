module CairoMakieExt

import CairoMakie
import ClimaAnalysis
import ClimaAnalysis: Visualize

MakiePlace = Union{CairoMakie.Figure, CairoMakie.GridLayout}

"""
    heatmap2D!(fig::CairoMakie.Figure,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs)
    heatmap2D!(grid_layout::CairoMakie.GridLayout,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs)


Plot a heatmap of the given 2D `var`iable in the given place and location.
The place can be a `Figure` or a `GridLayout`.

The plot comes with labels, units, and a colorbar.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.

"""
function Visualize.heatmap2D!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
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

    axis_kwargs = get(more_kwargs, :axis, Dict())
    plot_kwargs = get(more_kwargs, :plot, Dict())
    cb_kwargs = get(more_kwargs, :cb, Dict())

    title = get(axis_kwargs, :title, var.attributes["long_name"])
    xlabel = get(axis_kwargs, :xlabel, "$dim1_name [$dim1_units]")
    ylabel = get(axis_kwargs, :ylabel, "$dim2_name [$dim2_units]")

    # dim_on_y is only supported by plot_line1D. We remove it here to ensure that we can a
    # consistent entry point between plot_line1D and heatmap2D. It we left it here, it would
    # be passed down and lead to a unknown argument error.
    #
    # TODO: Refactor: We shouldn't have to deal with dim_on_y if we don't use it!
    if haskey(axis_kwargs, :dim_on_y)
        axis_kwargs_dict = Dict(axis_kwargs)
        pop!(axis_kwargs_dict, :dim_on_y)
        axis_kwargs = pairs(axis_kwargs_dict)
    end

    CairoMakie.Axis(place[p_loc...]; title, xlabel, ylabel, axis_kwargs...)

    plot = CairoMakie.heatmap!(dim1, dim2, var.data; plot_kwargs...)

    p_loc_cb = Tuple([p_loc[1], p_loc[2] + 1])
    CairoMakie.Colorbar(
        place[p_loc_cb...],
        plot,
        label = colorbar_label;
        cb_kwargs...,
    )
end

"""
Private function to define `sliced_` functions.

It slices a given variable and applies `func` to it.
"""
function _sliced_plot_generic(
    func,
    fig,
    var,
    cut;
    p_loc,
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
)
    isnothing(cut) && (cut = Dict())

    var_sliced = var

    for (dim_name, val) in cut
        var_sliced = ClimaAnalysis._slice_general(var_sliced, val, dim_name)
    end

    func(fig, var_sliced; p_loc, more_kwargs)
end

"""
Private function to define `plot!` functions.

It composes a `cut` from given `kwargs`. Used with `sliced` functions.
"""
function _plot_generic_kwargs(
    func,
    fig,
    var;
    p_loc,
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
    kwargs...,
)
    cut = Dict("$k" => v for (k, v) in kwargs)
    length(cut) == 0 && (cut = nothing)
    return func(fig, var, cut; p_loc, more_kwargs)
end

"""
    sliced_heatmap!(fig::CairoMakie.Figure,
                    var::ClimaAnalysis.OutputVar,
                    cut::Union{Nothing, AbstractDict{String, <: Real}};
                    p_loc = (1,1),
                    more_kwargs,
                    )
    sliced_heatmap!(grid_layout::CairoMakie.GridLayout,
                    var::ClimaAnalysis.OutputVar,
                    cut::Union{Nothing, AbstractDict{String, <: Real}};
                    p_loc = (1,1),
                    more_kwargs,
                    )

Take a `var`iable, slice as directed, and plot a 2D heatmap in the given place and
location.

The place can be a `Figure` or a `GridLayout`.

The plot comes with labels, units, and a colorbar.

Arguments
=========

If the variable is not 2D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` heatmap at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.sliced_heatmap!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
)
    return _sliced_plot_generic(
        Visualize.heatmap2D!,
        place,
        var,
        cut;
        p_loc,
        more_kwargs,
    )
end

"""
    heatmap!(place::MakiePlace,
             var::ClimaAnalysis.OutputVar;
             p_loc = (1,1),
             more_kwargs,
             kwargs...
            )

Syntactic sugar for `sliced_heatmap` with `kwargs` instead of `cut`.

Example
=======

`heatmap!(fig, var, time = 100, lat = 70)` plots a heatmap by slicing `var` along
the time nearest to 100 and latitude nearest 70.

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.heatmap!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_heatmap!,
        place,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

"""
    line_plot1D!(place::CairoMakie.Figure,
                 var::ClimaAnalysis.OutputVar;
                 p_loc = (1,1),
                 more_kwargs
                 )
    line_plot1D!(place::CairoMakie.GridLayout,
                 var::ClimaAnalysis.OutputVar;
                 p_loc = (1,1),
                 more_kwargs
                 )

Plot a line plot of the given 1D `var`iable in the given place and location.
The place can be a `Figure` or a `GridLayout`.

The plot comes with labels, units.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.

A special argument that can be passed to `:axis` is `:dim_on_y`, which puts the dimension on
the y axis instead of the variable. This is useful to plot columns with `z` on the vertical
axis instead of the horizontal one.

"""
function Visualize.line_plot1D!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :axis => Dict()),
)
    length(var.dims) == 1 || error("Can only plot 1D variables")

    dim_name = var.index2dim[]
    dim = var.dims[dim_name]

    units = var.attributes["units"]
    short_name = var.attributes["short_name"]

    dim_units = var.dim_attributes[dim_name]["units"]

    axis_kwargs = get(more_kwargs, :axis, Dict())
    plot_kwargs = get(more_kwargs, :plot, Dict())

    title = get(axis_kwargs, :title, var.attributes["long_name"])
    xlabel = get(axis_kwargs, :xlabel, "$dim_name [$dim_units]")
    ylabel = get(axis_kwargs, :ylabel, "$short_name [$units]")

    x, y = dim, var.data

    if get(axis_kwargs, :dim_on_y, false)
        xlabel, ylabel = ylabel, xlabel
        x, y = y, x
        # dim_on_y is not a real keyword for Axis, so we have to remove it from the
        # arguments. Since axis_kwargs is a Pairs object, we have to go through its
        # underlying dictionary first
        axis_kwargs_dict = Dict(axis_kwargs)
        pop!(axis_kwargs_dict, :dim_on_y)
        axis_kwargs = pairs(axis_kwargs_dict)
    end

    CairoMakie.Axis(place[p_loc...]; title, xlabel, ylabel, axis_kwargs...)
    CairoMakie.lines!(x, y; plot_kwargs...)
end

"""
    sliced_line_plot!(place::CairoMakie.Figure,
                      var::ClimaAnalysis.OutputVar,
                      cut::Union{Nothing, AbstractDict{String, <: Real}};
                      p_loc = (1,1),
                      more_kwargs
                      )
    sliced_line_plot!(place::CairoMakie.GridLayout,
                      var::ClimaAnalysis.OutputVar,
                      cut::Union{Nothing, AbstractDict{String, <: Real}};
                      p_loc = (1,1),
                      more_kwargs
                      )

Take a `var`iable, slice as directed, and plot a 1D line plot in the given place and
location. The place can be a `Figure` or a `GridLayout`.

The plot comes with labels, and units.

Arguments
=========

If the variable is not 1D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` heatmap at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)


Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.sliced_line_plot!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :axis => Dict()),
)
    return _sliced_plot_generic(
        Visualize.line_plot1D!,
        place,
        var,
        cut;
        p_loc,
        more_kwargs,
    )
end

"""
    line_plot!(place::CairoMakie.Figure,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs,
               kwargs...
               )
    line_plot!(place::CairoMakie.GridLayout,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs,
               kwargs...
               )

Syntactic sugar for `sliced_line_plot` with `kwargs` instead of `cut`.

Example
=======

`line_plot!(fig, var, time = 100, lat = 70)` plots a line plot by slicing `var` along
the time nearest to 100 and latitude nearest 70.

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.line_plot!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_line_plot!,
        place,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

"""
    sliced_plot!(place::CairoMakie.Figure,
                 var::ClimaAnalysis.OutputVar,
                 cut::Union{Nothing, AbstractDict{String, <: Real}};
                 p_loc = (1,1),
                 more_kwargs
                 )
    sliced_plot!(place::CairoMakie.GridLayout,
                 var::ClimaAnalysis.OutputVar,
                 cut::Union{Nothing, AbstractDict{String, <: Real}};
                 p_loc = (1,1),
                 more_kwargs
                 )

Take a `var`iable, slice as directed, and plot a 1D line plot or 2D heatmap in the given place and
location. The place can be a `Figure` or a `GridLayout`.

The plot comes with labels, and units (and possibly a colorbar).

Arguments
=========

If the variable is not 1D/2D, `cut` has to be a dictionary that maps the dimension that has to
be sliced and the value where to cut.

For example, if `var` has four dimensions: `time`, `long`, `lat`, `z`, this function can be
used to plot a `lat-long` heatmap at fixed `time` and `z`. Assuming we want to plot
time `100.` and altitude `50.`, `cut` should be `Dict("time" => 100., "z" => 50.)`.

This function assumes that the following attributes are available:
- long_name
- short_name
- units (also for the dimensions)

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.sliced_plot!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
)
    initial_dim = length(var.dims)
    removed_dims = isnothing(cut) ? 0 : length(cut)
    final_dim = initial_dim - removed_dims

    if final_dim == 1
        fun = Visualize.line_plot1D!
    elseif final_dim == 2
        fun = Visualize.heatmap2D!
    else
        error("Sliced variable has $final_dim dimensions (needed 1 or 2)")
    end

    return _sliced_plot_generic(fun, place, var, cut; p_loc, more_kwargs)
end


"""
    plot!(place::CairoMakie.Figure,
          var::ClimaAnalysis.OutputVar;
          p_loc = (1,1),
          more_kwargs,
          kwargs...
          )
    plot!(place::CairoMakie.GridLayout,
          var::ClimaAnalysis.OutputVar;
          p_loc = (1,1),
          more_kwargs,
          kwargs...
          )

Syntactic sugar for `sliced_plot` with `kwargs` instead of `cut`.

Example
=======

`line_plot!(fig, var, time = 100, lat = 70)` plots a line plot or a heatmap by slicing
`var` along the time nearest to 100 and latitude nearest 70.

Additional arguments to the plotting and axis functions
=======================================================

`more_kwargs` can be a dictionary that maps symbols to additional options for:
- the axis (`:axis`)
- the plotting function (`:plot`)
- the colorbar (`:cb`)

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.
"""
function Visualize.plot!(
    place::MakiePlace,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_plot!,
        place,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

end
