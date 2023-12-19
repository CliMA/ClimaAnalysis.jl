module CairoMakieExt

import CairoMakie
import ClimaAnalysis
import ClimaAnalysis: Visualize

"""
    heatmap2D!(fig::CairoMakie.Figure,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs)

Plot a heatmap of the given 2D `var`iable in the given `fig`ure and location.

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
    fig::CairoMakie.Figure,
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

    title = var.attributes["long_name"]
    xlabel = "$dim1_name [$dim1_units]"
    ylabel = "$dim2_name [$dim2_units]"

    # TODO: (FIX)
    #
    # Override title/xlabel/ylabel

    CairoMakie.Axis(
        fig[p_loc...];
        title,
        xlabel,
        ylabel,
        get(more_kwargs, :axis, Dict())...,
    )

    plot = CairoMakie.heatmap!(
        dim1,
        dim2,
        var.data;
        get(more_kwargs, :plot, Dict())...,
    )

    p_loc_cb = Tuple([p_loc[1], p_loc[2] + 1])
    CairoMakie.Colorbar(
        fig[p_loc_cb...],
        plot,
        label = colorbar_label;
        get(more_kwargs, :cb, Dict())...,
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
        var_sliced = ClimaAnalysis.slice_general(var_sliced, val, dim_name)
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
    sliced_heatmap!(
                        fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar,
                        cut::Union{Nothing, AbstractDict{String, <: Real}};
                        p_loc = (1,1),
                        more_kwargs,
                        )

Take a `var`iable, slice as directed, and plot a 2D heatmap in the given `fig`ure and
location.

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
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
)
    return _sliced_plot_generic(
        Visualize.heatmap2D!,
        fig,
        var,
        cut;
        p_loc,
        more_kwargs,
    )
end

"""
    heatmap!(
                fig::CairoMakie.Figure,
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
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_heatmap!,
        fig,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

"""
    line_plot1D!(
                 fig::CairoMakie.Figure,
                 var::ClimaAnalysis.OutputVar;
                 p_loc = (1,1),
                 more_kwargs
                 )

Plot a heatmap of the given 2D `var`iable in the given `fig`ure and location.

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

The values are splatted in the relevant functions. Populate them with a
Dictionary of `Symbol`s => values to pass additional options.

"""
function Visualize.line_plot1D!(
    fig::CairoMakie.Figure,
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

    title = var.attributes["long_name"]
    xlabel = "$dim_name [$dim_units]"
    ylabel = "$short_name [$units]"

    CairoMakie.Axis(
        fig[p_loc...];
        title,
        xlabel,
        ylabel,
        get(more_kwargs, :axis, Dict())...,
    )
    CairoMakie.lines!(dim, var.data; title, get(more_kwargs, :plot, Dict())...)
end

"""
    sliced_line_plot!(
                        fig::CairoMakie.Figure,
                        var::ClimaAnalysis.OutputVar,
                        cut::Union{Nothing, AbstractDict{String, <: Real}};
                        p_loc = (1,1),
                        more_kwargs
                        )

Take a `var`iable, slice as directed, and plot a 1D line plot in the given `fig`ure and
location.

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
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar,
    cut::Union{Nothing, AbstractDict{String, <:Real}} = nothing;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :axis => Dict()),
)
    return _sliced_plot_generic(
        Visualize.line_plot1D!,
        fig,
        var,
        cut;
        p_loc,
        more_kwargs,
    )
end

"""
    line_plot!(
                fig::CairoMakie.Figure,
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
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_line_plot!,
        fig,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

"""
    sliced_plot!(
                 fig::CairoMakie.Figure,
                 var::ClimaAnalysis.OutputVar,
                 cut::Union{Nothing, AbstractDict{String, <: Real}};
                 p_loc = (1,1),
                 more_kwargs
                )

Take a `var`iable, slice as directed, and plot a 1D line plot or 2D heatmap in the given `fig`ure and
location.

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
    fig::CairoMakie.Figure,
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

    return _sliced_plot_generic(fun, fig, var, cut; p_loc, more_kwargs)
end


"""
    plot!(
          fig::CairoMakie.Figure,
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
    fig::CairoMakie.Figure,
    var::ClimaAnalysis.OutputVar;
    p_loc = (1, 1),
    more_kwargs = Dict(:plot => Dict(), :cb => Dict(), :axis => Dict()),
    kwargs...,
)
    _plot_generic_kwargs(
        Visualize.sliced_plot!,
        fig,
        var;
        p_loc,
        more_kwargs,
        kwargs...,
    )
end

end
