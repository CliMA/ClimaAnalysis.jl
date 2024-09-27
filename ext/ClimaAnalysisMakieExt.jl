module ClimaAnalysisMakieExt

import Makie
import ClimaAnalysis
import ClimaAnalysis: Visualize

MakiePlace = Union{Makie.Figure, Makie.GridLayout}

"""
    heatmap2D!(fig::Makie.Figure,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs)
    heatmap2D!(grid_layout::Makie.GridLayout,
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

    units = ClimaAnalysis.units(var)
    short_name = var.attributes["short_name"]
    colorbar_label = "$short_name [$units]"
    dim1_units = var.dim_attributes[dim1_name]["units"]
    dim2_units = var.dim_attributes[dim2_name]["units"]

    axis_kwargs = get(more_kwargs, :axis, Dict())
    plot_kwargs = get(more_kwargs, :plot, Dict())
    cb_kwargs = get(more_kwargs, :cb, Dict())

    var.attributes["long_name"] =
        ClimaAnalysis.Utils.warp_string(var.attributes["long_name"])

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

    Makie.Axis(place[p_loc...]; title, xlabel, ylabel, axis_kwargs...)

    plot = Makie.heatmap!(dim1, dim2, var.data; plot_kwargs...)

    p_loc_cb = Tuple([p_loc[1], p_loc[2] + 1])
    Makie.Colorbar(
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
        var_sliced = ClimaAnalysis.Var._slice_general(var_sliced, val, dim_name)
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
    sliced_heatmap!(fig::Makie.Figure,
                    var::ClimaAnalysis.OutputVar,
                    cut::Union{Nothing, AbstractDict{String, <: Real}};
                    p_loc = (1,1),
                    more_kwargs,
                    )
    sliced_heatmap!(grid_layout::Makie.GridLayout,
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
    line_plot1D!(place::Makie.Figure,
                 var::ClimaAnalysis.OutputVar;
                 p_loc = (1,1),
                 more_kwargs
                 )
    line_plot1D!(place::Makie.GridLayout,
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

    units = ClimaAnalysis.units(var)
    short_name = var.attributes["short_name"]

    dim_units = var.dim_attributes[dim_name]["units"]

    axis_kwargs = get(more_kwargs, :axis, Dict())
    plot_kwargs = get(more_kwargs, :plot, Dict())

    var.attributes["long_name"] =
        ClimaAnalysis.Utils.warp_string(var.attributes["long_name"])

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

    Makie.Axis(place[p_loc...]; title, xlabel, ylabel, axis_kwargs...)
    Makie.lines!(x, y; plot_kwargs...)
end

"""
    sliced_line_plot!(place::Makie.Figure,
                      var::ClimaAnalysis.OutputVar,
                      cut::Union{Nothing, AbstractDict{String, <: Real}};
                      p_loc = (1,1),
                      more_kwargs
                      )
    sliced_line_plot!(place::Makie.GridLayout,
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
    line_plot!(place::Makie.Figure,
               var::ClimaAnalysis.OutputVar;
               p_loc = (1,1),
               more_kwargs,
               kwargs...
               )
    line_plot!(place::Makie.GridLayout,
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
    sliced_plot!(place::Makie.Figure,
                 var::ClimaAnalysis.OutputVar,
                 cut::Union{Nothing, AbstractDict{String, <: Real}};
                 p_loc = (1,1),
                 more_kwargs
                 )
    sliced_plot!(place::Makie.GridLayout,
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
    plot!(place::Makie.Figure,
          var::ClimaAnalysis.OutputVar;
          p_loc = (1,1),
          more_kwargs,
          kwargs...
          )
    plot!(place::Makie.GridLayout,
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

"""
    _to_unitrange(x::Number, lo::Number, hi::Number)

Linearly transform x ∈ [lo, hi] to [0, 1].
"""
_to_unitrange(x::Number, lo::Number, hi::Number) = (x - lo) / (hi - lo)

"""
    _constrained_cmap(cols::Vector, lo, hi; mid = 0, categorical = false, rev = false)
    _constrained_cmap(cols::Makie.ColorScheme, lo, hi; mid = 0, categorical = false, rev = false)

Constrain a colormap to a given range.

Given a colormap implicitly defined in `± maximum(abs, (lo, hi))`, constrain it to the range
[lo, hi]. This is useful to ensure that a colormap which is desired to diverge
symmetrically around zero maps the same color intensity to the same magnitude.

# Arguments
- `cols`: a vector of colors, or a ColorScheme
- `lo`: lower bound of the range
- `hi`: upper bound of the range

# Keyword Arguments
- `mid`: midpoint of the range  # TODO: test `mid` better
- `categorical`: flag for whether returned colormap should be categorical or continuous
- `rev`: flag for whether to reverse the colormap before constraining cmap

# Returns
- `cmap::Makie.ColorGradient`: a colormap
"""
function Visualize._constrained_cmap(
    cols::Vector,
    lo,
    hi;
    mid = 0,
    categorical = false,
    rev = false,
)
    Visualize._constrained_cmap(
        Makie.ColorScheme(cols),
        lo,
        hi;
        mid,
        categorical,
        rev,
    )
end

function Visualize._constrained_cmap(
    cols::Makie.ColorScheme,
    lo,
    hi;
    mid = 0,
    categorical = false,
    rev = false,
)
    # Reverse colorscheme if requested, don't reverse below in `cgrad`
    rev && (cols = reverse(cols))
    absmax = maximum(abs, (lo, hi) .- mid)

    # Map lo, hi ∈ [-absmax, absmax] onto [0,1] to sample their corresponding colors
    lo_m, hi_m = _to_unitrange.((lo, hi) .- mid, -absmax, absmax)

    # Values on [0,1] where each color in cols is defined
    colsvals = range(0, 1; length = length(cols))

    # Filter colsvals, keep only values in [lo_m, hi_m] + the endpoints lo_m and hi_m
    filter_colsvals =
        filter(x -> lo_m <= x <= hi_m, unique([lo_m; colsvals; hi_m]))

    # Select colors in filtered range; interpolate new low and hi colors
    newcols = Makie.get(cols, filter_colsvals)

    # Values on [0,1] where the new colors are defined
    new_colsvals = _to_unitrange.(filter_colsvals, lo_m, hi_m)
    cmap = Makie.cgrad(newcols, new_colsvals; categorical, rev = false)
    return cmap
end

"""
    Visualize.plot_boxplot!(fig,
                            rmse_var::ClimaAnalysis.RMSEVariable;
                            model_names = ["CliMA"],
                            ploc = (1, 1),
                            best_and_worst_category_name = "ANN",
                            legend_text_width = 10)

Plot a Tukey style boxplot for each category in `rmse_var`.

The best and worst single models are found for the category `best_and_worst_category_name`
and are plotted on the boxplot. When finding the best and worst single models, any models in
`model_names` will be excluded. Additionally, any model in `model_names` will also be
plotted on the boxplot.

The parameter `ploc` determines where to place the plot on the figure. The parameter
`legend_text_width` determines the number of characters on each line in the legend.
"""
function Visualize.plot_boxplot!(
    fig,
    rmse_var::ClimaAnalysis.RMSEVariable;
    model_names = ["CliMA"],
    ploc = (1, 1),
    best_and_worst_category_name = "ANN",
    legend_text_width = 10,
)
    # Unit checking
    ClimaAnalysis.Leaderboard._unit_check(rmse_var)

    num_cats = length(rmse_var.category2index)
    units = values(rmse_var.units) |> collect |> first

    # Title and labels for x-axis and y-axis
    ax = Makie.Axis(
        fig[ploc...],
        ylabel = "$(rmse_var.short_name) [$units]",
        xticks = (1:num_cats, ClimaAnalysis.category_names(rmse_var)),
        title = "Global RMSE $(rmse_var.short_name) [$units]",
    )

    # Set up for box plot
    cats = reduce(
        vcat,
        [
            fill(cat_val, length(rmse_var.model2index)) for
            cat_val in 1:length(rmse_var.category2index)
        ],
    )
    vals = reduce(vcat, rmse_var.RMSEs)

    # Filter out NaNs because we can't plot with NaNs
    not_nan_idices = findall(!isnan, vals)
    cats = cats[not_nan_idices]
    vals = vals[not_nan_idices]

    # Add box plot
    Makie.boxplot!(
        ax,
        cats,
        vals,
        whiskerwidth = 1,
        width = 0.35,
        mediancolor = :black,
        color = :gray,
        whiskerlinewidth = 1,
    )

    # Delete any model in model_names to exclude them when finding best and worst models
    rmse_var_delete =
        ClimaAnalysis.Leaderboard._delete_model(rmse_var, model_names...)

    # Find and plot best and worst models
    absolute_worst_values, absolute_worst_model_name =
        ClimaAnalysis.find_worst_single_model(
            rmse_var_delete,
            category_name = best_and_worst_category_name,
        )
    absolute_best_values, absolute_best_model_name =
        ClimaAnalysis.find_best_single_model(
            rmse_var_delete,
            category_name = best_and_worst_category_name,
        )
    best_pt = Makie.scatter!(
        ax,
        1:num_cats,
        absolute_best_values,
        label = absolute_best_model_name,
    )
    worst_pt = Makie.scatter!(
        ax,
        1:num_cats,
        absolute_worst_values,
        label = absolute_worst_model_name,
    )

    # Plotting the median model
    median_pt = Makie.scatter!(
        ax,
        1:num_cats,
        ClimaAnalysis.median(rmse_var),
        label = "Median",
        color = :black,
        marker = :hline,
        markersize = 10,
        visible = false,
    )

    # Keep track of points plotted by scatter! and names for plotting the legend
    # later
    pts_on_boxplot = [median_pt, best_pt, worst_pt]
    names_on_legend = vcat(
        ["Median", absolute_best_model_name, absolute_worst_model_name],
        model_names,
    )
    # Plot CliMA model and other models
    for model_name in model_names
        ClimaAnalysis.Leaderboard._model_name_check(rmse_var, model_name)
        if model_name == "CliMA"
            Makie.scatter!(
                ax,
                1:num_cats,
                rmse_var[model_name],
                label = model_name,
                marker = :star5,
                markersize = 20,
                color = :green,
            ) |> pt -> push!(pts_on_boxplot, pt)
        else
            Makie.scatter!(
                ax,
                1:num_cats,
                rmse_var[model_name],
                label = model_name,
                markersize = 20,
                color = :red,
            ) |> pt -> push!(pts_on_boxplot, pt)
        end
    end

    # Add a new line character every `legend_text_width` characters to prevent legend from
    # overflowing onto the box plot
    # Soln from
    # https://stackoverflow.com/questions/40545980/insert-a-newline-character-every-10-characters-in-a-string-using-julia
    names_on_legend = map(
        name ->
            replace(name, Regex("(.{$legend_text_width})") => s"\1\n") |>
            (name -> rstrip(name, '\n')),
        names_on_legend,
    )
    # Hack to make legend appear better
    Makie.axislegend(ax, pts_on_boxplot, names_on_legend)
    Makie.scatter!(ax, [num_cats + 2.5], [0.1], markersize = 0.01)
end

"""
    Visualize.plot_leaderboard!(fig,
                                rmse_vars::ClimaAnalysis.RMSEVariable...;
                                ploc = (1, 1),
                                model_names = ["CliMA"],
                                best_category_name = "ANN")

Plot a heatmap over the categories and models. The models that appear is the best model
as found for the category `best_category_name` and any other models in `model_names`. The
root mean squared errors for each variable of interest is normalized by dividing over the
median root mean squared error of each variable.

The parameter `ploc` determines where to place the plot on the figure.
"""
function Visualize.plot_leaderboard!(
    fig,
    rmse_vars::ClimaAnalysis.RMSEVariable...;
    ploc = (1, 1),
    model_names = ["CliMA"],
    best_category_name = "ANN",
)
    # Check if rmse_model_vars all have the same categories
    categories_names = ClimaAnalysis.category_names.(rmse_vars)
    categories_same = length(unique(categories_names)) == 1
    categories_same ||
        error("Categories are not all the same across the RMSEVariable")

    rmse_var = first(rmse_vars)
    categ_names = ClimaAnalysis.category_names(rmse_var)
    num_variables = length(rmse_vars)
    num_boxes = length(categ_names) # number of categories
    num_models = 1 + length(model_names)  # best model plus the other models in model_names

    # Initialize variables we need for storing RMSEs for plotting and short names for axis
    rmse_normalized_arr = zeros(num_boxes * num_models, num_variables)
    short_names = String[]

    for (idx, var) in enumerate(reverse(rmse_vars))
        # Get all the short name of the rmse_vars
        push!(short_names, var.short_name)

        # Compute median and best values for RMSE
        med_vals = ClimaAnalysis.median(var)
        best_vals, _ = ClimaAnalysis.find_best_single_model(
            var,
            category_name = best_category_name,
        )

        # Find normalized values for the models we are interested in and the normalized best
        # value and store them
        normalized_vals = [var[model] ./ med_vals for model in model_names]
        normalized_vals = reduce(vcat, normalized_vals)
        rmse_normalized_arr[:, idx] =
            vcat(normalized_vals, best_vals ./ med_vals)'
    end

    # Finding the midpoint for placing labels
    start_x_tick = div(num_boxes, 2, RoundUp)

    ax_bottom_and_left = Makie.Axis(
        fig[ploc...],
        yticks = (1:length(short_names), short_names),
        xticks = (
            [start_x_tick, start_x_tick + num_boxes],
            vcat(model_names, ["Best model"]),
        ),
        aspect = num_boxes * num_models,
        xgridvisible = false,
        ygridvisible = false,
    )
    ax_top = Makie.Axis(
        fig[ploc...],
        xaxisposition = :top,
        xticks = (0.5:1.0:length(categ_names), categ_names),
        aspect = num_boxes * num_models,
        xgridvisible = false,
        ygridvisible = false,
    )
    Makie.hidespines!(ax_top)
    Makie.hideydecorations!(ax_top)

    colormap = Makie.Reverse(:RdYlGn)

    # Filter out NaNs here because we need to take the maximum and extrema for the
    # colorrange and limits
    rmse_no_nan_vec = rmse_normalized_arr |> vec |> filter(!isnan)
    Makie.heatmap!(
        ax_bottom_and_left,
        rmse_normalized_arr,
        colormap = colormap,
        # Trick to exclude the zeros
        lowclip = :white,
        colorrange = (1e-10, maximum(rmse_no_nan_vec)),
    )
    for idx in eachindex(model_names)
        Makie.vlines!(ax_top, num_boxes * idx, color = :black, linewidth = 3.0)
    end
    row, col = ploc
    col += 1
    Makie.Colorbar(
        fig[row, col],
        limits = extrema(rmse_no_nan_vec),
        label = "RMSE/median(RMSE)",
        colormap = colormap,
    )
end

end
