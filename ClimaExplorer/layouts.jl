"""
    layout(path)

Construct most of the webapp (and do error handling).
"""
function layout(path)
    if isdir(path)
        simdir = ClimaAnalysis.SimDir(path)
        if isempty(simdir)
            return Card(
                DOM.h4("⚠️ $path does not contain ClimaDiagnostics data"),
            )
        else
            return layout_valid(simdir)
        end
    else
        return layout_invalid(path)
    end
end

function layout_invalid(path)
    return Card(DOM.h4("⚠️ $path is not a valid path!"))
end

function _prepare_menus(simdir::ClimaAnalysis.SimDir)
    variables = ClimaAnalysis.available_vars(simdir) |> collect |> sort
    vars_menu = Dropdown(variables)
    reductions =
        ClimaAnalysis.available_reductions(
            simdir;
            short_name = vars_menu.value.val,
        ) |> collect
    reductions_menu = Dropdown(reductions)
    periods =
        ClimaAnalysis.available_periods(
            simdir;
            short_name = vars_menu.value.val,
            reduction = reductions_menu.value.val,
        ) |> collect
    periods_menu = Dropdown(periods)
    return vars_menu, reductions_menu, periods_menu
end

"""
    layout_layout(simdir::ClimaAnalysis.SimDir)

Construct most of the webapp, knowing that the `simdir` is valid.
"""
function layout_valid(simdir::ClimaAnalysis.SimDir)
    # Prepare the menus.
    vars_menu, reductions_menu, periods_menu = _prepare_menus(simdir)

    on(vars_menu.value) do short_name
        reductions =
            ClimaAnalysis.available_reductions(simdir; short_name) |> collect
        # setindex! triggers on(reductions_menu.value)
        setindex!(reductions_menu.options, reductions)
    end

    on(reductions_menu.value) do reduction
        periods =
            ClimaAnalysis.available_periods(
                simdir;
                short_name = vars_menu.value.val,
                reduction,
            ) |> collect
        # setindex! triggers on(periods_menu.value)
        setindex!(periods_menu.options, periods)
    end

    # We have to do this little trick to support Observable vars because
    # ClimaAnalysis.get() returns OutputVars with different types, so we
    # cannot simply map(observable) because it does not work when types are
    # not uniform. We just have to unpack the value when using the var
    generic_ref = Any[get(
        simdir;
        short_name = vars_menu.value.val,
        reduction = reductions_menu.value.val,
        period = periods_menu.value.val,
    )]
    var = Observable(generic_ref)

    on(periods_menu.value) do period
        short_name = vars_menu.value.val
        reduction = reductions_menu.value.val
        setindex!(var, [get(simdir; short_name, reduction = reduction, period)])
    end

    header_dom = Card(
        Grid(
            Centered(DOM.div("Short name:", vars_menu)),
            Centered(DOM.div("Reduction:", reductions_menu)),
            Centered(DOM.div("Period:", periods_menu)),
            columns = "1fr 1fr 1fr",
        ),
    )

    # This is the main body. It contains the side menu and the figure(s).
    content_dom = map(var) do var_complete
        content_layout(var_complete[])
    end
    return DOM.div(header_dom, content_dom)
end

"""
    content_layout(simdir::ClimaAnalysis.OutputVar)

Prepare the body of the webapp.

This function dispatches the layout so that boxes/sphere/columns/surfaces and so
on can have their customize layout.
"""
function content_layout(var_complete::ClimaAnalysis.OutputVar)
    # TODO: Add other content_layouts (ie, for box and column)
    if ClimaAnalysis.has_altitude(var_complete)
        return sphere_content_layout(var_complete)
    else
        return spherical_shell_content_layout(var_complete)
    end
end


"""
    play_button(time_slider, var_complete::ClimaAnalysis.OutputVar)

Add a button that steps the `time_slider`

TODO: Add a pause button
"""
function play_button(time_slider, var_complete::ClimaAnalysis.OutputVar)
    play_button = Button("▶")
    on(play_button.value) do _
        play_times = filter(
            t -> t >= time_slider.value.val,
            ClimaAnalysis.times(var_complete),
        )
        for t in play_times
            setindex!(time_slider, t)
        end
    end
    return play_button
end

"""
    time_units(var_complete::ClimaAnalysis.OutputVar)

Return the units of time in `var_complete`.
"""
function time_units(var_complete::ClimaAnalysis.OutputVar)
    return var_complete.dim_attributes[ClimaAnalysis.time_name(var_complete)]["units"]
end

"""
    altitude_units(var_complete::ClimaAnalysis.OutputVar)

Return the units of altitude in `var_complete`.
"""
function altitude_units(var_complete::ClimaAnalysis.OutputVar)
    return var_complete.dim_attributes[ClimaAnalysis.altitude_name(var_complete)]["units"]
end

"""
    sphere_content_layout(var::ClimaAnalysis.OutputVar)

Set up the page content_layout and content for the given var, assuming it is a 3D sphere.

A sphere content_layout has two columns:

Left column
- A slider for the time
- A play button
- A slider for the level

Right column:
- Figure with the 2D map at given time and level
"""
function sphere_content_layout(var_complete::ClimaAnalysis.OutputVar)
    time_slider = Slider(ClimaAnalysis.times(var_complete))
    level_slider = Slider(ClimaAnalysis.altitudes(var_complete))

    on(time_slider) do _
        # Trigger level_slider
        setindex!(level_slider, level_slider.value.val)
    end

    alt_name = ClimaAnalysis.altitude_name(var_complete)

    var = map(level_slider) do level
        # We need this because variables might have different altitude dimensions (e.g., z
        # vs z_reference)
        kwargs = (Symbol(alt_name) => level,)
        ClimaAnalysis.slice(
            var_complete;
            time = time_slider.value.val,
            kwargs...,
        )
    end

    fig = map(var) do sliced_var
        fig = WGLMakie.Figure(size = (1200, 760), fontsize = 30)
        ClimaAnalysis.Visualize.contour2D_on_globe!(fig, sliced_var)
        fig
    end

    controls = Centered(
        Grid(
            DOM.div("Time:  ", time_slider.value, " ", time_units(var_complete)),
            time_slider,
            play_button(time_slider, var_complete),
            DOM.div("Altitude:  ", level_slider.value, " ",  altitude_units(var_complete)),
            level_slider;
        ),
    )

    return Grid(
        Card(controls),
        Card(Centered(DOM.div(fig)));
        columns = "1fr 5fr",
    )
end

using WGLMakie
"""
    spherical_shell_content_layout(var::ClimaAnalysis.OutputVar)

Set up the page content_layout and content for the given var, assuming it is a (2D) spherical shell.

A spherical sphere content_layout has two columns:

Left column
- A slider for the time
- A play button

Right column:
- Figure with the 2D map at given time and level
"""
function spherical_shell_content_layout(var_complete::ClimaAnalysis.OutputVar)
    time_slider = Slider(ClimaAnalysis.times(var_complete))
    fig = WGLMakie.Figure(size = (1200, 760), fontsize = 30)
    ax = GeoMakie.GeoAxis(fig[1, 1])

    #time = time_slider.value # time is an Observable
    #var = @lift(ClimaAnalysis.slice(var_complete; $time)) # var is an Observable

    var = map(time_slider) do time
        ClimaAnalysis.slice(var_complete; time)
    end

	ClimaAnalysis.Visualize.contour2D_on_globe!(fig, var, ax) 
 
    controls = Centered(
        Grid(
            DOM.div("Time:  ", time_slider.value, " $(time_units(var_complete))"),
            time_slider,
            play_button(time_slider, var_complete),
        ),
    )

    return Grid(
        Card(controls),
        Card(Centered(DOM.div(fig)));
        columns = "1fr 5fr",
    )
end
