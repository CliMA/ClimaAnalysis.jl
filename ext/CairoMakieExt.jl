module CairoMakieExt

import CairoMakie
import ClimaAnalysis: Visualize

"""
    contour_plot!(fig; X, Y, Z)

Create and save a Makie plot object, which consists of either a line-plot (X,Y)
or a surface plot ((X,Y),Z).
"""
function Visualize.contour_plot!(
    fig::CairoMakie.Figure;
    X = nothing,
    Y = nothing,
    Z = nothing,
    p_loc::Tuple = (1, 1),
    title = "",
    xlabel = "",
    ylabel = "",
    colorbar_label = "",
    xscale = identity,
    yscale = identity,
    linewidth = 6,
)
    if isnothing(Z)
        generic_axis = fig[p_loc[1], p_loc[2]] = CairoMakie.GridLayout()
        CairoMakie.Axis(
            generic_axis[1, 1];
            title,
            xlabel,
            ylabel,
            xscale,
            yscale,
        )
        CairoMakie.CairoMakie.lines!(
            X,
            Y;
            title,
            linewidth,
            label = colorbar_label,
        )
    else
        generic_axis = fig[p_loc[1], p_loc[2]] = CairoMakie.GridLayout()
        CairoMakie.Axis(generic_axis[1, 1]; title, xlabel, ylabel, yscale)

        # custom_levels is a workaround for plotting constant fields with CairoMakie
        custom_levels =
            minimum(Z) ≈ maximum(Z) ? (minimum(Z):0.1:(minimum(Z) + 0.2)) : 25

        generic_plot = CairoMakie.contourf!(X, Y, Z, levels = custom_levels)
        CairoMakie.Colorbar(
            generic_axis[1, 2],
            generic_plot,
            label = colorbar_label,
        )
    end
end

end
