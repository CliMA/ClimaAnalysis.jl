module Visualize

export plot!

function _constrained_cmap end

function oceanmask end

function landmask end

function heatmap2D! end

function sliced_heatmap! end

function heatmap! end

function line_plot1D! end

function sliced_line_plot! end

function line_plot! end

function sliced_plot! end

function plot! end

function contour2D_on_globe! end

function heatmap2D_on_globe! end

function plot_bias_on_globe! end

function plot_boxplot! end

function plot_leaderboard! end

extension_fns = [
    :Makie => [
        :heatmap2D!,
        :sliced_heatmap!,
        :heatmap!,
        :line_plot1D!,
        :sliced_line_plot!,
        :line_plot!,
        :sliced_plot!,
        :plot!,
        :plot_boxplot!,
        :plot_leaderboard!,
        :_constrained_cmap,
    ],
    :GeoMakie => [
        :oceanmask,
        :landmask,
        :heatmap2D_on_globe!,
        :contour2D_on_globe!,
        :plot_bias_on_globe!,
    ],
]

"""
    is_pkg_loaded(pkg::Symbol)

Check if `pkg` is loaded or not.
"""
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

function __init__()
    # Register error hint if a package is not loaded
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(
            MethodError,
        ) do io, exc, _argtypes, _kwargs
            for (pkg, fns) in extension_fns
                if Symbol(exc.f) in fns && !is_pkg_loaded(pkg)
                    if pkg == :Makie
                        print(
                            io,
                            "\nImport one of the Makie backends (GLMakie, CairoMakie, WGLMakie, RPRMakie, etc.) to enable `$(exc.f)`.";
                        )
                    elseif pkg == :GeoMakie
                        print(io, "\nImport GeoMakie to enable `$(exc.f)`.";)
                    end
                end
            end
        end
    end
end

end
