module Visualize

function heatmap2D! end

function sliced_heatmap! end

function heatmap! end

function line_plot1D! end

function sliced_line_plot! end

function line_plot! end

function sliced_plot! end

function plot! end

function heatmap2D_on_globe! end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f in [heatmap2D!, sliced_heatmap!, heatmap!, line_plot1D!, sliced_line_plot!, line_plot!, sliced_plot!, plot!, heatmap2D_on_globe!]
            if isempty(methods(exc.f))
                print(io, "\n$(exc.f) has no methods, yet. Makie has to be loaded for the plotting extension to be activated. Run `using Makie`, `using CairoMakie`, `using GLMakie` or any other package that also loads Makie.")
            end
        end
    end
end

end
