module ClimaAnalysis
import Reexport: @reexport

include("Utils.jl")
import .Utils

include("Numerics.jl")

include("Var.jl")
@reexport using .Var
include("Sim.jl")
@reexport using .Sim
include("Leaderboard.jl")
@reexport using .Leaderboard

include("Visualize.jl")
@reexport using .Visualize

include("Atmos.jl")

include("Template.jl")

# In case we want to turn Unitful into an extension
include("../ext/ClimaAnalysisUnitfulExt.jl")

end # module ClimaAnalysis
