using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "Docstrings" begin @time include("doctest.jl") end
@safetestset "Format" begin @time include("format.jl") end

@safetestset "Utils" begin @time include("test_Utils.jl") end
@safetestset "SimDir" begin @time include("test_Sim.jl") end
@safetestset "Atmos" begin @time include("test_Atmos.jl") end
@safetestset "OutputVar" begin @time include("test_Var.jl") end
@safetestset "CairoMakieExt" begin @time include("test_CairoMakieExt.jl") end
@safetestset "GeoMakieExt" begin @time include("test_GeoMakieExt.jl") end
#! format: on

nothing
