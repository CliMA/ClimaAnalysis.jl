using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "Docstrings" begin @time include("doctest.jl") end
@safetestset "Format" begin @time include("format.jl") end

@safetestset "Utils" begin @time include("Utils.jl") end
@safetestset "SimDir" begin @time include("SimDir.jl") end
@safetestset "OutputVar" begin @time include("OutputVar.jl") end
@safetestset "CairoMakieExt" begin @time include("CairoMakieExt.jl") end
#! format: on

nothing
