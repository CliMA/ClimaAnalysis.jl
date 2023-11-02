using SafeTestsets
using Test

#! format: off
@safetestset "Aqua" begin @time include("aqua.jl") end
@safetestset "SimDir" begin @time include("SimDir.jl") end
#! format: on

nothing
