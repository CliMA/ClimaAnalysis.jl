<h1 align="center">
  <img src="logo.svg" width="100px"> <br>
ClimaAnalysis.jl
</h1>
<p align="center">
  <strong>Analyzing and visualizing ClimaAtmos simulations</strong>
</p>
[![Test](https://github.com/Sbozzolo/ClimaAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Sbozzolo/ClimaAnalysis.jl/actions/workflows/CI.yml)
[![Docs](https://github.com/Sbozzolo/ClimaAnalysis.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/Sbozzolo/ClimaAnalysis.jl/actions/workflows/Documentation.yml)

`ClimaAnalysis.jl` is a Julia library to post-process and visualize `ClimaAtmos`
simulations.

Check out the documentation for more information and tutorials.

## Features

- Read, organize, and process NetCDF files
- Visualize heatmaps and 1D profiles with `CairoMakie`
- Apply averages and other reductions to the output variables.
- Slice variables along a given value (e.g., take the slice with altitude of 500 meters)
