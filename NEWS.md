ClimaAnalysis.jl Release Notes
===============================

main
------
- Add support for evaluating `OutputVar`s onto arbitrary target points (with
  multilinear interpolation).
- `average` operations now ignore `NaN`s by default.
- Add `has_*` methods to query whether a `Var` has a given dimension (e.g., `z`).
- Support `Makie` backends besides `CairoMakie`.
- Add methods to get the range and units of a given dimension in `Var`. 

v0.5.6
------
- Fix finding variables with name like `clwup_1m_40s_inst.nc` (composed period).
- Add support for weighted averages in `average_lat`.

v0.5.5
------
- Fix reading `NetCDF` files with dimensions in incorrect order.

v0.5.4
------
- Added support for extraction dimension from functions, such as `times`.
- Reorganized internal modules so that each file is a module.

v0.5.3
------
- Add `Visualize.contour2D_on_globe!` for discrete contours.
