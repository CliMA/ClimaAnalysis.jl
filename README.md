<h1 align="center">
  <img src="logo.svg" width="100px"> <br>
ClimaAnalysis.jl
</h1>
<p align="center">
  <strong>Analyzing and visualizing ClimaAtmos simulations</strong>
</p>

[![CI](https://github.com/CliMA/ClimaAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/CliMA/ClimaAnalysis.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://CliMA.github.io/ClimaAnalysis.jl)
[![codecov](https://codecov.io/gh/CliMA/ClimaAnalysis.jl/graph/badge.svg?token=tXO3LzS8v9)](https://codecov.io/gh/CliMA/ClimaAnalysis.jl)

`ClimaAnalysis.jl` is a Julia library to post-process and visualize `ClimaAtmos`
simulations.

Check out the [documentation](https://CliMA.github.io/ClimaAnalysis.jl) for more information and tutorials.

## Features

- Read, organize, and process NetCDF files
- Visualize heatmaps and 1D profiles with `CairoMakie`
- Visualize heatmaps on a globe with `GeoMakie`
- Apply averages and other reductions to the output variables.
- Slice variables along a given value (e.g., take the slice with altitude of 500 meters)
- Window variables within given ranges (e.g., select times between 10 and 100 days)
- Perform mathematical operations between output variables.

## ClimaAnalysis.jl Developer Guidelines

These guidelines aim to ensure consistent code quality, maintainability, and a
smooth collaborative workflow for `ClimaAnalysis.jl`. Please, read these
guidelines even if you are familiar with other CliMA packages as there may be
some differences.

### Tests and environments

We prioritize well-tested code to guarantee `ClimaAnalysis.jl` functions
reliably. Here are some principles we follow:

#### Tests are collected in the `test` folder and are exclusively there

This means that all the tests can be run with `Pkg.test()`.

#### There are no checked `Manifest.toml` files

While checking in `Manifest.toml` files ensures reproducibility, it also
introduces some nuisance, including:
- lot of git/repository noise just for "up deps";
- multiple environments that have to be managed;
- busywork to keep the manifests updated.

In this repository, we have two environments:
- project,
- documentation.

The project environment defines the test dependencies in its `extras` (to reduce
the number of environments and to avoid the "cannot merge projects" problem).

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

#### Running tests

`ClimaAnalysis.jl` defines the test dependencies directly in the main
`Project.toml`. This means that the package can be tested simply by running `]
test` in a Julia REPL, as shown below:

Start a Julia session in the `ClimaAnalysis` directory:
``` sh
julia --project
```
Enter `Pkg` mode by typing `]`. This will change the prompt. Run `test`.

When doing so, `Julia` will start a new temporary environment where the tests
are run in isolation. Tests are running checking for in-bounds and for
deprecations, and this can result in code invalidation and new precompilation.

Note, the project environment does not contain the test dependencies. Therefore,
you will find that some dependencies are missing if you try "manually" run the
test in a REPL. To solve this problem, use
[TestEnv](https://github.com/JuliaTesting/TestEnv.jl). Install `TestEnv` in your
base environment (`julia -e 'using Pkg; Pkg.add("TestEnv")'`). Then, when you
want to use the test dependencies, activate it from your REPL with `using
TestEnv; TestEnv.activate()`. This will bump you to an environment where the
test dependencies are available.

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

#### Code Formatting with `JuliaFormatter.jl`

One of the tests consists in checking that the code is uniformly formatted. We
use [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl) to achieve
consistent formatting. Here's how to use it:

You can either install in your base environment with
``` sh
julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
```
or use it from within the `TestEnv` (or base) environments (see previous section).

Then, you can format the package running:
``` julia
using JuliaFormatter; format(".")
```
or just with `format(".")` if the package is already imported.

The rules for formatting are defined in the `.JuliaFormatter.toml`.

If you are used to formatting from the command line instead of the REPL, you can
install `JuliaFormatter` in your base environment and call
``` sh
julia -e 'using JuliaFormatter; format(".")'
```
You could also define a shell alias
``` sh
alias julia_format_here="julia -e 'using JuliaFormatter; format(\".\")'"
```

> :note: Please, open an issue if you find workflow problems/friction with this
> system.

### Documentation

Documentation is generated with
[Documenter.jl](https://documenter.juliadocs.org/stable/). We strive to have
complete and up-to-date information.

To generate documentation, run
``` sh
julia --project=docs docs/make.jl
```

Please, update the documentation if you add new features or change the behavior
of existing ones.

We encourage using `jldoctest` to add and test examples in docstrings.

### Pull Request (PR) and commits

Here's how to structure your contributions effectively:

- Descriptive Title: Briefly summarize the changes your PR introduces. Commit
  titles should preferably be under 50 characters, start with a capital latter,
  and use imperative verbs (e.g., "Remove superfluous function call").
- Detailed Description: Explain the purpose of your changes. Focus on the
  intent.
- Breaking Changes: If your PR introduces breaking changes, highlight them
  clearly in the description.

Your pull request can contain one or multiple commits. In either cases, it is
important that each commit is atomic (meaning that each commit represents a
single logical change).

Please, squash commits that represent a single logical change (e.g., do not have
two commits when the second just fixes the first).

Pull requests are not merged, but _rebased_, ensuring a linear history (this is
handled automatically by GitHub).
