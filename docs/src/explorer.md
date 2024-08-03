`ClimaExplorer`
=============

`ClimaExplorer` is an interactive utility to visualize the output of `CliMA`
simulations. `ClimaExplorer` uses `ClimaAnalysis.jl` and `Bonito.jl` to start a
local server and dynamically serve plots.

How to install
===============

To install `ClimaExplorer`, first clone the `ClimaAnalysis.jl` repo
```sh
git clone https://github.com/CliMA/ClimaAnalysis.jl.git
```
Then, instantiate the `climaexplorer` environment
```sh
julia --project=ClimaAnalysis.jl/climaexplorer -e 'using Pkg; Pkg.instantiate()'
```

How to use
==========

Once you installed `ClimaExplorer`, you can start the server with
```sh
julia --project=ClimaAnalysis.jl/climaexplorer ClimaAnalysis.jl/climaexplorer/explorer.jl PATH_OF_OUTPUT 
```
where `PATH_OF_OUTPUT` is the path where the output files live.
The first time you do so, this will take a little while.

If everything went correctly, you should see a message like
```
Server:
  isrunning: true
  listen_url: http://localhost:8080
  online_url: http://localhost:8080
  http routes: 1
    / => App
  websocket routes: 0
```
This means that the server has started. Visit `http://localhost:8080` to start
inspecting the output. Note that the port number `8080` might be different.
To close the server, press `Ctrl-C`.
