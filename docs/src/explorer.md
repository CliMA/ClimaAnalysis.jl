`ClimaExplorer`
=============

`ClimaExplorer` is an interactive utility to visualize the output of `CliMA`
simulations. `ClimaExplorer` uses `ClimaAnalysis.jl` and `Bonito.jl` to start a
local server and dynamically serve plots.

How to install
===============

To install `ClimaExplorer`, first clone the `ClimaAnalysis.jl` repository
```sh
git clone https://github.com/CliMA/ClimaAnalysis.jl.git
```
Then, instantiate the `climaexplorer` environment
```sh
julia --project=ClimaAnalysis.jl/ClimaExplorer -e 'using Pkg; Pkg.instantiate()'
```

How to use
==========

Once you installed `ClimaExplorer`, you can start the server with
```sh
julia --project=ClimaAnalysis.jl/ClimaExplorer ClimaAnalysis.jl/ClimaExplorer/explorer.jl 
```
This might take a while.

You can also pass a path as an argument to `explorer.jl`. This will start `ClimaExplorer` with  
that folder (by default, the current working directory is used).

If everything went correctly, you should see a message like
```
Server:
  isrunning: true
  listen_url: http://localhost:8080
  online_url: http://localhost:9384
  http routes: 1
    / => App
  websocket routes: 0
```
This means that the server has started. Visit `http://localhost:8080` to start
inspecting the output. Note that the port number `8080` might be different.
To close the server, press `Ctrl-C`.

How to use on remote clusters
=============================

Often data resides on remote clusters. You can run `ClimaExplorer` on the
machine where the data lives and use an ssh tunnel to access the server from
your computer.

To do so, first start the server on the remote cluster (as described above). 
Next, on your local machine, forward the remote port:
```bash
ssh -f username@host -L 9384:localhost:8080 -N
```
instead of `username@host`, put your `usarename` and `hostname` (e.g.,
`pippo@login.hpc.caltech.edu`)

Next, open your browser and access `ClimaExplorer` visiting `localhost:9394`.
