include(joinpath(@__DIR__, "ClimaExplorer.jl"))

import Bonito: route!, wait, Server

path = length(ARGS) >= 1 ? ARGS[1] : "."

# Create server
IPa = "127.0.0.1"
port = 8080
server = Server(IPa, port; proxy_url = "http://localhost:8080")

app = ClimaExplorer.BonitoApp(path)
route!(server, "/" => app)
println(server)
wait(server)
