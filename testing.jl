using Interpolations
agrid = [0.0, 1.0]
A = [[1.0, 2.0]  [3.0, 4.0]]
vals = eachslice(A, dims = 2) |> collect
itp = interpolate((agrid,), vals, Gridded(Linear()))
