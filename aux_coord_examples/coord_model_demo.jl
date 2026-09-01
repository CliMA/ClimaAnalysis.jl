import NCDatasets

# ---------------------------------------------------------------------------
# The data model: ONE coordinate type for both dimension and auxiliary coords.
# ---------------------------------------------------------------------------

struct Coordinate{T, N, A <: AbstractArray{T, N}}
    name::String
    dims::NTuple{N, String}   # which axes of the PARENT VARIABLE it spans
    values::A                 # array of exactly that shape
    attribs::Dict{String, Any}
end

struct Var{T, N, A <: AbstractArray{T, N}}
    name::String
    dims::NTuple{N, String}              # axis names, storage order
    data::A
    coords::Dict{String, Coordinate}
    attribs::Dict{String, Any}
end

"A dimension coordinate is the special case: 1D and spans the axis of its own name."
is_dimension_coord(c::Coordinate) = c.dims == (c.name,)

"""
The core lookup contract: given a full index `I` into `var.data`, return the
value of coordinate `cname` at that point — no matter how many axes it spans.
"""
function coordvalue(var::Var, cname::String, I::Tuple)
    c = var.coords[cname]
    # position of each coord dim within the variable's dims
    axpos = map(d -> findfirst(==(d), var.dims), c.dims)
    return c.values[map(p -> I[p], axpos)...]
end

"Label-based selection is only well-defined for dimension coordinates."
function nearest_index(var::Var, cname::String, target)
    c = var.coords[cname]
    is_dimension_coord(c) || error(
        "cannot select by '$cname': it spans $(c.dims), not just itself. " *
        "Interpolate/regrid first.",
    )
    return argmin(abs.(c.values .- target))
end

# ---------------------------------------------------------------------------
# Loader with the dims ⊆ var.dims rule doing the discovery work
# ---------------------------------------------------------------------------

function load_var(path, varname)
    NCDatasets.NCDataset(path) do nc
        v = nc[varname]
        vdims = Tuple(NCDatasets.dimnames(v))
        coords = Dict{String, Coordinate}()
        for (name, cv) in nc
            name == varname && continue
            cdims = Tuple(NCDatasets.dimnames(cv))
            # THE rule: a coordinate's dims must be a subset of the var's dims
            issubset(cdims, vdims) || continue
            coords[name] = Coordinate(
                name, cdims, Array(cv), Dict{String, Any}(cv.attrib),
            )
        end
        return Var(
            varname, vdims, Array(v), coords, Dict{String, Any}(v.attrib),
        )
    end
end

# ---------------------------------------------------------------------------
# Demo 1: ClimaDiagnostics output with topography
# ---------------------------------------------------------------------------

println("="^70)
ta = load_var("/home/kphan2/aux_coord_examples/ta_with_topography.nc", "ta")
println("ta dims: ", ta.dims, "  size: ", size(ta.data))
for (n, c) in sort(collect(ta.coords); by = first)
    kind = is_dimension_coord(c) ? "dimension" : "auxiliary"
    println(rpad(n, 12), " spans ", rpad(string(c.dims), 36), " -> ", kind)
end

# One point, all its coordinates:
I = (1, 10, 5, 2)   # (time, lon, lat, z_reference) indices
println("\npoint at index ", I, ":")
for n in ["lon", "lat", "z_reference", "z_physical"]
    println("  ", rpad(n, 12), " = ", round(coordvalue(ta, n, I); digits = 1))
end

# Selection works on dimension coords, refuses aux coords:
println("\nnearest z_reference to 3000 m -> index ",
    nearest_index(ta, "z_reference", 3000.0))
try
    nearest_index(ta, "z_physical", 3000.0)
catch e
    println("z_physical selection -> ERROR: ", e.msg)
end

# ---------------------------------------------------------------------------
# Demo 2: rasm.nc, curvilinear grid — same code, zero changes
# ---------------------------------------------------------------------------

println("\n", "="^70)
tair = load_var("/home/kphan2/aux_coord_examples/rasm.nc", "Tair")
println("Tair dims: ", tair.dims, "  size: ", size(tair.data))
for (n, c) in sort(collect(tair.coords); by = first)
    kind = is_dimension_coord(c) ? "dimension" : "auxiliary"
    println(rpad(n, 12), " spans ", rpad(string(c.dims), 36), " -> ", kind)
end

J = (100, 150, 3)   # (x, y, time) indices — x/y are bare index axes!
println("\npoint at index ", J, ": lon = ",
    round(coordvalue(tair, "xc", J); digits = 2), ", lat = ",
    round(coordvalue(tair, "yc", J); digits = 2), ", Tair = ",
    round(tair.data[J...]; digits = 2))
