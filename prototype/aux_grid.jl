# Prototype of the axes-vs-coordinates model (GRID_DESIGN.md Phase 5, HANDOFF.md §5):
# a Grid holds one 1D dimension coordinate per data axis, plus any number of
# auxiliary coordinates, each recording which data axes it spans.
#
# Deliberately simplified vs the real design (no units/bounds fields, plain Dicts,
# no conventional-name resolution, printing inside library code): the goal is to
# exercise every operation whose semantics auxiliary coordinates could break.

module AuxGrid

using NCDatasets

export Coord, Grid, Var, read_var, subset, slice, average, weighted_average,
    coordinate_field, coordvalue, nearest_point, interp_to_coord_levels,
    axis_names, coord_names, coord, is_dimension_coord

# ------------------------------------------------------------------ types

"""
A coordinate: values laid over one or more axes of the parent variable.
A *dimension coordinate* is the special case `spans == (name,)` with 1D values.
`spans == ()` (0-dim values) is a scalar coordinate left behind after integer
indexing dropped all the axes it spanned.
"""
struct Coord{N, A <: AbstractArray, T}
    name::String
    spans::NTuple{N, String}
    values::A
    attributes::Dict{String, T}
end

function Coord(
    name::AbstractString,
    spans::NTuple{N, <:AbstractString},
    values::AbstractArray;
    attributes = Dict{String, Any}(),
) where {N}
    ndims(values) == N || error(
        "coordinate '$name': values have $(ndims(values)) dims but spans lists $N axes",
    )
    Coord(String(name), NTuple{N, String}(spans), values, attributes)
end

is_dimension_coord(c::Coord) = c.spans == (c.name,) && ndims(c.values) == 1

"""
`dims`: one 1D dimension coordinate per data axis, in data-axis order.
`coords`: auxiliary coordinates; each spans a subset of the axis names.
"""
struct Grid{D <: Tuple{Vararg{Coord}}, C <: Tuple{Vararg{Coord}}}
    dims::D
    coords::C

    function Grid(dims::Tuple{Vararg{Coord}}, coords::Tuple{Vararg{Coord}} = ())
        all(is_dimension_coord, dims) ||
            error("grid axes must be 1D dimension coordinates (spans == (name,))")
        names = [d.name for d in dims]
        allunique(vcat(names, [c.name for c in coords])) ||
            error("coordinate names must be unique")
        for c in coords
            issubset(c.spans, names) ||
                error("auxiliary coordinate '$(c.name)' spans $(c.spans) ⊄ axes $names")
            expected = Tuple(length(dims[findfirst(==(s), names)].values) for s in c.spans)
            size(c.values) == expected || error(
                "auxiliary coordinate '$(c.name)': size $(size(c.values)) ≠ $expected from its axes $(c.spans)",
            )
        end
        new{typeof(dims), typeof(coords)}(dims, coords)
    end
end

axis_names(g::Grid) = map(d -> d.name, g.dims)
coord_names(g::Grid) = map(c -> c.name, g.coords)
Base.ndims(g::Grid) = length(g.dims)
Base.size(g::Grid) = map(d -> length(d.values), g.dims)

function axis_index(g::Grid, name::AbstractString)
    i = findfirst(d -> d.name == name, g.dims)
    i === nothing && error("no axis '$name'; axes: $(join(axis_names(g), ", "))")
    i
end

"Look up any coordinate (dimension or auxiliary) by name."
function coord(g::Grid, name::AbstractString)
    for d in g.dims
        d.name == name && return d
    end
    for c in g.coords
        c.name == name && return c
    end
    error(
        "no coordinate '$name'; axes: $(join(axis_names(g), ", ")); " *
        "auxiliary: $(join(coord_names(g), ", "))",
    )
end

Base.getindex(g::Grid, name::AbstractString) = coord(g, name).values

struct Var{G <: Grid, A <: AbstractArray, T}
    attributes::Dict{String, T}
    grid::G
    data::A

    function Var(attributes::Dict{String, T}, grid::Grid, data::AbstractArray) where {T}
        size(data) == size(grid) ||
            error("data size $(size(data)) ≠ grid size $(size(grid))")
        new{typeof(grid), typeof(data), T}(attributes, grid, data)
    end
end

coord(var::Var, name) = coord(var.grid, name)
axis_index(var::Var, name) = axis_index(var.grid, name)

"Value of any coordinate at a full data index, however many axes it spans."
function coordvalue(var::Var, name, I::CartesianIndex)
    c = coord(var.grid, name)
    pos = map(s -> axis_index(var.grid, s), c.spans)
    c.values[map(p -> I[p], pos)...]
end

# ------------------------------------------------------------------ show

_span_str(v::AbstractVector) = isempty(v) ? "(empty)" : "$(first(v)) … $(last(v))"

function Base.show(io::IO, ::MIME"text/plain", var::Var)
    println(io, "Var $(size(var.data)) ::$(eltype(var.data))")
    println(io, "  dimension coordinates:")
    for d in var.grid.dims
        u = get(d.attributes, "units", "")
        println(io, "    ", rpad(d.name, 13), lpad(length(d.values), 5), "  ",
            rpad(_span_str(d.values), 34), u)
    end
    isempty(var.grid.coords) || println(io, "  auxiliary coordinates:")
    for c in var.grid.coords
        u = get(c.attributes, "units", "")
        println(io, "    ", rpad(c.name, 13), " spans ", rpad(string(c.spans), 32),
            "size ", size(c.values), "  ", u)
    end
end

# ------------------------------------------------- indexing / subsetting

"""
    subset(var; lon = 3, time = 1:10, ...)

Positional indexing by axis name: integer drops the axis, range/vector keeps it
(Base-array semantics, as in GRID_DESIGN structural indexing). Each auxiliary
coordinate is subset along the axes it shares and loses axes dropped by integer
indices (3D → 2D → … → 0-dim scalar coordinate). A dropped axis leaves its own
coordinate value behind as a scalar coordinate (xarray semantics) — the
principled version of today's `attributes["slice_lon"]` bookkeeping.
"""
function subset(var::Var; kw...)
    g = var.grid
    inds = Any[(:) for _ in 1:ndims(g)]
    for (k, v) in kw
        inds[axis_index(g, String(k))] = v
    end
    newdims = Tuple(
        Coord(d.name, d.spans, d.values[inds[i]], d.attributes)
        for (i, d) in enumerate(g.dims) if !(inds[i] isa Integer)
    )
    scalars = Tuple(
        Coord(d.name, (), fill(d.values[inds[i]]), d.attributes)
        for (i, d) in enumerate(g.dims) if inds[i] isa Integer
    )
    newcoords = (map(c -> _subset_coord(g, c, inds), g.coords)..., scalars...)
    Var(copy(var.attributes), Grid(newdims, newcoords), var.data[inds...])
end

function _subset_coord(g::Grid, c::Coord, inds)
    cinds = map(s -> inds[axis_index(g, s)], c.spans)
    newspans = Tuple(s for (s, ci) in zip(c.spans, cinds) if !(ci isa Integer))
    v = c.values[cinds...]
    v isa AbstractArray || (v = fill(v))  # all spanned axes dropped -> scalar coordinate
    Coord(c.name, newspans, v, c.attributes)
end

"""
    slice(var; lon = -87.4, ...)

Nearest-value selection — dimension coordinates only. Selection by an auxiliary
coordinate has no well-defined axis index; the error says what to use instead.
"""
function slice(var::Var; kw...)
    for (k, target) in kw
        c = coord(var.grid, String(k))
        is_dimension_coord(c) || error(
            "cannot slice by '$(c.name)': it spans $(c.spans), not just itself. " *
            "Use nearest_point (curvilinear) or interp_to_coord_levels (terrain-following) instead.",
        )
        i = argmin(abs.(c.values .- target))
        var = subset(var; NamedTuple{(Symbol(k),)}((i,))...)
    end
    var
end

# ------------------------------------------------------------- reductions

"""
    average(var; over = ("lon", "lat"))

Unweighted mean over the named axes, skipping `missing`. Reduced axes are
dropped; auxiliary coordinates spanning a reduced axis are dropped too
(xarray semantics — their values no longer describe any remaining point).
"""
average(var::Var; over) = weighted_average(var, ones(size(var.data)); over)

function weighted_average(var::Var, w::AbstractArray; over)
    over isa Union{AbstractString, Symbol} && (over = (over,))
    names = map(String, Tuple(over))
    axs = map(n -> axis_index(var.grid, n), names)
    present = map(x -> !ismissing(x) && !isnan(x), var.data)
    clean = ifelse.(present, float.(coalesce.(var.data, 0.0)), 0.0)
    num = sum(clean .* w; dims = axs)
    den = sum(w .* present; dims = axs)
    Var(copy(var.attributes), _drop_axes(var.grid, names), dropdims(num ./ den; dims = axs))
end

function _drop_axes(g::Grid, names)
    keep(c) = !any(in(names), c.spans)
    dropped = [c.name for c in g.coords if !keep(c)]
    isempty(dropped) || println(
        "  [reduce] dropping auxiliary coordinate(s) ",
        join(dropped, ", "), ": spanned a reduced axis",
    )
    Grid(
        Tuple(d for d in g.dims if d.name ∉ names),
        Tuple(c for c in g.coords if keep(c)),
    )
end

# ------------------------------------------------- permutedims / cat

"Reorder data axes. Auxiliary coordinate arrays are untouched: each carries its own axis order in `spans`."
function Base.permutedims(var::Var, perm)
    p = collect(map(n -> axis_index(var.grid, String(n)), perm))
    sort(p) == 1:ndims(var.grid) || error("perm must name every axis exactly once")
    Var(
        copy(var.attributes),
        Grid(Tuple(var.grid.dims[i] for i in p), var.grid.coords),
        permutedims(var.data, p),
    )
end

_close(a, b) = eltype(a) <: Number ? isapprox(a, b) : isequal(a, b)

"""
    cat(vars...; along = "time")

Concatenate along a named axis. Auxiliary coordinates spanning that axis are
concatenated along their own corresponding array axis; the rest must match.
"""
function Base.cat(vars::Var...; along::AbstractString)
    g1 = first(vars).grid
    k = axis_index(g1, along)
    all(v -> axis_names(v.grid) == axis_names(g1), vars) || error("axis names differ")
    newdims = ntuple(ndims(g1)) do i
        ds = [v.grid.dims[i] for v in vars]
        if i == k
            Coord(along, (along,), reduce(vcat, [d.values for d in ds]), ds[1].attributes)
        else
            all(d -> _close(d.values, ds[1].values), ds) ||
                error("dimension coordinate '$(ds[1].name)' differs between variables")
            ds[1]
        end
    end
    newcoords = map(g1.coords) do c
        cs = [coord(v.grid, c.name) for v in vars]
        if along in c.spans
            ax = findfirst(==(along), c.spans)
            Coord(c.name, c.spans, cat([x.values for x in cs]...; dims = ax), c.attributes)
        else
            all(x -> _close(x.values, c.values), cs) ||
                error("auxiliary coordinate '$(c.name)' differs between variables; refusing to cat")
            c
        end
    end
    Var(
        copy(first(vars).attributes),
        Grid(newdims, newcoords),
        cat([v.data for v in vars]...; dims = k),
    )
end

# ----------------------------------------------- coordinate consumption

"""
    coordinate_field(var, name)

The coordinate broadcast to the full data shape — one generic primitive that
works identically for dimension and auxiliary coordinates. Feeds weights,
masks, and verification.
"""
function coordinate_field(var::Var, name)
    c = coord(var.grid, name)
    pos = collect(map(s -> axis_index(var.grid, s), c.spans))
    vals = length(c.spans) > 1 ? permutedims(c.values, sortperm(pos)) : c.values
    shape = ntuple(i -> i in pos ? size(var.data, i) : 1, ndims(var.data))
    out = Array{eltype(vals)}(undef, size(var.data))
    out .= reshape(vals, shape)
    out
end

"""
    nearest_point(var; xc = 118.3, yc = 57.2)

Nearest-neighbor lookup over auxiliary coordinates sharing the same spans
(curvilinear grids). Selection *on* aux coords is possible — it is an explicit
search returning indices, not axis indexing. (Real version: KD-tree, as xoak.)
"""
function nearest_point(var::Var; kw...)
    cs = [coord(var.grid, String(k)) for (k, _) in kw]
    spans = cs[1].spans
    all(c -> c.spans == spans, cs) ||
        error("nearest_point coordinates must share the same spans")
    d2 = zeros(size(cs[1].values))
    for ((_, target), c) in zip(pairs(kw), cs)
        d2 .+= (c.values .- target) .^ 2
    end
    I = argmin(d2)
    subset(var; NamedTuple{map(Symbol, spans)}(Tuple(I))...)
end

"""
    interp_to_coord_levels(var, cname, targets; along)

Interpolate data onto fixed values of the auxiliary coordinate `cname`
(e.g. terrain-following levels → fixed physical altitudes): per column, 1D
linear interpolation of data against the coordinate along axis `along`.
That axis is replaced by a new dimension coordinate named `cname` holding
`targets`; the consumed auxiliary coordinate disappears (it became an axis).
NaN where a column does not reach the target (e.g. below terrain).

The analog of Atmos.to_pressure_coordinates / cf_xarray decode_vertical_coords.
"""
function interp_to_coord_levels(var::Var, cname, targets::AbstractVector; along::AbstractString)
    g = var.grid
    c = coord(g, cname)
    along in c.spans || error("'$cname' does not span axis '$along'")
    k = axis_index(g, along)
    nk = size(var.data, k)
    out = fill(NaN, ntuple(i -> i == k ? length(targets) : size(var.data, i), ndims(var.data)))
    xs = Vector{Float64}(undef, nk)
    ys = Vector{Float64}(undef, nk)
    for J in CartesianIndices(ntuple(i -> axes(var.data, i < k ? i : i + 1), ndims(var.data) - 1))
        for j in 1:nk
            I = _insert(J, k, j)
            xs[j] = coordvalue(var, cname, I)
            ys[j] = var.data[I]
        end
        for (t, target) in enumerate(targets)
            out[_insert(J, k, t)] = _lininterp(xs, ys, target)
        end
    end
    newdims = ntuple(
        i -> i == k ? Coord(String(cname), (String(cname),), float.(collect(targets)), c.attributes) : g.dims[i],
        ndims(g),
    )
    newcoords = Tuple(cc for cc in g.coords if cc.name != String(cname) && along ∉ cc.spans)
    Var(copy(var.attributes), Grid(newdims, newcoords), out)
end

_insert(J::CartesianIndex, k, j) =
    CartesianIndex(ntuple(i -> i < k ? J[i] : (i == k ? j : J[i - 1]), length(J) + 1))

function _lininterp(xs, ys, x)
    if first(xs) > last(xs)
        xs, ys = reverse(xs), reverse(ys)
    end
    (x < first(xs) || x > last(xs)) && return NaN
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    ys[i] + t * (ys[i + 1] - ys[i])
end

# ------------------------------------------------------------- discovery

const COORD_STANDARD_NAMES =
    Set(["longitude", "latitude", "time", "altitude", "height", "air_pressure", "depth"])
const COORD_UNITS =
    Set(["degrees_east", "degrees_north", "degree_east", "degree_north", "degreesE", "degreesN"])
const COORD_NAME_CONVENTIONS = Set([
    "lon", "long", "longitude", "lat", "latitude", "time", "t", "date",
    "z", "z_reference", "z_physical", "height", "pfull", "pressure_level",
])

"""
    read_var(path, varname)

Load a variable and discover its coordinates.

Dimension coordinates: the file variable with the same name as each data
dimension; a dimension without one becomes a bare index axis `1:n` (rasm x/y).

Auxiliary coordinates, by fallback chain:
  1. names listed in the variable's CF `coordinates` attribute;
  2. CF metadata signals: standard_name / axis attribute / coordinate-like units;
  3. name conventions — required for ClimaDiagnostics output, which links
     z_physical by neither `coordinates` nor standard_name.
Hard filters first: candidate dims ⊆ variable dims, not a bounds variable, not
the variable itself. Dims-subset alone is NOT a promotion signal: a sibling
data variable on the same axes must stay a data variable.
"""
function read_var(path, varname)
    NCDataset(path) do ds
        v = ds[varname]
        dnames = dimnames(v)
        dims = map(dnames) do dn
            if haskey(ds, dn)
                Coord(dn, (dn,), Array(ds[dn]); attributes = Dict{String, Any}(ds[dn].attrib))
            else
                Coord(dn, (dn,), 1:ds.dim[dn];
                    attributes = Dict{String, Any}("note" => "bare index axis (no coordinate variable in file)"))
            end
        end
        boundsvars = Set{String}()
        for n in keys(ds)
            b = get(ds[n].attrib, "bounds", nothing)
            b === nothing || push!(boundsvars, b)
        end
        listed = split(get(v.attrib, "coordinates", ""))
        aux = Coord[]
        for n in keys(ds)
            (n == varname || n in dnames) && continue
            if n in boundsvars
                println("  [discovery] skipping '$n': bounds variable")
                continue
            end
            cdims = dimnames(ds[n])
            issubset(cdims, dnames) || continue
            reason = _promotion_reason(n, ds[n], listed)
            if reason === nothing
                println("  [discovery] NOT promoting '$n': dims ⊆ var dims but no ",
                    "coordinate signal (a greedy subset-only rule would have taken it)")
                continue
            end
            println("  [discovery] promoting '$n' as auxiliary coordinate ($reason)")
            push!(aux, Coord(n, Tuple(cdims), Array(ds[n]);
                attributes = Dict{String, Any}(ds[n].attrib)))
        end
        Var(Dict{String, Any}(v.attrib), Grid(Tuple(dims), Tuple(aux)), Array(v))
    end
end

function _promotion_reason(name, ncvar, listed)
    name in listed && return "CF coordinates attribute"
    a = ncvar.attrib
    get(a, "standard_name", "") in COORD_STANDARD_NAMES && return "standard_name"
    haskey(a, "axis") && return "axis attribute"
    u = get(a, "units", "")
    (u in COORD_UNITS || startswith(u, "seconds since") || startswith(u, "days since")) &&
        return "coordinate-like units"
    lowercase(name) in COORD_NAME_CONVENTIONS && return "name convention"
    nothing
end

end # module
