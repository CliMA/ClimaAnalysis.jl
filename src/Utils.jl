module Utils

export match_nc_filename, squeeze, nearest_index, kwargs, seconds_to_prettystr

"""
    match_nc_filename(filename::String)

Return `short_name`, `period`, `reduction` extracted from the filename, if matching the
expected convention.

The convention is: `shortname_(period)_reduction.nc`, with `period` being optional.

Examples
=========

```jldoctest
julia> match_nc_filename("bob")
```

```jldoctest
julia> match_nc_filename("ta_1d_average.nc")
("ta", "1d", "average")
```

```jldoctest
julia> match_nc_filename("pfull_6.0min_max.nc")
("pfull", "6.0min", "max")
```

```jldoctest
julia> match_nc_filename("hu_inst.nc")
("hu", nothing, "inst")
```
"""
function match_nc_filename(filename::String)
    # Let's unpack this regular expression to find files names like "orog_inst.nc" or
    # "ta_3.0h_average.nc" and extract information from there.

    # ^ $: mean match the entire string
    # (\w+?): the first capturing group, matching any word non greedily
    # _: matches this literal character
    # (?>([a-zA-Z0-9\.]*)_)?: an optional group (it doesn't exist for _inst variables)
    #                         ?> means that we don't want to capture the outside group
    #                         the inside group is any combinations of letters/numbers,
    #                         and the literal character ., followed by the _.
    #                         We capture the combination of characters because that's
    #                         the reduction
    # (\w+): Again, any word
    # \.nc: file extension has to be .nc
    re = r"^(\w+?)_(?>([a-zA-Z0-9\.]*)_)?(\w*)\.nc$"
    m = match(re, filename)
    if !isnothing(m)
        # m.captures returns `SubString`s (or nothing). We want to have actual `String`s (or
        # nothing) so that we can assume we have `String`s everywhere.
        return Tuple(
            isnothing(cap) ? nothing : String(cap) for cap in m.captures
        )
    else
        return nothing
    end
end

"""
    squeeze(A :: AbstractArray; dims)

Return an array that has no dimensions with size 1.

When an iterable `dims` is passed, only try to squeeze the given `dim`ensions.

Examples
=========

```jldoctest
julia> A = [[1 2] [3 4]];

julia> size(A)
(1, 4)

julia> A_squeezed = squeeze(A);

julia> size(A_squeezed)
(4,)

julia> A_not_squeezed = squeeze(A; dims = (2, ));

julia> size(A_not_squeezed)
(1, 4)
```
"""
function squeeze(A::AbstractArray; dims = nothing)
    isnothing(dims) && (dims = Tuple(1:length(size(A))))

    # TODO: (Refactor)
    #
    # Find a cleaner way to identify `keepdims`

    dims_to_drop = Tuple(
        dim for (dim, len) in enumerate(size(A)) if dim in dims && len == 1
    )
    keepdims = Tuple(
        len for (dim, len) in enumerate(size(A)) if !(dim in dims_to_drop)
    )
    # We use reshape because of
    # https://stackoverflow.com/questions/52505760/dropping-singleton-dimensions-in-julia
    return reshape(A, keepdims)
end

"""
    nearest_index(A::AbstractArray, val)

Return the index in `A` closest to the given `val`.

Examples
=========

```jldoctest
julia> A = [-1, 0, 1, 2, 3, 4, 5];

julia> nearest_index(A, 3)
5

julia> nearest_index(A, 0.1)
2
```
"""
function nearest_index(A::AbstractArray, val)
    val < minimum(A) && return findmin(A)[2]
    val > maximum(A) && return findmax(A)[2]
    return findmin(A -> abs(A - val), A)[2]
end

"""
    kwargs(; kwargs...)

Convert keyword arguments in a dictionary that maps `Symbol`s to values.

Useful to pass keyword arguments to different constructors in a function.

Examples
=========

```jldoctest
julia> kwargs(a = 1)
pairs(::NamedTuple) with 1 entry:
  :a => 1
```
"""
kwargs(; kwargs...) = kwargs

"""
    seconds_to_prettystr(seconds::Real)

Convert the given `seconds` into a string with rich time information.

One year is defined as having 365 days.

Examples
=========

```jldoctest
julia> seconds_to_prettystr(10)
"10s"

julia> seconds_to_prettystr(600)
"10m"

julia> seconds_to_prettystr(86400)
"1d"

julia> seconds_to_prettystr(864000)
"10d"

julia> seconds_to_prettystr(864010)
"10d 10s"

julia> seconds_to_prettystr(24 * 60 * 60 * 365 + 1)
"1y 1s"
```
"""
function seconds_to_prettystr(seconds::Real)
    time = String[]

    years, rem_seconds = divrem(seconds, 24 * 60 * 60 * 365)
    days, rem_seconds = divrem(rem_seconds, 24 * 60 * 60)
    hours, rem_seconds = divrem(rem_seconds, 60 * 60)
    minutes, seconds = divrem(rem_seconds, 60)

    # At this point, days, hours, minutes, seconds have to be integers.
    # Let us force them to be such so that we can have a consistent string output.
    years, days, hours, minutes = map(Int, (years, days, hours, minutes))

    years > 0 && push!(time, "$(years)y")
    days > 0 && push!(time, "$(days)d")
    hours > 0 && push!(time, "$(hours)h")
    minutes > 0 && push!(time, "$(minutes)m")
    seconds > 0 && push!(time, "$(seconds)s")

    return join(time, " ")
end

end
