module Utils

export match_nc_filename,
    squeeze, nearest_index, kwargs, seconds_to_prettystr, warp_string

import Dates
import OrderedCollections: OrderedDict

"""
    match_nc_filename(filename::String)

Return `short_name`, `period`, `reduction`, `coord_type` extracted from the filename, if
matching the expected convention.

The convention is: `shortname_(period)_reduction_(coord_type).nc`, with `period` and
`coord_type` being optional.

Examples
=========

```jldoctest
julia> match_nc_filename("bob")
```

```jldoctest
julia> match_nc_filename("ta_1d_average.nc")
("ta", "1d", "average", nothing)
```

```jldoctest
julia> match_nc_filename("ta_1d_average_pressure.nc")
("ta", "1d", "average", "pressure")
```

```jldoctest
julia> match_nc_filename("pfull_6.0m_max.nc")
("pfull", "6.0m", "max", nothing)
```

```jldoctest
julia> match_nc_filename("hu_inst.nc")
("hu", nothing, "inst", nothing)
```
"""
function match_nc_filename(filename::String)
    # Let's unpack this regular expression to find files names like "orog_inst.nc" or
    # "ta_3.0h_average.nc" and extract information from there.

    # ^: Matches the beginning of the string

    # (\w+?): Matches one or more word characters (letters, numbers, or underscore)
    # non-greedily and captures it as the first group (variable name)

    # (?:_((?:\d+(?:\.\d+)?[mMdsyh]+(?:_?\d+(?:\.\d+)?[mMdsyh]+)*)))?:
    # Optionally matches an underscore followed by the period string and captures it as the
    # second group. Period format: digit(s), optional decimal, time unit letter(s) (m, M, d,
    # s, y, h), optionally repeated with or without underscores (e.g., "1m_40s", "2d1s")

    # _: Matches the underscore before the statistic

    # ([a-zA-Z0-9]+): Matches one or more alphanumeric characters and captures it as the
    # third group (statistic)

    # (?:_([a-zA-Z]+))?: Optionally matches an underscore followed by alphabetic characters
    # (coords suffix) and captures it as the fourth group (e.g. "pressure")

    # \.nc: Matches the literal ".nc" file extension

    # $: Matches the end of the string

    re =
        r"^(\w+?)(?:_((?:\d+(?:\.\d+)?[mMdsyh]+(?:_?\d+(?:\.\d+)?[mMdsyh]+)*)))?_([a-zA-Z0-9]+)(?:_([a-zA-Z]+))?\.nc$"
    m = match(re, filename)
    if !isnothing(m)
        # m.captures returns `SubString`s (or nothing). We want to have actual `String`s (or
        # nothing) so that we can assume we have `String`s everywhere. We also take care of
        # the case where the period is matched to an empty string and return nothing instead
        return Tuple(
            (isnothing(cap) || cap == "") ? nothing : String(cap) for
            cap in m.captures
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

julia> seconds_to_prettystr(0)
"0.0s"

julia> seconds_to_prettystr(-864010)
"-10d 10s"
```
"""
function seconds_to_prettystr(seconds::Real)
    iszero(seconds) && return "0.0s"

    sign_prefix = signbit(seconds) ? "-" : ""
    seconds = abs(seconds)

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

    return sign_prefix * join(time, " ")
end

"""
    warp_string(str::AbstractString; max_width = 70)

Return a string where each line is at most `max_width` characters or less
or at most one word.

Examples
=========

```jldoctest
julia> warp_string("space", max_width = 5)
"space"

julia> warp_string("space", max_width = 4)
"space"

julia> warp_string("\\tspace    ", max_width = 4)
"space"

julia> warp_string("space space", max_width = 5)
"space\\nspace"

julia> warp_string("space space", max_width = 4)
"space\\nspace"

julia> warp_string("\\n   space  \\n  space", max_width = 4)
"space\\nspace"
```
"""
function warp_string(str::AbstractString; max_width = 70)
    return_str = ""
    current_width = 0
    for word in split(str, isspace)
        word_width = length(word)
        if word_width + current_width <= max_width
            return_str *= "$word "
            current_width += word_width + 1
        else
            # Ensure that spaces never precede newlines
            return_str = rstrip(return_str)
            return_str *= "\n$word "
            current_width = word_width + 1
        end
    end
    # Remove new line character when the first word is longer than
    # `max_width` characters and remove leading and trailing
    # whitespace
    return strip(lstrip(return_str, '\n'))
end

"""
    split_by_season(dates::AbstractArray{<: Dates.DateTime};
                    seasons = ("MAM", "JJA", "SON", "DJF")))

Return four vectors with `dates` split by seasons.

The months of the seasons are March to May, June to August, September to November, and
December to February. The order of the tuple is MAM, JJA, SON, and DJF.

Examples
=========

```jldoctest
julia> import Dates

julia> dates = [Dates.DateTime(2024, 1, 1), Dates.DateTime(2024, 3, 1), Dates.DateTime(2024, 6, 1), Dates.DateTime(2024, 9, 1)];

julia> split_by_season(dates)
([Dates.DateTime("2024-03-01T00:00:00")], [Dates.DateTime("2024-06-01T00:00:00")], [Dates.DateTime("2024-09-01T00:00:00")], [Dates.DateTime("2024-01-01T00:00:00")])
```
"""
function split_by_season(
    dates::AbstractArray{<:Dates.DateTime};
    seasons = ("MAM", "JJA", "SON", "DJF"),
)
    all(season -> season in ("MAM", "JJA", "SON", "DJF"), seasons) ||
        error("$seasons does not contain only seasons")

    MAM, JJA, SON, DJF = Vector{Dates.DateTime}(),
    Vector{Dates.DateTime}(),
    Vector{Dates.DateTime}(),
    Vector{Dates.DateTime}()

    season_to_dates =
        Dict("MAM" => MAM, "JJA" => JJA, "SON" => SON, "DJF" => DJF)

    for date in dates
        if Dates.Month(3) <= Dates.Month(date) <= Dates.Month(5)
            push!(MAM, date)
        elseif Dates.Month(6) <= Dates.Month(date) <= Dates.Month(8)
            push!(JJA, date)
        elseif Dates.Month(9) <= Dates.Month(date) <= Dates.Month(11)
            push!(SON, date)
        else
            push!(DJF, date)
        end
    end

    return tuple((season_to_dates[season] for season in seasons)...)
end

"""
    split_by_season_across_time(dates::AbstractArray{<:Dates.DateTime})

Split `dates` into vectors representing seasons, arranged in chronological order. Each
vector corresponds to a single season and the ordering of the vectors is determined by the
dates of the season. The return type is a vector of vectors of dates.

If no dates are found for a particular season, then the vector will be empty. The first
vector is guaranteed to be non-empty.

This function differs from `split_by_season` as `split_by_season` splits dates into
different seasons and ignores that dates could come from seasons in different years. In
contrast, `split_by_season_across_time` splits dates into seasons for each year.

Examples
=========

```jldoctest
julia> import Dates

julia> dates = collect(Dates.DateTime(2010, i) for i in 1:12);

julia> split_by_season_across_time(dates)
5-element Vector{Vector{Dates.DateTime}}:
 [Dates.DateTime("2010-01-01T00:00:00"), Dates.DateTime("2010-02-01T00:00:00")]
 [Dates.DateTime("2010-03-01T00:00:00"), Dates.DateTime("2010-04-01T00:00:00"), Dates.DateTime("2010-05-01T00:00:00")]
 [Dates.DateTime("2010-06-01T00:00:00"), Dates.DateTime("2010-07-01T00:00:00"), Dates.DateTime("2010-08-01T00:00:00")]
 [Dates.DateTime("2010-09-01T00:00:00"), Dates.DateTime("2010-10-01T00:00:00"), Dates.DateTime("2010-11-01T00:00:00")]
 [Dates.DateTime("2010-12-01T00:00:00")]
```
"""
function split_by_season_across_time(dates::AbstractArray{<:Dates.DateTime})
    # Dates are not necessarily sorted
    dates = sort(dates)

    # Empty case
    isempty(dates) && return Vector{Vector{eltype(dates)}}[]

    # Find the first date of the season that first(dates) belongs in
    (first_season, first_year) = find_season_and_year(first(dates))
    season_to_month = Dict("MAM" => 3, "JJA" => 6, "SON" => 9, "DJF" => 12)
    # Because the year is determined from the second month, we need to handle the case when
    # the season is DJF
    first_year = first_season == "DJF" ? first_year - 1 : first_year
    first_date_of_season =
        Dates.DateTime(first_year, season_to_month[first_season], 1)

    # Create an ordered dict to map between season and year to vector of dates
    season_and_year2dates = OrderedDict{
        Tuple{typeof(first_season), typeof(first_year)},
        Vector{eltype(dates)},
    }()
    # Need to iterate because some seasons can be empty and we want empty vectors for that
    curr_date = first_date_of_season
    while curr_date <= dates[end]
        (season, year) = find_season_and_year(curr_date)
        season_and_year2dates[(season, year)] = typeof(curr_date)[]
        curr_date += Dates.Month(3) # season change every 3 months
    end

    # Add dates to the correct vectors in season_and_year2dates
    for date in dates
        (season, year) = find_season_and_year(date)
        push!(season_and_year2dates[(season, year)], date)
    end
    return collect(values(season_and_year2dates))
end

"""
    split_by_month(dates::AbstractArray{<:Dates.DateTime})

Split `dates` into vectors representing months, arranged in chronological order. Each
vector corresponds to a single month and the ordering of the vectors is determined by the
dates of the season. The return type is a vector of vectors of dates.

If no dates are found for a particular season, then the vector will be empty.

Examples
=========

```jldoctest
julia> import Dates

julia> dates = vcat(collect(Dates.DateTime(2010, i) for i in 1:11), [Dates.DateTime(2011, 4)]);

julia> split_by_month(dates)
12-element Vector{Vector{Dates.DateTime}}:
 [Dates.DateTime("2010-01-01T00:00:00")]
 [Dates.DateTime("2010-02-01T00:00:00")]
 [Dates.DateTime("2010-03-01T00:00:00")]
 [Dates.DateTime("2010-04-01T00:00:00"), Dates.DateTime("2011-04-01T00:00:00")]
 [Dates.DateTime("2010-05-01T00:00:00")]
 [Dates.DateTime("2010-06-01T00:00:00")]
 [Dates.DateTime("2010-07-01T00:00:00")]
 [Dates.DateTime("2010-08-01T00:00:00")]
 [Dates.DateTime("2010-09-01T00:00:00")]
 [Dates.DateTime("2010-10-01T00:00:00")]
 [Dates.DateTime("2010-11-01T00:00:00")]
 []
```
"""
function split_by_month(dates::AbstractArray{<:Dates.DateTime})
    # Dates are not necessarily sorted
    dates = sort(dates)

    # Empty case
    isempty(dates) && return Vector{Vector{eltype(dates)}}[]

    # Each integer represent a month
    month2dates = OrderedDict{Int, Vector{eltype(dates)}}()
    for i in 1:12
        month2dates[i] = []
    end

    for date in dates
        month = Dates.month(date)
        push!(month2dates[month], date)
    end

    return collect(values(month2dates))
end

"""
    find_season_and_year(date::Dates.DateTime)

Return a tuple of the year and season belong to `date`. The variable `year` is
an integer and `season` is a string.

The months of the seasons are March to May, June to August, September to
November, and December to February.

The year is determined by the second month of the season. For example, if the
months are December 2009, January 2010, and February 2010, then the year of the
season is 2010.
"""
function find_season_and_year(date::Dates.DateTime)
    if Dates.Month(3) <= Dates.Month(date) <= Dates.Month(5)
        return ("MAM", Dates.year(date))
    elseif Dates.Month(6) <= Dates.Month(date) <= Dates.Month(8)
        return ("JJA", Dates.year(date))
    elseif Dates.Month(9) <= Dates.Month(date) <= Dates.Month(11)
        return ("SON", Dates.year(date))
    else
        corrected_year =
            Dates.month(date) == 12 ? Dates.year(date) + 1 : Dates.year(date)
        return ("DJF", corrected_year)
    end
end

"""
    _isequispaced(arr::Vector)

Return whether the array is equispaced or not.

Examples
=========

```jldoctest
julia> Utils._isequispaced([1.0, 2.0, 3.0])
true

julia> Utils._isequispaced([0.0, 2.0, 3.0])
false
```
"""
function _isequispaced(arr::Vector)
    return all(diff(arr) .≈ arr[begin + 1] - arr[begin])
end

"""
    time_to_date(start_date::Dates.DateTime, time::AbstractFloat)

Convert the given time to a calendar date starting from `start_date`.

Examples
=========

```jldoctest
julia> import Dates

julia> Utils.time_to_date(Dates.DateTime(2013, 7, 1, 12), 86400.0)
2013-07-02T12:00:00

julia> Utils.time_to_date(Dates.DateTime(2013, 7, 1, 12), 3600.0)
2013-07-01T13:00:00

julia> Utils.time_to_date(Dates.DateTime(2013, 7, 1, 12), 60.0)
2013-07-01T12:01:00

julia> Utils.time_to_date(Dates.DateTime(2013, 7, 1, 12), 1.0)
2013-07-01T12:00:01
```
"""
function time_to_date(start_date::Dates.DateTime, time::AbstractFloat)
    # We go through milliseconds to allow fractions of a second (otherwise, Second(0.8)
    # would fail). Milliseconds is the level of resolution that one gets when taking the
    # difference between two DateTimes. In addition to this, we add a round to account for
    # floating point errors. If the floating point error is small enough, round will correct
    # it.
    milliseconds = Dates.Millisecond.(round.(1_000 * time))
    return start_date + milliseconds
end

"""
    date_to_time(reference_date::Dates.DateTime, date::Dates.DateTime)

Convert the given calendar date to a time (in seconds) where t=0 is `reference_date`.

Examples
=========

```jldoctest
julia> import Dates

julia> Utils.date_to_time(Dates.DateTime(2013, 7, 1, 12), Dates.DateTime(2013, 7, 2, 12))
86400.0

julia> Utils.date_to_time(Dates.DateTime(2013, 7, 1, 12), Dates.DateTime(2013, 7, 1, 13))
3600.0

julia> Utils.date_to_time(Dates.DateTime(2013, 7, 1, 12), Dates.DateTime(2013, 7, 1, 12, 1))
60.0

julia> Utils.date_to_time(Dates.DateTime(2013, 7, 1, 12), Dates.DateTime(2013, 7, 1, 12, 0, 1))
1.0
```
"""
function date_to_time(reference_date::Dates.DateTime, date::Dates.DateTime)
    return period_to_seconds_float(date - reference_date)
end

"""
    period_to_seconds_float(period::Dates.Period)

Convert the given `period` to seconds in Float64.

Examples
=========

```jldoctest
julia> import Dates

julia> Utils.period_to_seconds_float(Dates.Millisecond(1))
0.001

julia> Utils.period_to_seconds_float(Dates.Second(1))
1.0

julia> Utils.period_to_seconds_float(Dates.Minute(1))
60.0

julia> Utils.period_to_seconds_float(Dates.Hour(1))
3600.0

julia> Utils.period_to_seconds_float(Dates.Day(1))
86400.0

julia> Utils.period_to_seconds_float(Dates.Week(1))
604800.0

julia> Utils.period_to_seconds_float(Dates.Month(1))
2.629746e6

julia> Utils.period_to_seconds_float(Dates.Year(1))
3.1556952e7
```
"""
function period_to_seconds_float(period::Dates.Period)
    # See https://github.com/JuliaLang/julia/issues/55406
    period isa Dates.OtherPeriod &&
        (period = Dates.Second(Dates.Day(1)) * Dates.days(period))
    return period / Dates.Second(1)
end

"""
    _data_at_dim_vals(data, dim_arr, dim_idx, vals)

Return a view of `data` by slicing along `dim_idx`. The slices are indexed by the indices
corresponding to values in `dim_arr` closest to the values in `vals`.

Examples
=========

```jldoctest
julia> data = [[1, 4, 7]  [2, 5, 8]  [3, 6, 9]];

julia> dim_arr = [1.0, 2.0, 4.0];

julia> dim_idx = 2;

julia> vals = [1.1, 4.0];

julia> Utils._data_at_dim_vals(data, dim_arr, dim_idx, vals)
3×2 view(::Matrix{Int64}, :, [1, 3]) with eltype Int64:
 1  3
 4  6
 7  9
```
"""
function _data_at_dim_vals(data, dim_arr, dim_idx, vals)
    nearest_indices = map(val -> nearest_index(dim_arr, val), vals)
    return selectdim(data, dim_idx, nearest_indices)
end

"""
    _recursive_merge(x::AbstractDict...)

Recursively merge nested dictionaries together. In the case of keys duplicated among the
arguments, the rightmost argument that owns the key gets its value stored in the result.

The function is helpful when dealing with keyword arguments when plotting and you want to
have default keyword arguments.

See: https://discourse.julialang.org/t/multi-layer-dict-merge/27261/6
"""
_recursive_merge(x::AbstractDict...) = merge(_recursive_merge, x...)
_recursive_merge(x...) = x[end]
end
