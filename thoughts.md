Future syntax:
var |> split(Season()) |> apply(mean) |> merge()
The result of var |> split(Season()) is a SplitedVar
Then, |> apply(mean) internally does a mean
Finally, merge combine all the results together
Add support for dropdims for singleton dimensions

Questions:
How to be efficient with memory since I won't know the size of things before I
do it unless I give myself some restrictions?
What will I support? Only reductions?
I will not support map operations and only reduction
(map operations can be done as if)
Any dimension can be splitted, but will be support multiple splitting (e.g.
split along time and lon, apply reduction that along those two dimensions, and
combine the final result)
New type called SplittedVar (goal is to reduce memory usage, so will try not
to allocate at this step).
Reductions like mean, sum, var, etc. will be supported in the sense and any
function of the form f(data; dims = ...) will be supported
New dimension after splitting will take the first element of the groups
The final result should be a single OutputVar after split, apply, and merge

TODO: Not sure if I should support operations like cumsum or diff (seem like
it is not needed)


Comments

Implementation notes:
#
1. Piping mechanics: split(splitter), apply(f), and combine() must all return
   functions so |> works. E.g., split(Season()) returns OutputVar -> SplittedVar.
#
2. Base.split conflict: dispatching on AbstractSplitter avoids ambiguity with
   Base.split (for strings), but do not re-export split carelessly.
#
3. SplittedVar must carry:
   - groups::Vector{OutputVar}  (non-empty only)
   - group_first_coords::Vector (first dim value per group; used as new dim coords)
   - dim_name::String           (which dimension was split on)
#
4. apply must NOT squeeze the split dimension — keep it as a size-1 axis so
   combine() can cat along it. Call f(group.data; dims = dim_idx) and rebuild
   the OutputVar with a 1-element dim array.
#
5. Time-based splitters (Season, Month, Year) require var.attributes["start_date"]
   to convert Float64 seconds -> DateTime via Utils.time_to_date. Error clearly
   if missing. ByDimension on a non-time dim has no such requirement.
#
6. Grouping logic for Season/Month/Year can mirror Utils.split_by_season_across_time
   and Utils.split_by_month. Utils.find_season_and_year is the key helper for Season.
   Convert DateTime groups back to seconds via Utils.date_to_time for indexing.
#
7. _split_along_dim(var, dim_name, split_vectors) in Var.jl is the core splitting
   primitive. It expects split_vectors::Vector{Vector{<eltype of dim>}}.
   Access as Var._split_along_dim since it is not exported.
#
8. combine() pattern (mirrors average_season_across_time in Var.jl ~line 2087):
   ret_data = cat(g.data..., dims = dim_idx)
   ret_dims[dim_name] = group_first_coords
   remake(first_group, data = ret_data, dims = ret_dims)
   Error clearly if groups have mismatched shapes on non-split dimensions.
#
9. ByDimension grouping: apply grouping_fn(element) -> label to each element of
   the chosen dim array; group by label in first-occurrence order; group_first_coords
   are the first original dim values (not the labels).
#
10. Imports needed:
    import ..Var: OutputVar, remake, _split_along_dim
    import ..Utils: time_to_date, date_to_time, find_season_and_year
    import Dates, OrderedCollections: OrderedDict
#
11. Regression anchor: var |> split(Season()) |> apply(nanmean) |> combine()
    should equal average_season_across_time(var).
#
12. apply kwargs passthrough: apply(f; kwargs...) splats into f(data; dims=..., kwargs...)
    so users don't need closures for e.g. var(data; dims, corrected=true).
#
GroupBy might be a better word here?
GroupBy internally assigns each index in the dimension a label and use that
for determing how to do things



end
