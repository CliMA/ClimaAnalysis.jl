module Sim

import Base: get

export SimDir, available_vars, available_reductions, available_periods

import ..Utils
import ..Var: read_var

import OrderedCollections: OrderedDict

"""
    SimDir(simulation_path::String)

Object that describes all the `ClimaAtmos` output found in the given `simulation_path`.
"""
struct SimDir{DV, DOV}

    "Path where the output data is stored"
    simulation_path::String

    indices::DV

    component_to_indices::DOV

    filepaths::Vector{Vector{String}}
end

# TODO: Add a cache system to this

function SimDir(simulation_path::String)
    # TODO: Nothing case is annoying to handle
    indices = OrderedDict{NTuple{3, Union{String, Nothing}}, Vector{String}}()
    component_to_indices =
        Dict{String, Dict{Union{String, Nothing}, Set{Int64}}}(
            component => Dict() for
            component in ["short_name", "period", "reduction"]
        )

    for (root, _, files) in walkdir(simulation_path)
        for file in files
            m = Utils.match_nc_filename(file)
            if !isnothing(m)
                short_name, period, reduction = m

                full_path = joinpath(root, file)
                filepaths = get!(indices, (short_name, period, reduction), [])
                push!(filepaths, full_path)
            end
        end
    end

    filepaths = Vector{String}[]
    for (i, ((short_name, period, reduction), filepath)) in enumerate(indices)
        push!(get!(component_to_indices["short_name"], short_name, Set()), i)
        push!(get!(component_to_indices["reduction"], reduction, Set()), i)
        push!(get!(component_to_indices["period"], period, Set()), i)
        push!(filepaths, filepath)
    end
    return SimDir(simulation_path, indices, component_to_indices, filepaths)
end

"""
    available_vars(simdir::SimDir)

Return the short names of the variables found in the given `simdir`.
"""
available_vars(simdir::SimDir) =
    Set{String}(keys(simdir.component_to_indices["short_name"]))

"""
    available_reductions(simdir::SimDir, short_name::String)

Return the reductions available for the given variable in the given `simdir`.
"""
function available_reductions(simdir::SimDir; short_name::String)
    (; filepaths, component_to_indices) = simdir
    filenames =
        filepaths[collect(component_to_indices["short_name"][short_name])] .|>
        first .|>
        basename
    components_of_files = Utils.match_nc_filename.(filenames)
    # TODO: reduction can be empty
    reductions = Set{String}((
        reduction for
        (_, _, reduction) in components_of_files if !isnothing(reduction)
    ))
    return reductions
end

"""
    available_periods(simdir::SimDir, short_name::String, reduction::String)

Return the periods associated to the given variable and reduction.
"""
function available_periods(
    simdir::SimDir;
    short_name::String,
    reduction::String, # TODO: There is nothing special about short_name and reduction, so we privledge those?
)
    (; filepaths, component_to_indices) = simdir
    filenames_with_short_name = component_to_indices["short_name"][short_name]
    filenames_with_reduction = component_to_indices["reduction"][reduction]
    indices =
        collect(intersect(filenames_with_short_name, filenames_with_reduction))
    components_of_files =
        filepaths[indices] .|> first .|> basename .|> Utils.match_nc_filename
    periods = Set((period for (_, period, _) in components_of_files))
    return periods
end

function Base.summary(io::IO, simdir::SimDir)
    print(io, "Output directory: $(simdir.simulation_path)\n")
    print(io, "Variables:")
    for short_name in available_vars(simdir)
        print(io, "\n- $short_name")
        for reduction in available_reductions(simdir; short_name)
            print(io, "\n    $reduction")
            periods = available_periods(simdir; short_name, reduction)
            print(io, " (", join(periods, ", "), ")")
        end
    end
end

"""
    Base.show(io::IO, simdir::SimDir)

Pretty print the contents of `simdir`.

Print the output directory and the periods associated to the given variable and reduction.
"""
Base.show(io::IO, simdir::SimDir) = summary(io, simdir)

"""
    get(simdir::SimDir;
        short_name,
        reduction = nothing,
        period = nothing)

Return a `OutputVar` for the corresponding combination of `short_name`, `reduction`,
and `period` (if it exists).

The variable is read only once and saved into the `simdir`.

Keyword arguments
==================

When passing `nothing` to `reduction` and `period`, `ClimaAnalysis` will try to
automatically deduce the value. An error will be thrown if this is not possible.

For instance, if the simulation has only one `ta`, you do not need to specify `short_name`,
`reduction`, and `period` (`short_name` is enough). Similarly, if there is only one
`ta_average` (ie, not multiple periods), `short_name` and `reduction` will be enough.
"""
function get(simdir::SimDir; kwargs...)
    # TODO: Remove anything with nothing to maintain backward compatibility
    kwargs = Dict(String(k) => v for (k, v) in kwargs if !isnothing(v))
    (; filepaths, component_to_indices) = simdir
    matching_indices = reduce(
        intersect,
        (component_to_indices[key][val] for (key, val) in kwargs);
        init = Set(eachindex(filepaths)),
    )
    if length(matching_indices) > 1
        error("Not enough specialization")
    end
    return read_var(first(filepaths[collect(matching_indices)]))
end

"""
    get(simdir::SimDir, short_names...)

If only one reduction and period exist for the given `short_name`s, return the corresponding
`OutputVar`s.
"""
function get(simdir::SimDir, short_names...)
    results =
        map(short_name -> get(simdir, short_name = short_name), short_names)
    return length(results) == 1 ? first(results) : results
end

"""
    isempty(simdir::SimDir)

Check if the given SimDir contains OutputVars.
"""
Base.isempty(simdir::SimDir) = isempty(simdir.indices)

end
