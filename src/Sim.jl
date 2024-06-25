module Sim

import Base: get

export SimDir, available_vars, available_reductions, available_periods

import ..Utils
import ..Var: read_var

"""
    SimDir(simulation_path::String)

Object that describes all the `ClimaAtmos` output found in the given `simulation_path`.
"""
struct SimDir{DV <: Dict, DOV <: Dict}

    "Path where the output data is stored"
    simulation_path::String

    "Dictionary of dictionaries that maps variables/reductions/periods to files"
    variable_paths::DV

    "Dictionary of dictionaries that maps variables/reductions/periods to OutputVars"
    vars::DOV

    "List of files that ClimaAtmos knows how to process"
    allfiles::Set{String}
end

function SimDir(simulation_path::String)

    variable_paths = Dict()
    vars = Dict()
    allfiles = Set{String}()

    foreach(readdir(simulation_path)) do path
        m = Utils.match_nc_filename(path)
        if !isnothing(m)
            short_name, period, reduction = m

            full_path = joinpath(simulation_path, path)

            # Get the dictionary variable_paths["short_name"] if it exists, otherwise make it
            # a new dictionary.
            variable_reduction = get!(variable_paths, short_name, Dict())

            variable_reduction_period =
                get!(variable_reduction, reduction, Dict())
            push!(variable_reduction_period, period => full_path)

            # Do the same for `vars`
            vars_reduction = get!(vars, short_name, Dict())
            vars_reduction_period = get!(vars_reduction, reduction, Dict())
            push!(vars_reduction_period, period => nothing)

            # Add to allfiles
            push!(allfiles, full_path)

            # At the end we have three layers of dictionaries.
            #
            # The first layer maps variable short names to a dictionary, which maps
            # reductions to a dictionary, which maps periods to the path of the file that
            # contains that chain
            #
            # Example: variable_paths = {"ta" : {"max": {"6.0h" => file.nc}}}
        end
    end
    return SimDir(simulation_path, variable_paths, vars, allfiles)
end

"""
    available_vars(simdir::SimDir)

Return the short names of the variables found in the given `simdir`.
"""
available_vars(simdir::SimDir) = keys(simdir.vars) |> Set

"""
    available_reductions(simdir::SimDir, short_name::String)

Return the reductions available for the given variable in the given `simdir`.
"""
function available_reductions(simdir::SimDir; short_name::String)
    if !(short_name in available_vars(simdir))
        error(
            "Variable $short_name not found. Available: $(available_vars(simdir))",
        )
    end

    return keys(simdir.vars[short_name]) |> Set
end

"""
    available_periods(simdir::SimDir, short_name::String, reduction::String)

Return the periods associated to the given variable and reduction.
"""
function available_periods(
    simdir::SimDir;
    short_name::String,
    reduction::String,
)
    if !(reduction in available_reductions(simdir; short_name))
        error(
            "Reduction $reduction not available for $short_name. Available: $(available_reductions(simdir; short_name))",
        )
    end

    return keys(simdir.vars[short_name][reduction]) |> Set
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
function get(
    simdir::SimDir;
    short_name::String,
    reduction::Union{String, Nothing} = nothing,
    period::Union{String, Nothing} = nothing,
)
    if isnothing(reduction)
        reductions = available_reductions(simdir; short_name)
        length(reductions) == 1 || error(
            "Found multiple reductions for $short_name: $reductions. You have to specify it.",
        )
        reduction = pop!(reductions)
    end

    if isnothing(period)
        periods = available_periods(simdir; short_name, reduction)
        length(periods) == 1 || error(
            "Found multiple periods for $short_name: $periods. You have to specify it.",
        )
        period = pop!(periods)
    else
        if !(period in available_periods(simdir; short_name, reduction))
            error(
                "Period $period not available for $short_name and reduction $reduction. " *
                "Available: $(available_periods(simdir; short_name, reduction))",
            )
        end
    end

    # Variable has not been read before. Read it now.
    if isnothing(simdir.vars[short_name][reduction][period])
        file_path = simdir.variable_paths[short_name][reduction][period]
        simdir.vars[short_name][reduction][period] = read_var(file_path)
    end
    return simdir.vars[short_name][reduction][period]
end

"""
    get(simdir::SimDir, short_name)

If only one reduction and period exist for `short_name`, return the corresponding
`OutputVar`.
"""
function get(simdir::SimDir, short_name::String)
    return get(simdir; short_name, reduction = nothing, period = nothing)
end

"""
    isempty(simdir::SimDir)

Check if the given SimDir contains OutputVars.
"""
Base.isempty(simdir::SimDir) = isempty(simdir.vars)

end
