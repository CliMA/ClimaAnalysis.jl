module Sim

import Base: get

export SimDir,
    available_vars,
    available_reductions,
    available_periods,
    available_coord_types

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

    for (root, _, files) in walkdir(simulation_path)
        for file in files
            m = Utils.match_nc_filename(file)
            if !isnothing(m)
                short_name, period, reduction, coord_type = m

                full_path = joinpath(root, file)

                # Get the dictionary variable_paths["short_name"] if it exists, otherwise make it
                # a new dictionary.
                variable_reduction = get!(variable_paths, short_name, Dict())

                variable_reduction_period =
                    get!(variable_reduction, reduction, Dict())

                variable_reduction_period_coordtype =
                    get!(variable_reduction_period, period, Dict())

                # Store file paths as a vector of strings
                if haskey(variable_reduction_period_coordtype, coord_type)
                    # Do not sort because walkdir will give the files in top-down order
                    push!(
                        variable_reduction_period_coordtype[coord_type],
                        full_path,
                    )
                else
                    push!(
                        variable_reduction_period_coordtype,
                        coord_type => [full_path],
                    )
                end

                # Do the same for `vars`
                vars_reduction = get!(vars, short_name, Dict())
                vars_reduction_period = get!(vars_reduction, reduction, Dict())
                vars_reduction_period_coordtype =
                    get!(vars_reduction_period, period, Dict())
                push!(vars_reduction_period_coordtype, coord_type => nothing)

                # Add to allfiles
                push!(allfiles, full_path)

                # At the end we have four layers of dictionaries.
                #
                # Each layer corresponds to a component of the filename:
                # short name, reduction, period, and coordinate type.
                #
                # The first layer maps variable short names to a dictionary,
                # which maps reductions to a dictionary, which maps periods to a
                # dictionary, which maps coordinate types to a vector of the
                # path(s) of the file(s).
                #
                # Example: variable_paths = {"ta" : {"max": {"6.0h": {"pfull" =>
                # ["file.nc"]}}}}
                # Example: variable_paths = {"ta" : {"max": {"6.0h" => {nothing
                # => ["output_0001/file.nc", "output_0002/file.nc"]}}}}
            end
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

"""
    available_coord_types(
        simdir::SimDir;
        short_name::String,
        reduction::String,
        period::Union{String, Nothing}
    )

Return the available coordinate types associated with the given variable,
reduction, and period.

Note that `nothing` is reserved for diagnostics on their native space.

!!! note "Compatibility"
    This function is available in versions of ClimaAnalysis after v0.5.20.
"""
function available_coord_types(
    simdir::SimDir;
    short_name::String,
    reduction::String,
    period::Union{String, Nothing},
)
    if !(period in available_periods(simdir; short_name, reduction))
        error(
            "Period $period not available for $short_name. Available: $(available_periods(simdir; short_name, reduction))",
        )
    end

    return keys(simdir.vars[short_name][reduction][period]) |> Set
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
        period = nothing,
        coord_type = missing)

Return a `OutputVar` for the corresponding combination of `short_name`, `reduction`,
`period` (if it exists), and `coord_type` (if it exists).

The variable is read only once and saved into the `simdir`.

Keyword arguments
==================

When passing `nothing` to `reduction` or `period` or `missing` to `coord_type`,
`ClimaAnalysis` will try to automatically deduce the value. An error will be thrown if this
is not possible.

For instance, if the simulation has only one `ta`, you do not need to specify the
`reduction`, and `period` (`short_name` is enough). Similarly, if there is only one
`ta_average` (ie, not multiple periods), `short_name` and `reduction` will be enough.

!!! note "`coord_type` keyword argument"
    The `coord_type` keyword argument is available in versions of ClimaAnalysis after
    v0.5.20.
"""
function get(
    simdir::SimDir;
    short_name::String,
    reduction::Union{String, Nothing} = nothing,
    period::Union{String, Nothing} = nothing,
    coord_type::Union{String, Nothing, Missing} = missing,
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

    if ismissing(coord_type)
        coords = available_coord_types(simdir; short_name, reduction, period)
        length(coords) == 1 || error(
            "Found multiple coordinates for $short_name: $coords. You have to specify it.",
        )
        coord_type = pop!(coords)
    else
        if !(
            coord_type in
            available_coord_types(simdir; short_name, reduction, period)
        )
            error(
                "Coordinates $coord_type not available for $short_name, reduction $reduction, and period $period. " *
                "Available: $(available_coord_types(simdir; short_name, reduction, period))",
            )
        end
    end

    # Variable has not been read before. Read it now.
    if isnothing(simdir.vars[short_name][reduction][period][coord_type])
        file_path =
            simdir.variable_paths[short_name][reduction][period][coord_type]
        simdir.vars[short_name][reduction][period][coord_type] =
            read_var(file_path)
    end
    return simdir.vars[short_name][reduction][period][coord_type]
end

"""
    get(simdir::SimDir, short_names...)

If only one reduction, period, and coordinate type exist for the given `short_name`s, return
the corresponding `OutputVar`s.
"""
function get(simdir::SimDir, short_names...)
    results = map(
        short_name -> get(
            simdir;
            short_name,
            reduction = nothing,
            period = nothing,
            coord_type = missing,
        ),
        short_names,
    )
    return length(results) == 1 ? first(results) : results
end

"""
    isempty(simdir::SimDir)

Check if the given SimDir contains OutputVars.
"""
Base.isempty(simdir::SimDir) = isempty(simdir.vars)

end
