"""
    SimDir(simulation_path::String)

Object that describes all the `ClimaAtmos` output found in the given `simulation_path`.
"""
struct SimDir{DV <: Dict}

    "Path where the output data is stored"
    simulation_path::String

    "Dictionary of dictionaries that maps variables/reductions/periods to files"
    variables_paths::DV

    "List of files that ClimaAtmos knows how to process"
    allfiles::Set{String}
end

function SimDir(simulation_path::String)

    variables_paths = Dict()
    allfiles = Set{String}()

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

    foreach(readdir(simulation_path)) do path
        m = match(re, path)
        if !isnothing(m)
            short_name, period, reduction = m.captures

            full_path = joinpath(simulation_path, path)

            # Get the dictionary variables_paths["short_name"] if it exists, otherwise make it
            # a new dictionary.
            variable_reduction = get!(variables_paths, short_name, Dict())

            if isnothing(period)
                push!(variable_reduction, reduction => full_path)
            else
                # We add another layer when we also have the period
                variable_reduction_period =
                    get!(variable_reduction, reduction, Dict())
                push!(variable_reduction_period, period => full_path)
            end

            # Add to allfiles
            push!(allfiles, full_path)

            # At the end we have three layers of dictionaries.
            #
            # The first layer maps variable short names to a dictionary, which maps
            # reductions to a dictionary, which maps periods to the path of the file that
            # contains that chain
        end
    end
    return SimDir(simulation_path, variables_paths, allfiles)
end

"""
    available_vars(simdir::SimDir)

Return the short names of the variables found in the given `simdir`.
"""
available_vars(simdir::SimDir) = keys(simdir.variables_paths) |> Set
