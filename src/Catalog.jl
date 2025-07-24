module Catalog

import NCDatasets

import Base: get
import ..Sim: available_vars
import ..Var: OutputVar

export NCCatalog, add_file!, available_vars

"""
    Organize external NetCDF files to easily initialize `OutputVar`s.
"""
struct NCCatalog
    "Dictionary that maps aliases to tuples of file path and variable name in the NetCDF
    file. The identifier can either be the variable name or an alias for the variable name
    in the file."
    identifier_to_path_varname::Dict{String, Tuple{String, String}}

    "Dictionary that maps variable names to filepaths. This is a fallback if a variable is
    requested and cannot be found in `identifier_to_path_varname`."
    varname_to_filepath::Dict{String, String}
end

"""
    NCCatalog()

Construct an empty `NCCatalog`.

Files can be added with [`ClimaAnalysis.Catalog.add_file!`](@ref).
"""
function NCCatalog()
    identifier_to_path_varname = Dict{String, Tuple{String, String}}()
    varname_to_filepath = Dict{String, String}()
    return NCCatalog(identifier_to_path_varname, varname_to_filepath)
end

"""
    add_file!(catalog::NCCatalog, filepath, short_names...)

Add the file at `filepath` to `catalog` with `short_names`.

The argument `short_names` can either be a short name (e.g., `"rsdt"`) or a pair of strings
(e.g., `"solar_mon" => "rsdt"`). If a pair is used, then the first name is the variable name
in the NetCDF file and the second name is the alias which will be used for `get`.

!!! note "Multiple files with the same short names"
    If there are multiple files that contain the same short names, then you should pass in
    the desired file first or pass in a short name or a pair of short name to alias for
    `short_names` to make the variable in the file retrievable with `get`.
"""
function add_file!(catalog::NCCatalog, filepath, short_names...)
    isfile(filepath) || throw(ArgumentError("$filepath is not a file"))
    endswith(filepath, ".nc") || throw(
        ArgumentError(
            "$filepath is not a NetCDF file, because the file path does not end with .nc",
        ),
    )

    maybe_short_name_pairs = short_names
    NCDatasets.NCDataset(filepath) do nc
        unordered_dims = NCDatasets.dimnames(nc)
        varnames = setdiff(keys(nc), unordered_dims)

        foreach(
            varname -> get!(catalog.varname_to_filepath, varname, filepath),
            varnames,
        )

        for maybe_short_name_pair in maybe_short_name_pairs
            maybe_short_name_pair isa Pair || (
                maybe_short_name_pair =
                    maybe_short_name_pair => maybe_short_name_pair
            )
            var_name = maybe_short_name_pair[1]
            identifier = maybe_short_name_pair[2]
            if !(var_name in varnames)
                @warn "Cannot find $var_name in file. Variable will not be added"
            elseif identifier in keys(catalog.identifier_to_path_varname)
                @warn "Identifier $identifier is already in catalog. Variable will not be added"
            else
                catalog.identifier_to_path_varname[identifier] =
                    (filepath, var_name)
            end
        end
    end
    return nothing
end

"""
    get(catalog::NCCatalog, short_name; var_kwargs)

Return a `OutputVar` with `short_name` from `catalog`.

The keyword argument `var_kwargs` is a tuple of keyword arguments passed to `OutputVar`. See
[`ClimaAnalysis.Var.OutputVar`](@ref) for a list of keyword arguments that can be passed.

Any short names or aliases specified when adding files are searched first to identify which
variable to load. If `short_name` cannot be found, then files are searched in the order that
the files are added to identify which variable to load.
"""
function get(catalog::NCCatalog, short_name; var_kwargs = ())
    filepath_varname_tuple =
        get(catalog.identifier_to_path_varname, short_name, nothing)
    if isnothing(filepath_varname_tuple)
        filepath = get(catalog.varname_to_filepath, short_name, nothing)
        isnothing(filepath) &&
            error("Cannot find variable with $short_name in catalog")
        varname = short_name
    else
        filepath, varname = filepath_varname_tuple
    end
    var = OutputVar(filepath, varname; var_kwargs...)
    var.attributes["short_name"] = short_name
    return var
end

"""
    available_vars(catalog::NCCatalog)

Return the short names of the variables found in `catalog`.

Both aliases and variable names in the NetCDF files are returned.
"""
function available_vars(catalog::NCCatalog)
    varnames_from_files = keys(catalog.varname_to_filepath)
    identifiers = keys(catalog.identifier_to_path_varname)
    return union(identifiers, varnames_from_files)
end

"""
    Base.isempty(catalog::NCCatalog)

Check if `NCCatalog` contains any files.
"""
function Base.isempty(catalog::NCCatalog)
    return isempty(catalog.varname_to_filepath)
end

end
