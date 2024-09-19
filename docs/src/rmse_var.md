# `RMSEVariable`s

`RMSEVariable`s contain all the information needed to process and compare root mean squared
errors (RMSEs) between different models and categories (e.g., seasons) for a single variable
of interest.

`ClimaAnalysis` provides several constructors for making a `RMSEVariable`. For all
constructors, a short name and a vector of model names must be provided. If units are not
provided, then each model will have no unit which denotes the missing unit. See the examples
below where the constructor can take in a short name, a vector of model names, a vector of
categories, and a dictionary mapping model names to units or a string of the name of the
unit.

```@example rmse_var
import ClimaAnalysis

rmse_var = ClimaAnalysis.RMSEVariable("ta", ["ACCESS-CM2", "ACCESS-ESM1-5"])
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"),
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"),
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    ones(2, 5),
    Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"),
)
# Convenience functions if models all share the same unit
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    "K",
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    "K",
)
rmse_var = ClimaAnalysis.RMSEVariable(
    "ta",
    ["ACCESS-CM2", "ACCESS-ESM1-5"],
    ["DJF", "MAM", "JJA", "SON", "ANN"],
    ones(2, 5),
    "K",
)

nothing # hide
```

The `RMSEVariable` can be inspected using `model_names`, `category_names`, and `rmse_units`
which provide the model names, the category names, and the units respectively.

```@repl rmse_var
ClimaAnalysis.model_names(rmse_var)
ClimaAnalysis.category_names(rmse_var)
ClimaAnalysis.rmse_units(rmse_var)
```

## Reading RMSEs from CSV file

Typically, the root mean squared errors (RMSEs) of different models across different
categories are stored in a different file and need to be loaded in. `ClimaAnalysis` can load
this information from a CSV file and store it in a `RMSEVariable`. The format of the CSV
file should have a header consisting of the entry "model_name" (or any other text as it is
ignored by the function) and rest of the entries should be the category names. Each row
after the header should start with the model name and the root mean squared errors for each
category for that model. The entries of the CSV file should be separated by commas.

See the example below using `read_rmses` where data is loaded from `test_csv.csv` and a
short name of `ta` is provided. One can also pass in a dictionary mapping model names to
units for `units` or a string if the units are the same for all the models.

```@example rmse_var
rmse_var = ClimaAnalysis.read_rmses("./data/test_csv.csv", "ta")
rmse_var = ClimaAnalysis.read_rmses(
    "./data/test_csv.csv",
    "ta",
    units = Dict("ACCESS-CM2" => "K", "ACCESS-ESM1-5" => "K"), # passing units as a dictionary
)
rmse_var = ClimaAnalysis.read_rmses(
    "./data/test_csv.csv",
    "ta",
    units = "K", # passing units as a string
)

nothing # hide
```

## Indexing

After loading the data, one may want to inspect, change, or manipulate the data. This is
possible by the indexing functionality that `RMSEVariable` provides. Indexing into a
`RMSEVariable` is similar, but not the same as indexing into an array. Indexing by
integer or string is supported, but linear indexing (e.g. `rmse_var[1]`) is not supported.
integer or string is supported, but linear indexing (e.g., `rmse_var[1]`) is not supported.

```@repl rmse_var
rmse_var[:, :]
rmse_var["ACCESS-CM2"]
rmse_var[:, "MAM"]
rmse_var["ACCESS-CM2", ["ANN", "DJF", "MAM"]]
rmse_var[2,5] = 11.2;
rmse_var[:, :]
```

## Adding categories, models, and units

It may be the case that the CSV file does not contain all the models you want to analyze, or
you want to consider another category but do not want to go in and manually edit the CSV
file to add it. `ClimaAnalysis` provides `add_category`, `add_model`, and `add_unit!` for
adding categories, models, and units respectively. Multiple model or categories can be
provided (e.g., `add_model(rmse_var, "model1", "model2")`) in the functions. For adding
multiple units, one can pass in a dictionary mapping model names to units. See the example
below using this functionality.

```@julia rmse_var
rmse_var2 = ClimaAnalysis.add_category(rmse_var, "Jan") # can take in more than one category
rmse_var = ClimaAnalysis.add_model(rmse_var, "CliMA") # can take in more than one model name
ClimaAnalysis.add_unit!(rmse_var, "CliMA", "K")
ClimaAnalysis.add_unit!(rmse_var, Dict("CliMA" => "K")) # for adding multiple units
```

```@repl rmse_var
ClimaAnalysis.category_names(rmse_var2)
ClimaAnalysis.model_names(rmse_var)
ClimaAnalysis.rmse_units(rmse_var)
rmse_var[:,:]
```