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