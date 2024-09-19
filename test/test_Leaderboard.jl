using Test
import ClimaAnalysis

@testset "Constructors and helper functions" begin
    # Testing constructor using short name and model names
    rmse_var = ClimaAnalysis.RMSEVariable("ta", ["model1", "model2"])
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["ANN", "DJF", "MAM", "JJA", "SON"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "", "model2" => "")
    @test rmse_var.short_name == "ta"
    @test isnan(rmse_var.RMSEs[1, 1])

    #  Testing constructor using short name, model names, and provided units
    rmse_var = ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        Dict("model1" => "units1", "model2" => "units2"),
    )
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["ANN", "DJF", "MAM", "JJA", "SON"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "units1", "model2" => "units2")
    @test rmse_var.short_name == "ta"
    @test isnan(rmse_var.RMSEs[1, 1])

    #  Testing constructor using short name, model names, category_names, and provided units
    rmse_var = ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        ["cat1", "cat2"],
        Dict("model1" => "units1", "model2" => "units2"),
    )
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2"]
    @test ClimaAnalysis.category_names(rmse_var) == ["cat1", "cat2"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "units1", "model2" => "units2")
    @test rmse_var.short_name == "ta"
    @test isnan(rmse_var.RMSEs[1, 1])

    # Testing constructor using short name, model names, category names, RMSEs, and units
    rmse_var = ClimaAnalysis.RMSEVariable(
        "short_name",
        ["model1", "model2", "model3"],
        ["cat1", "cat2"],
        ones(3, 2),
        Dict("model1" => "units"),
    )
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2", "model3"]
    @test ClimaAnalysis.category_names(rmse_var) == ["cat1", "cat2"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "units", "model2" => "", "model3" => "")
    @test rmse_var.short_name == "short_name"
    @test rmse_var.RMSEs == ones(3, 2)

    # Testing constructors for passing units as a string
    rmse_var =
        ClimaAnalysis.RMSEVariable("ta", ["model1", "model2"], "same_unit")
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["ANN", "DJF", "MAM", "JJA", "SON"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "same_unit", "model2" => "same_unit")
    @test rmse_var.short_name == "ta"
    @test isnan(rmse_var.RMSEs[1, 1])

    rmse_var = ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        ["cat1", "cat2"],
        "same_unit",
    )
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2"]
    @test ClimaAnalysis.category_names(rmse_var) == ["cat1", "cat2"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "same_unit", "model2" => "same_unit")
    @test rmse_var.short_name == "ta"
    @test isnan(rmse_var.RMSEs[1, 1])

    rmse_var = ClimaAnalysis.RMSEVariable(
        "short_name",
        ["model1", "model2", "model3"],
        ["cat1", "cat2"],
        ones(3, 2),
        "units",
    )
    @test ClimaAnalysis.model_names(rmse_var) == ["model1", "model2", "model3"]
    @test ClimaAnalysis.category_names(rmse_var) == ["cat1", "cat2"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("model1" => "units", "model2" => "units", "model3" => "units")
    @test rmse_var.short_name == "short_name"
    @test rmse_var.RMSEs == ones(3, 2)


    # Error handling
    # Duplicate model names and category names
    @test_throws ErrorException ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model1"],
        ["cat1"],
        ones(2, 1),
        Dict("model1" => ""),
    )
    @test_throws ErrorException ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        ["cat1", "cat1"],
        ones(2, 2),
        Dict("model1" => ""),
    )
    # Negative RMSE
    @test_throws ErrorException ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        ["cat1", "cat1"],
        -1.0 .* ones(2, 2),
        Dict("model1" => ""),
    )
    # Shape of RMSE array
    @test_throws ErrorException ClimaAnalysis.RMSEVariable(
        "ta",
        ["model1", "model2"],
        ["cat1", "cat1"],
        ones(1, 2),
        Dict("model1" => ""),
    )
end
