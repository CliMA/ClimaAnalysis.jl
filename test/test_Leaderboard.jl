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

@testset "Reading RMSEs from CSV file" begin
    # Testing constructor using CSV file
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    @test ClimaAnalysis.model_names(rmse_var) == ["ACCESS-CM2", "ACCESS-ESM1-5"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["DJF", "MAM", "JJA", "SON", "ANN"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("ACCESS-CM2" => "", "ACCESS-ESM1-5" => "")
    @test rmse_var.short_name == "ta"
    @test rmse_var.RMSEs[1, 1] == 11.941
    @test isnan(rmse_var.RMSEs[2, 5])

    # Testing constructor using CSV file with units provided
    rmse_var = ClimaAnalysis.read_rmses(
        csv_file_path,
        "ta",
        units = Dict("ACCESS-ESM1-5" => "m", "wacky" => "weird"),
    )
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("ACCESS-CM2" => "", "ACCESS-ESM1-5" => "m")

    # Testing constructor using CSV file with units being a string
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta", units = "m")
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("ACCESS-CM2" => "m", "ACCESS-ESM1-5" => "m")
end

@testset "Indexing" begin
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")

    rmse_var[1, 1] = 100.0
    @test rmse_var.RMSEs[1, 1] == 100.0

    rmse_var[1, 1:3] = [100.0, 110.0, 120.0]
    @test rmse_var.RMSEs[1, 1:3] == [100.0, 110.0, 120.0]

    rmse_var[1:2, 1] = [100.0, 110.0]'
    @test rmse_var.RMSEs[1:2, 1] == [100.0, 110.0]

    rmse_var["ACCESS-ESM1-5", "DJF"] = 200.0
    @test rmse_var["ACCESS-ESM1-5", "DJF"] == 200.0

    rmse_var["ACCESS-ESM1-5", ["DJF", "MAM", "ANN"]] = [200.0, 210.0, 220.0]
    @test rmse_var["ACCESS-ESM1-5", [1, 2, 5]] == [200.0, 210.0, 220.0]

    rmse_var["ACCESS-ESM1-5"] = [120.0, 130.0, 140.0, 150.0, 160.0]
    @test rmse_var["ACCESS-ESM1-5"] == [120.0, 130.0, 140.0, 150.0, 160.0]

    # Check error handling
    @test_throws ErrorException rmse_var[5, 5] = 100.0
    @test_throws ErrorException rmse_var["do not exist"] = 100.0
    @test_throws ErrorException rmse_var["do not exist", "test"] = 100.0
    @test_throws ErrorException rmse_var["ACCESS-ESM1-5", "test"] = 100.0
    @test_throws ErrorException rmse_var["model1", "ANN"] = 100.0
    @test_throws DimensionMismatch rmse_var["ACCESS-ESM1-5"] =
        [120.0, 130.0, 140.0]
    @test_throws DimensionMismatch rmse_var[1, :] = [120.0, 130.0, 140.0]
end
