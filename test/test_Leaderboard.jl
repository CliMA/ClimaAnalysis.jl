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

@testset "Adding model and category" begin
    # Add single model
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "new model")
    @test ClimaAnalysis.model_names(rmse_var) ==
          ["ACCESS-CM2", "ACCESS-ESM1-5", "new model"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["DJF", "MAM", "JJA", "SON", "ANN"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("ACCESS-CM2" => "", "ACCESS-ESM1-5" => "", "new model" => "")
    @test all(isnan.(rmse_var["new model"]))

    # Add single category
    rmse_var = ClimaAnalysis.add_category(rmse_var, "new cat")
    @test ClimaAnalysis.model_names(rmse_var) ==
          ["ACCESS-CM2", "ACCESS-ESM1-5", "new model"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["DJF", "MAM", "JJA", "SON", "ANN", "new cat"]
    @test ClimaAnalysis.rmse_units(rmse_var) ==
          Dict("ACCESS-CM2" => "", "ACCESS-ESM1-5" => "", "new model" => "")
    @test all(isnan.(rmse_var[:, "new cat"]))

    # Add multiple models
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "new model", "new model 2")
    @test ClimaAnalysis.model_names(rmse_var) ==
          ["ACCESS-CM2", "ACCESS-ESM1-5", "new model", "new model 2"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["DJF", "MAM", "JJA", "SON", "ANN"]
    @test ClimaAnalysis.rmse_units(rmse_var) == Dict(
        "ACCESS-CM2" => "",
        "ACCESS-ESM1-5" => "",
        "new model" => "",
        "new model 2" => "",
    )
    @test all(isnan.(rmse_var["new model"]))
    @test all(isnan.(rmse_var["new model 2"]))

    # Add multiple categories
    rmse_var = ClimaAnalysis.add_category(rmse_var, "new cat", "new cat 2")
    @test ClimaAnalysis.model_names(rmse_var) ==
          ["ACCESS-CM2", "ACCESS-ESM1-5", "new model", "new model 2"]
    @test ClimaAnalysis.category_names(rmse_var) ==
          ["DJF", "MAM", "JJA", "SON", "ANN", "new cat", "new cat 2"]
    @test ClimaAnalysis.rmse_units(rmse_var) == Dict(
        "ACCESS-CM2" => "",
        "ACCESS-ESM1-5" => "",
        "new model" => "",
        "new model 2" => "",
    )
    @test all(isnan.(rmse_var[:, "new cat"]))
    @test all(isnan.(rmse_var[:, "new cat 2"]))
end

@testset "Removing model" begin
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var = ClimaAnalysis.add_model(rmse_var, "new model", "new model 2")
    rmse_var[2, 5] = 10.0
    rmse_var_delete = ClimaAnalysis.Leaderboard._delete_model(
        rmse_var,
        "ACCESS-CM2",
        "new model",
    )
    @test rmse_var_delete.short_name == "ta"
    @test ClimaAnalysis.model_names(rmse_var_delete) ==
          ["ACCESS-ESM1-5", "new model 2"]
    @test ClimaAnalysis.category_names(rmse_var_delete) ==
          ClimaAnalysis.category_names(rmse_var)
    @test rmse_var_delete[1, :] == rmse_var[2, :]
    @test all(isnan.(rmse_var_delete[2, :]))

    # Delete all models
    rmse_var_delete = ClimaAnalysis.Leaderboard._delete_model(
        rmse_var_delete,
        "ACCESS-ESM1-5",
        "new model 2",
    )
    @test rmse_var_delete.short_name == "ta"
    @test isempty(ClimaAnalysis.model_names(rmse_var_delete))
    @test ClimaAnalysis.category_names(rmse_var_delete) ==
          ClimaAnalysis.category_names(rmse_var)
    @test rmse_var_delete[:, :] |> size == (0, 5)
end

@testset "Adding units" begin
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")

    # Adding a single unit
    rmse_var = ClimaAnalysis.add_model(rmse_var, "hi")
    ClimaAnalysis.add_unit!(rmse_var, "hi", "test")
    @test ClimaAnalysis.rmse_units(rmse_var)["hi"] == "test"

    # Adding multiple units
    rmse_var = ClimaAnalysis.add_model(rmse_var, "hello1", "hello2")
    ClimaAnalysis.add_unit!(
        rmse_var,
        Dict("hello1" => "units1", "hello2" => "units2"),
    )
    @test ClimaAnalysis.rmse_units(rmse_var)["hello1"] == "units1"
    @test ClimaAnalysis.rmse_units(rmse_var)["hello2"] == "units2"
end

@testset "Finding best, worst, and median model" begin
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var[:, :] = [[1.0 2.0 3.0 4.0 5.0]; [6.0 7.0 8.0 9.0 10.0]]
    ClimaAnalysis.add_unit!(rmse_var, "ACCESS-CM2", "units")
    ClimaAnalysis.add_unit!(rmse_var, "ACCESS-ESM1-5", "units")
    val, model_name =
        ClimaAnalysis.find_best_single_model(rmse_var, category_name = "ANN")
    @test model_name == "ACCESS-CM2"
    @test val == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test val |> size == (5,)

    val, model_name =
        ClimaAnalysis.find_worst_single_model(rmse_var, category_name = "ANN")
    @test model_name == "ACCESS-ESM1-5"
    @test val == [6.0, 7.0, 8.0, 9.0, 10.0]
    @test val |> size == (5,)

    val = ClimaAnalysis.median(rmse_var)
    @test val == [7.0, 9.0, 11.0, 13.0, 15.0] ./ 2.0
    @test val |> size == (5,)

    # Test with NaN in RMSE array
    rmse_var = ClimaAnalysis.add_model(rmse_var, "for adding NaN")
    ClimaAnalysis.add_unit!(rmse_var, "for adding NaN", "units")
    val, model_name =
        ClimaAnalysis.find_best_single_model(rmse_var, category_name = "ANN")
    @test model_name == "ACCESS-CM2"
    @test val == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test val |> size == (5,)

    val, model_name =
        ClimaAnalysis.find_worst_single_model(rmse_var, category_name = "ANN")
    @test model_name == "ACCESS-ESM1-5"
    @test val == [6.0, 7.0, 8.0, 9.0, 10.0]
    @test val |> size == (5,)

    val = ClimaAnalysis.median(rmse_var)
    @test val == [7.0, 9.0, 11.0, 13.0, 15.0] ./ 2.0
    @test val |> size == (5,)
end

@testset "Reordering categories" begin
    csv_file_path = joinpath(@__DIR__, "sample_data/test_csv.csv")
    rmse_var = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var[:, :] = [[1.0 2.0 3.0 4.0 5.0]; [6.0 7.0 8.0 9.0 10.0]]

    # Identity reordering
    rmse_var_reordered = ClimaAnalysis.reorder_categories(
        rmse_var,
        ["DJF", "MAM", "JJA", "SON", "ANN"],
    )
    @test rmse_var_reordered.RMSEs == rmse_var.RMSEs
    @test ClimaAnalysis.category_names(rmse_var_reordered) ==
          ClimaAnalysis.category_names(rmse_var)
    @test ClimaAnalysis.model_names(rmse_var_reordered) ==
          ClimaAnalysis.model_names(rmse_var)

    # Shift everything by one
    rmse_var_reordered = ClimaAnalysis.reorder_categories(
        rmse_var,
        ["ANN", "DJF", "MAM", "JJA", "SON"],
    )
    @test rmse_var_reordered.RMSEs == rmse_var.RMSEs[:, [5, 1, 2, 3, 4]]
    @test ClimaAnalysis.category_names(rmse_var_reordered) ==
          ClimaAnalysis.category_names(rmse_var)[[5, 1, 2, 3, 4]]
    @test ClimaAnalysis.model_names(rmse_var_reordered) ==
          ClimaAnalysis.model_names(rmse_var)

    # Matching the order of two different RMSE vars
    rmse_var1 = ClimaAnalysis.read_rmses(csv_file_path, "ta")
    rmse_var1[:, :] = [[11.0 12.0 13.0 14.0 15.0]; [16.0 17.0 18.0 19.0 20.0]]

    # Order are the same
    rmse_var_reordered = ClimaAnalysis.match_category_order(rmse_var, rmse_var1)
    @test rmse_var_reordered.RMSEs == rmse_var.RMSEs
    @test ClimaAnalysis.category_names(rmse_var_reordered) ==
          ClimaAnalysis.category_names(rmse_var)
    @test ClimaAnalysis.model_names(rmse_var_reordered) ==
          ClimaAnalysis.model_names(rmse_var)

    # Order are different
    rmse_var1 = ClimaAnalysis.reorder_categories(
        rmse_var,
        ["ANN", "DJF", "MAM", "JJA", "SON"],
    )
    rmse_var_reordered = ClimaAnalysis.match_category_order(rmse_var, rmse_var1)
    @test rmse_var_reordered.RMSEs == rmse_var.RMSEs[:, [5, 1, 2, 3, 4]]
    @test ClimaAnalysis.category_names(rmse_var_reordered) ==
          ClimaAnalysis.category_names(rmse_var)[[5, 1, 2, 3, 4]]
    @test ClimaAnalysis.model_names(rmse_var_reordered) ==
          ClimaAnalysis.model_names(rmse_var)

    # Test for error handling
    @test_throws ErrorException ClimaAnalysis.reorder_categories(
        rmse_var,
        ["This", "should", "not", "work!"],
    )
    rmse_var_diff_cats = ClimaAnalysis.Leaderboard.RMSEVariable(
        rmse_var.short_name,
        ClimaAnalysis.model_names(rmse_var),
        ["!", "work", "not", "should", "This"],
        rmse_var.RMSEs,
        rmse_var.units,
    )
    @test_throws ErrorException ClimaAnalysis.match_category_order(
        rmse_var,
        rmse_var_diff_cats,
    )
end
