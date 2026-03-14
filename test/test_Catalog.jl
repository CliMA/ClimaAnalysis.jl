using Test
import ClimaAnalysis
import Dates

@testset "NCCatalog" begin
    catalog = ClimaAnalysis.NCCatalog()
    @test isempty(catalog)

    gpp_path = joinpath(@__DIR__, "sample_nc/test_gpp.nc")
    ClimaAnalysis.add_file!(catalog, gpp_path, "gpp" => "cool_gpp")
    @test !isempty(catalog)
    var_names = ["cool_gpp", "gpp"]
    available_vars_in_catalog = ClimaAnalysis.available_vars(catalog)
    for var_name in var_names
        @test var_name in available_vars_in_catalog
    end

    pr_path = joinpath(@__DIR__, "sample_nc/test_pr.nc")
    ClimaAnalysis.add_file!(catalog, pr_path, "precip" => "pr")
    @test !isempty(catalog)
    var_names = ["pr", "precip", "cool_gpp", "gpp"]
    available_vars_in_catalog = ClimaAnalysis.available_vars(catalog)
    for var_name in var_names
        @test var_name in available_vars_in_catalog
    end

    pr_path2 = joinpath(@__DIR__, "sample_nc/test_pr2.nc")
    ClimaAnalysis.add_file!(catalog, pr_path2)
    @test !isempty(catalog)
    var_names = ["pr", "precip", "cool_gpp", "gpp"]
    available_vars_in_catalog = ClimaAnalysis.available_vars(catalog)
    for var_name in var_names
        @test var_name in available_vars_in_catalog
    end

    # Test errors
    not_nc_path = joinpath(@__DIR__, "test_Catalog.jl")
    @test_throws ArgumentError ClimaAnalysis.add_file!(catalog, not_nc_path)
    @test_throws ArgumentError ClimaAnalysis.add_file!(catalog, @__DIR__)

    # Test warnings
    @test_logs (:warn, r"Cannot find") ClimaAnalysis.add_file!(
        catalog,
        gpp_path,
        "rad_gpp",
    )
    @test_logs (:warn, r"Cannot find") ClimaAnalysis.add_file!(
        catalog,
        gpp_path,
        "rad_gpp" => "more_rad_gpp",
    )
    @test_logs (:warn, r"is already in catalog") ClimaAnalysis.add_file!(
        catalog,
        gpp_path,
        "gpp" => "cool_gpp",
    )

    gpp = get(catalog, "gpp")
    @test gpp.attributes["short_name"] == "gpp"

    # Both OutputVars should be the same since pr_path was added before pr_path2
    precip = get(catalog, "precip")
    pr = get(catalog, "pr")
    @test precip.data == pr.data
    @test precip.dims == precip.dims
    @test all(
        precip.attributes[key] == pr.attributes[key] for
        key in union(keys(precip.attributes), keys(pr.attributes)) if
        key != "short_name"
    )
    @test precip.dim_attributes == pr.dim_attributes
    @test precip.attributes["short_name"] == "precip"

    # Test keyword arguments work
    pr = get(catalog, "pr", var_kwargs = (shift_by = Dates.lastdayofmonth,))
    @test ClimaAnalysis.dates(pr) ==
          [Dates.DateTime(1979, 1, 31), Dates.DateTime(2023, 5, 31)]
end

@testset "Add multiple files to concatenate for NCCatalog" begin
    catalog = ClimaAnalysis.NCCatalog()
    ncpath1 = joinpath(@__DIR__, "sample_nc/test_pr.nc")
    ncpath2 = joinpath(@__DIR__, "sample_nc/test_pr2.nc")

    ClimaAnalysis.add_file!(
        catalog,
        [ncpath1, ncpath2],
        "precip" => "all_precip",
    )
    ClimaAnalysis.add_file!(catalog, [ncpath1], "precip" => "pr1")
    ClimaAnalysis.add_file!(catalog, [ncpath2], "precip" => "pr2")
    pr_var_all = get(catalog, "all_precip")
    pr_var1 = get(catalog, "pr1")
    pr_var2 = get(catalog, "pr2")

    # NCCatalog use the alias (e.g. pr1 and pr2) for the short name so we update
    # it to match
    pr_var1.attributes["short_name"] = "all_precip"

    @test pr_var_all.attributes == pr_var1.attributes
    @test ClimaAnalysis.dates(pr_var_all) ==
          vcat(ClimaAnalysis.dates(pr_var1), ClimaAnalysis.dates(pr_var2))
    @test pr_var_all.dims["lon"] == pr_var1.dims["lon"]
    @test pr_var_all.dims["lat"] == pr_var1.dims["lat"]
    @test pr_var_all.data == cat(pr_var1.data, pr_var2.data, dims = 3)
    @test pr_var_all.dim_attributes == pr_var1.dim_attributes

    # Times are not sorted
    @test_throws ErrorException ClimaAnalysis.add_file!(
        catalog,
        [ncpath2, ncpath1],
    )
end
