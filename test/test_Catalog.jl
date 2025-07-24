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
