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

@testset "Merge NCCatalog" begin
    gpp_path = joinpath(@__DIR__, "sample_nc/test_gpp.nc")
    pr_path = joinpath(@__DIR__, "sample_nc/test_pr.nc")
    pr_path2 = joinpath(@__DIR__, "sample_nc/test_pr2.nc")

    # Disjoint catalogs: the merged catalog contains the variables of both
    catalog_gpp = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_gpp, gpp_path, "gpp" => "cool_gpp")
    catalog_pr = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_pr, pr_path, "precip" => "pr")

    gpp_vars_before = ClimaAnalysis.available_vars(catalog_gpp)
    pr_vars_before = ClimaAnalysis.available_vars(catalog_pr)

    merged = merge(catalog_gpp, catalog_pr)
    merged_vars = ClimaAnalysis.available_vars(merged)
    for var_name in ("cool_gpp", "gpp", "pr", "precip")
        @test var_name in merged_vars
    end

    # merge does not modify the input catalogs
    @test ClimaAnalysis.available_vars(catalog_gpp) == gpp_vars_before
    @test ClimaAnalysis.available_vars(catalog_pr) == pr_vars_before
    @test !("pr" in ClimaAnalysis.available_vars(catalog_gpp))
    @test !("cool_gpp" in ClimaAnalysis.available_vars(catalog_pr))

    # Rightmost wins on a duplicate identifier. test_pr.nc and test_pr2.nc both
    # contain "precip" but cover different times, so the loaded data differs.
    catalog_a = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_a, pr_path, "precip" => "pr")
    catalog_b = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_b, pr_path2, "precip" => "pr")

    merged_ab = merge(catalog_a, catalog_b)
    @test get(merged_ab, "pr") == get(catalog_b, "pr")
    @test get(merged_ab, "pr") != get(catalog_a, "pr")

    # Rightmost wins on the varname_to_filepath fallback (no aliases)
    catalog_a_bare = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_a_bare, pr_path)
    catalog_b_bare = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_b_bare, pr_path2)

    merged_bare = merge(catalog_a_bare, catalog_b_bare)
    @test get(merged_bare, "precip") == get(catalog_b_bare, "precip")

    # Merging a single catalog yields an equivalent, independent catalog
    single = merge(catalog_a)
    @test ClimaAnalysis.available_vars(single) ==
          ClimaAnalysis.available_vars(catalog_a)

    # Merging three catalogs: rightmost precedence holds across all of them
    catalog_c = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_c, pr_path2, "precip" => "pr")
    merged_abc = merge(catalog_gpp, catalog_a, catalog_c)
    @test "cool_gpp" in ClimaAnalysis.available_vars(merged_abc)
    @test get(merged_abc, "pr") == get(catalog_c, "pr")

    # merge! mutates the first catalog in place and returns it
    catalog_dst = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_dst, pr_path, "precip" => "pr")
    catalog_src = ClimaAnalysis.NCCatalog()
    ClimaAnalysis.add_file!(catalog_src, gpp_path, "gpp" => "cool_gpp")
    ClimaAnalysis.add_file!(catalog_src, pr_path2, "precip" => "pr")

    result = merge!(catalog_dst, catalog_src)
    @test result === catalog_dst
    @test "cool_gpp" in ClimaAnalysis.available_vars(catalog_dst)
    # The duplicate "pr" identifier now resolves to catalog_src's file
    @test get(catalog_dst, "pr") == get(catalog_src, "pr")
end
