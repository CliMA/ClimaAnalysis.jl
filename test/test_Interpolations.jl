using Test
import ClimaAnalysis

@testset "Get indices and sign" begin
    indices = ((1, 2),)
    @test ClimaAnalysis.Interpolations.get_indices(indices, 0) == (2,)
    @test ClimaAnalysis.Interpolations.get_indices(indices, 1) == (1,)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 0) ==
          (1,)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 1) ==
          (2,)
    @test ClimaAnalysis.Interpolations.get_sign(indices, 0) == -1
    @test ClimaAnalysis.Interpolations.get_sign(indices, 1) == 1

    indices = ((1, 2), (3, 4))
    @test ClimaAnalysis.Interpolations.get_indices(indices, 0) == (2, 4)
    @test ClimaAnalysis.Interpolations.get_indices(indices, 1) == (1, 4)
    @test ClimaAnalysis.Interpolations.get_indices(indices, 2) == (2, 3)
    @test ClimaAnalysis.Interpolations.get_indices(indices, 3) == (1, 3)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 0) ==
          (1, 3)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 1) ==
          (2, 3)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 2) ==
          (1, 4)
    @test ClimaAnalysis.Interpolations.get_complement_indices(indices, 3) ==
          (2, 4)
    @test ClimaAnalysis.Interpolations.get_sign(indices, 0) == 1
    @test ClimaAnalysis.Interpolations.get_sign(indices, 1) == -1
    @test ClimaAnalysis.Interpolations.get_sign(indices, 2) == -1
    @test ClimaAnalysis.Interpolations.get_sign(indices, 3) == 1
end

@testset "Find indices for cell" begin
    val1 = 5
    ax1 = [0, 10]
    @test ClimaAnalysis.Interpolations.find_cell_indices_for_ax(val1, ax1) ==
          (1, 2)

    val2 = 6
    ax2 = [0, 4, 10]
    @test ClimaAnalysis.Interpolations.find_cell_indices_for_ax(val2, ax2) ==
          (2, 3)

    @test ClimaAnalysis.Interpolations.find_cell_indices_for_axes(
        (val1, val2),
        (ax1, ax2),
    ) == ((1, 2), (2, 3))
end

@testset "Extrapolation conditions" begin
    throw = ClimaAnalysis.Interpolations.extp_cond_throw()
    flat = ClimaAnalysis.Interpolations.extp_cond_flat()
    periodic = ClimaAnalysis.Interpolations.extp_cond_periodic()

    ax = [0, 1, 2, 3]
    @test_throws ErrorException throw.get_val_for_point(10, ax)
    @test_throws ErrorException throw.get_val_for_point(-1, ax)
    @test throw.get_val_for_point(1.5, ax) == 1.5

    @test flat.get_val_for_point(10, ax) == 3
    @test flat.get_val_for_point(-1, ax) == 0
    @test flat.get_val_for_point(1.5, ax) == 1.5

    @test periodic.get_val_for_point(10, ax) == 1
    @test periodic.get_val_for_point(-1, ax) == 2
    @test periodic.get_val_for_point(1.5, ax) == 1.5
    @test periodic.get_val_for_point(3, ax) == 3
end

@testset "Extrapolate to new point" begin
    throw = ClimaAnalysis.Interpolations.extp_cond_throw()
    flat = ClimaAnalysis.Interpolations.extp_cond_flat()
    periodic = ClimaAnalysis.Interpolations.extp_cond_periodic()

    ax1 = [0, 1, 2, 3]
    @test ClimaAnalysis.Interpolations.extp_to_point((1,), (ax1,), (throw,)) ==
          (1,)

    ax2 = [4, 5, 6, 7]
    @test ClimaAnalysis.Interpolations.extp_to_point(
        (-1, 8),
        (ax1, ax2),
        (flat, periodic),
    ) == (0, 5)
end

@testset "Interpolation" begin
    throw = ClimaAnalysis.Interpolations.extp_cond_throw()
    flat = ClimaAnalysis.Interpolations.extp_cond_flat()
    periodic = ClimaAnalysis.Interpolations.extp_cond_periodic()

    # 1D case
    axes = ([1.0, 2.0, 3.0],)
    data = [3.0, 1.0, 0.0]

    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.0,),
        axes,
        data,
        (throw,),
    ) == 3.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (3.0,),
        axes,
        data,
        (throw,),
    ) == 0.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5,),
        axes,
        data,
        (throw,),
    ) == 2.0

    # 1D case with extrapolation conditions
    @test_throws ErrorException ClimaAnalysis.Interpolations.linear_interpolate(
        (0.0,),
        axes,
        data,
        (throw,),
    )
    @test_throws ErrorException ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0,),
        axes,
        data,
        (throw,),
    )
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (0.0,),
        axes,
        data,
        (flat,),
    ) == 3.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0,),
        axes,
        data,
        (flat,),
    ) == 0.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (0.0,),
        axes,
        data,
        (periodic,),
    ) == 1.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0,),
        axes,
        data,
        (periodic,),
    ) == 1.0

    # 2D case
    axes = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
    data = reshape(1:9, (3, 3))

    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.0, 4.0),
        axes,
        data,
        (throw, throw),
    ) == 1.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (3.0, 6.0),
        axes,
        data,
        (throw, throw),
    ) == 9.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (2.0, 5.0),
        axes,
        data,
        (throw, throw),
    ) == 5.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5, 4.5),
        axes,
        data,
        (throw, throw),
    ) == 3.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5, 5.5),
        axes,
        data,
        (throw, throw),
    ) == 6.0

    # 2D cases with extrapolation conditions
    @test_throws ErrorException ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0, 5.0),
        axes,
        data,
        (throw, flat),
    )
    @test_throws ErrorException ClimaAnalysis.Interpolations.linear_interpolate(
        (2.0, 7.0),
        axes,
        data,
        (flat, throw),
    )
    @test_throws ErrorException ClimaAnalysis.Interpolations.linear_interpolate(
        (0.0, 8.0),
        axes,
        data,
        (throw, throw),
    )
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (0.0, 8.0),
        axes,
        data,
        (flat, flat),
    ) == 7.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0, 7.0),
        axes,
        data,
        (periodic, periodic),
    ) == 5.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (3.0, 6.0),
        axes,
        data,
        (periodic, periodic),
    ) == 9.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (4.0, 7.0),
        axes,
        data,
        (flat, periodic),
    ) == 6.0

    # 3D cases with extrapolation conditions
    axes = ([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0])
    data = reshape(1:27, (3, 3, 3))
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.0, 5.0, 7.0),
        axes,
        data,
        (throw, throw, throw),
    ) == 4.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5, 5.2, 7.5),
        axes,
        data,
        (throw, throw, throw),
    ) ≈ 9.6

    # Non equispaced
    axes = ([1.0, 3.0, 7.0], [4.0, 5.0, 7.0])
    data = reshape(1:9, (3, 3))
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (2.0, 4.5),
        axes,
        data,
        (throw, throw, throw),
    ) == 3.0
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (5.0, 6.0),
        axes,
        data,
        (throw, throw, throw),
    ) == 7.0

    # Different types
    # Axes have different types and inputs have different types
    axes = ([1.0f0, 2.0f0], [Float16(3.0), Float16(4.0)])
    data = [[1.0, 2.0] [3.0, 4.0]]
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5f0, 3.5),
        axes,
        data,
        (throw, throw),
    ) == 2.5
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        (1.5, 4.5f0),
        axes,
        data,
        (flat, flat),
    ) == 3.5
    @test ClimaAnalysis.Interpolations.linear_interpolate(
        [1.5, 4.5f0],
        axes,
        data,
        (flat, flat),
    ) == 3.5

    # Single number
    axes = ([1.0, 2.0, 3.0],)
    data = [3.0, 1.0, 0.0]

    @test ClimaAnalysis.Interpolations.linear_interpolate(
        1.0,
        axes,
        data,
        (throw,),
    ) == 3.0
end
