using Test
using BoundaryIntegral
import BoundaryIntegral as BI

include("../src/ScreenedOrbitalSolve.jl")
using .ScreenedOrbitalSolve

@testset "pair targets apply the graphene a1 shift" begin
    xs = [0.0, 0.5]
    ys = [0.0]
    zs = [0.0]
    weights = ones(length(xs), length(ys), length(zs))
    density_1 = reshape([1.0, 2.0], length(xs), length(ys), length(zs))
    density_2 = reshape([3.0, 4.0], length(xs), length(ys), length(zs))

    vs1 = BI.VolumeSource((xs, ys, zs), weights, density_1)
    vs2 = BI.VolumeSource((xs, ys, zs), weights, density_2)
    pairs = pair_targets(vs1, vs2; a1 = (2.465, 0.0, 0.0))

    @test keys(pairs) == (:U_00, :U_01, :U_02, :U_03)
    @test pairs.U_00.positions ≈ vs1.positions atol = 1e-12
    @test pairs.U_01.positions ≈ vs2.positions atol = 1e-12
    @test pairs.U_02.positions ≈ (vs1.positions .+ reshape([2.465, 0.0, 0.0], 3, 1)) atol = 1e-12
    @test pairs.U_03.positions ≈ (vs2.positions .+ reshape([2.465, 0.0, 0.0], 3, 1)) atol = 1e-12
end

@testset "target potential integration uses orbital density weights" begin
    xs = [0.0, 1.0]
    ys = [0.0]
    zs = [0.0]
    weights = ones(length(xs), length(ys), length(zs))
    density = reshape([2.0, 3.0], length(xs), length(ys), length(zs))
    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    potential = [4.0, 5.0]

    @test integrate_target_potential(vs, potential) ≈ (2.0 * 4.0 + 3.0 * 5.0) atol = 1e-12
end

@testset "volume potential matches direct sum on synthetic data" begin
    xs = [-0.25, 0.25]
    ys = [-0.25, 0.25]
    zs = [-0.25, 0.25]
    weights = fill(0.25, length(xs), length(ys), length(zs))
    density = fill(1.0, length(xs), length(ys), length(zs))
    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    targets = [(-1.0, 0.0, 0.5), (1.0, 0.0, 0.5)]

    direct = [
        sum(vs.weights[j] * vs.density[j] * BI.laplace3d_pot((vs.positions[1, j], vs.positions[2, j], vs.positions[3, j]), target) for j in eachindex(vs.density))
        for target in targets
    ]
    approx = evaluate_volume_potential(vs, targets; tol = 1e-8, kmax = 30.0)

    @test approx ≈ direct rtol = 5e-3 atol = 1e-8
end

@testset "pair interaction evaluates the full target set in one call" begin
    xs = [0.0, 1.0, 2.0]
    ys = [0.0]
    zs = [0.0]
    weights = ones(length(xs), length(ys), length(zs))
    source_density = reshape([1.0, 1.0, 1.0], length(xs), length(ys), length(zs))
    target_density = reshape([2.0, 3.0, 4.0], length(xs), length(ys), length(zs))

    source_vs = BI.VolumeSource((xs, ys, zs), weights, source_density)
    target_vs = BI.VolumeSource((xs, ys, zs), weights, target_density)

    volume_calls = Ref(0)
    scatter_calls = Ref(0)
    volume_target_count = Ref(0)
    scatter_target_count = Ref(0)

    fake_volume(source, targets; tol = 1e-6, kmax = nothing) = begin
        volume_calls[] += 1
        volume_target_count[] = size(targets, 2)
        return fill(1.5, size(targets, 2))
    end

    fake_scatter(interface, sigma, targets; fmm_tol = 1e-6, hcubature_atol = 1e-6, range_factor = 5.0, include_edges_src = false) = begin
        scatter_calls[] += 1
        scatter_target_count[] = size(targets, 2)
        return fill(-0.25, size(targets, 2))
    end

    result = evaluate_screened_pair_interaction(
        source_vs,
        target_vs,
        nothing,
        [1.0];
        volume_evaluator = fake_volume,
        scatter_evaluator = fake_scatter,
    )

    target_weights = target_vs.weights .* target_vs.density
    @test volume_calls[] == 1
    @test scatter_calls[] == 1
    @test volume_target_count[] == length(target_vs.density)
    @test scatter_target_count[] == length(target_vs.density)
    @test result.u_int_raw ≈ 1.5 * sum(target_weights) atol = 1e-12
    @test result.u_scatter_raw ≈ -0.25 * sum(target_weights) atol = 1e-12
    @test result.u_total_raw ≈ 1.25 * sum(target_weights) atol = 1e-12
end
