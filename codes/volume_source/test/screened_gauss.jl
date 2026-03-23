using Test

include(joinpath(@__DIR__, "..", "scripts", "screened_gauss.jl"))

@testset "erf-screened Gaussian decay helpers" begin
    sigma = 1.0
    narrow = 0.2
    wide = 1.0

    @test erf_screen_function(-5.0, narrow) < 1.0e-6
    @test erf_screen_function(5.0, narrow) > 1.0 - 1.0e-6
    @test erf_screen_function(0.0, narrow) ≈ 0.5
    @test erf_screen_function(0.5, narrow) > erf_screen_function(0.5, wide)

    study = run_erf_screen_decay_study(;
        sigma,
        bandwidths = [0.05, narrow],
        half_width_sigmas = 4.0,
        npoints = 32,
        refinement_half_width_sigmas = 0.5,
        refinement_factor = 8,
        mode_extent = 64,
    )

    @test study.bandwidths == [0.05, narrow]
    @test length(study.curves) == 2
    @test all(length(curve.k_positive) == length(curve.magnitude_positive) for curve in study.curves)
    @test all(all(curve.magnitude_positive .>= 0.0) for curve in study.curves)
    @test first(study.curves[1].k_positive) ≈ 0.0
    @test issorted(study.curves[1].k_positive)
    @test study.curves[1].k_positive == study.curves[2].k_positive
    @test study.min_weight < study.coarse_weight
    @test study.step_curve.label == "step"
    @test study.step_curve.k_positive == study.curves[1].k_positive
    @test all(study.step_curve.magnitude_positive .>= 0.0)
    @test study.unscreened_curve.label == "unscreened"
    @test study.unscreened_curve.k_positive == study.curves[1].k_positive
    @test all(study.unscreened_curve.magnitude_positive .>= 0.0)
end

@testset "refined quadrature and weighted NUFFT" begin
    quad = screen_locally_refined_quadrature(;
        sigma = 1.0,
        half_width_sigmas = 6.0,
        npoints = 32,
        refinement_half_width_sigmas = 0.5,
        refinement_factor = 8,
    )

    @test length(quad.x) == length(quad.weights)
    @test all(quad.weights .> 0.0)
    @test sum(quad.weights) ≈ 12.0 atol = 1.0e-12
    @test minimum(quad.weights) < maximum(quad.weights)

    values = gaussian_screen_base.(quad.x, Ref(1.0))
    k, fhat = screen_continuous_nufft_type1_weighted(
        quad.x,
        values,
        quad.weights,
        quad.interval_length;
        mode_extent = 32,
    )

    positive = findall(k .>= 0.0)
    reference = gaussian_screen_base_ft.(k[positive][1:5], Ref(1.0))

    @test maximum(abs.(real.(fhat[positive][1:5]) .- reference)) < 1.0e-3
    @test maximum(abs.(imag.(fhat[positive][1:5]))) < 1.0e-3
end

@testset "regularized line convolution" begin
    xmin = -2.0
    xmax = 2.0
    x0 = 0.25
    delta = 0.1

    result = regularized_line_convolution(
        x -> 1.0,
        x0;
        xmin,
        xmax,
        delta,
        rtol = 1.0e-10,
        atol = 1.0e-12,
    )

    expected = log(((x0 - xmin) * (xmax - x0)) / delta^2)

    @test result.value ≈ expected atol = 1.0e-8 rtol = 1.0e-8
    @test result.error >= 0.0
end

@testset "screened line convolution study" begin
    study = run_screened_line_convolution_study(;
        sigma = 1.0,
        bandwidths = [0.2],
        source_half_width_sigmas = 4.0,
        target_half_width_sigmas = 1.0,
        target_npoints = 9,
        exclusion_sigmas = 0.05,
        rtol = 1.0e-6,
        atol = 1.0e-8,
    )

    @test study.step_curve.label == "step"
    @test length(study.targets) == 9
    @test study.step_curve.targets == study.targets
    @test length(study.curves) == 1
    @test study.curves[1].targets == study.targets
    @test all(isfinite, study.step_curve.potential)
    @test all(isfinite, study.curves[1].potential)
    @test all(isfinite, study.curves[1].difference_from_step)
end
