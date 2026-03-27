using Test
using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie
using LaTeXStrings

include("../src/ScreenedDensityAnalysis.jl")
using .ScreenedDensityAnalysis

function synthetic_linear_datagrid(nx::Int, ny::Int, nz::Int)
    origin = (-0.5, -0.5, -0.5)
    A = (1.0, 0.0, 0.0)
    B = (0.0, 1.0, 0.0)
    C = (0.0, 0.0, 1.0)
    values = Array{Float64}(undef, nx, ny, nz)
    datagrid = (nx = nx, ny = ny, nz = nz, origin = origin, A = A, B = B, C = C, values = values)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        x, y, z = BI.grid_point(datagrid, ix, iy, iz)
        values[ix, iy, iz] = x + 2y + 3z
    end
    return datagrid
end

function synthetic_point_datagrid(nx::Int, ny::Int, nz::Int, active_index::NTuple{3, Int})
    origin = (0.0, 0.0, 0.0)
    A = (1.0, 0.0, 0.0)
    B = (0.0, 1.0, 0.0)
    C = (0.0, 0.0, 1.0)
    values = zeros(Float64, nx, ny, nz)
    values[active_index...] = 3.0
    return (nx = nx, ny = ny, nz = nz, origin = origin, A = A, B = B, C = C, values = values)
end

function synthetic_periodic_z_datagrid(nx::Int, ny::Int, nz::Int)
    origin = (-0.5, -0.5, -0.5)
    A = (1.0, 0.0, 0.0)
    B = (0.0, 1.0, 0.0)
    C = (0.0, 0.0, 1.0)
    values = Array{Float64}(undef, nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        θ = 2π * (iz - 1) / nz
        values[ix, iy, iz] = ix - 0.5 * iy + cos(θ) + 0.25 * sin(2θ)
    end
    return (nx = nx, ny = ny, nz = nz, origin = origin, A = A, B = B, C = C, values = values)
end

@testset "screened density samples on Cartesian points" begin
    datagrid = synthetic_linear_datagrid(4, 4, 4)
    points = [(0.0, 0.0, 0.0), (0.0, 0.2, -0.1)]
    values = screened_density_at_points(
        datagrid,
        points,
        BI.SharpScreening();
        bounds = ((-0.25, -0.25, -0.25), (0.25, 0.25, 0.25)),
        eps_in = 4.0,
        eps_out = 2.0,
        tol = 1e-12,
    )

    @test values[1] ≈ 0.0 atol = 1e-12
    @test values[2] ≈ ((0.0 + 2 * 0.2 + 3 * -0.1) / 4.0) atol = 1e-10
end

@testset "shifted datagrid is respected by Cartesian sampling" begin
    datagrid = synthetic_linear_datagrid(4, 4, 4)
    shifted = shift_datagrid(datagrid, (1.0, 0.0, 0.0))
    values = screened_density_at_points(
        shifted,
        [(1.0, 0.0, 0.0)],
        BI.SharpScreening();
        bounds = ((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0)),
        eps_in = 1.0,
        eps_out = 1.0,
        tol = 1e-12,
    )

    @test values[1] ≈ 0.0 atol = 1e-10
end

@testset "centering shift moves density centroid to origin" begin
    datagrid = synthetic_point_datagrid(4, 4, 4, (3, 2, 4))
    shift = centering_shift(datagrid; tol = 0.0)
    shifted = shift_datagrid(datagrid, shift)
    vs = BI.VolumeSource(shifted, tol = 0.0)
    centroid = density_centroid(vs)
    active_point = BI.grid_point(datagrid, 3, 2, 4)

    @test collect(centroid) ≈ [0.0, 0.0, 0.0] atol = 1e-12
    @test collect(shift) ≈ [-active_point[1], -active_point[2], -active_point[3]] atol = 1e-12
end

@testset "kz zero-mode reduction matches direct quadrature" begin
    xs = [0.0, 0.5]
    ys = [0.0, 0.5]
    zs = [-0.5, 0.0, 0.5]
    weights = ones(length(xs), length(ys), length(zs))
    density = Array{Float64}(undef, length(xs), length(ys), length(zs))
    for ix in eachindex(xs), iy in eachindex(ys), iz in eachindex(zs)
        density[ix, iy, iz] = 1 + zs[iz]
    end

    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    kz_values = [0.0, 0.75, 1.5]
    spectrum = kz_zero_mode_spectrum(vs, vs.density, kz_values)

    direct = [
        sum(vs.weights .* vs.density .* exp.(-im * kz .* vs.positions[3, :]))
        for kz in kz_values
    ]

    @test spectrum ≈ direct atol = 1e-12
end

@testset "TKM3D NUFFT kz-slab max spectrum is constant for a single active source" begin
    xs = [0.0, 1.0]
    ys = [0.0, 1.0]
    zs = [0.0, 1.0]
    weights = ones(length(xs), length(ys), length(zs))
    density = zeros(length(xs), length(ys), length(zs))
    density[1, 1, 1] = 1.0

    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    kz_values, slab_max = tkm3d_kz_slab_max_spectrum(vs, vs.density; kmax = 2.5, eps = 1e-12)

    @test !isempty(kz_values)
    @test first(kz_values) ≈ 0.0 atol = 1e-12
    @test all(kz_values .>= 0.0)
    @test slab_max ≈ ones(length(slab_max)) atol = 1e-10
end

@testset "global z-only NUFFT upsampling preserves original z planes" begin
    datagrid = synthetic_periodic_z_datagrid(2, 3, 4)
    factor = 3
    refined = ScreenedDensityAnalysis.global_nufft_upsample_datagrid_z(datagrid; factor = factor, eps = 1e-12)

    @test refined.nx == datagrid.nx
    @test refined.ny == datagrid.ny
    @test refined.nz == datagrid.nz * factor
    @test refined.values[:, :, 1:factor:end] ≈ datagrid.values atol = 1e-10

    for ix in 1:datagrid.nx, iy in 1:datagrid.ny, iz in 1:datagrid.nz
        original_point = BI.grid_point(datagrid, ix, iy, iz)
        refined_point = BI.grid_point(refined, ix, iy, 1 + (iz - 1) * factor)
        @test collect(refined_point) ≈ collect(original_point) atol = 1e-12
    end
end

@testset "3D NUFFT z-target resampling preserves original z planes" begin
    datagrid = synthetic_periodic_z_datagrid(2, 3, 4)
    factor = 2
    refined_values = ScreenedDensityAnalysis.nufft_resample_datagrid_z_via_3d(datagrid; factor = factor, eps = 1e-12)

    @test size(refined_values) == (datagrid.nx, datagrid.ny, datagrid.nz * factor)
    @test refined_values[:, :, 1:factor:end] ≈ datagrid.values atol = 1e-10
end

@testset "z-upsampled TKM3D decay extends kz coverage for refined factors" begin
    datagrid = synthetic_periodic_z_datagrid(2, 2, 6)
    bounds = ((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0))
    curves = ScreenedDensityAnalysis.z_upsampled_tkm3d_decay_curves(
        datagrid,
        bounds,
        1.0,
        1.0,
        BI.SharpScreening();
        upsample_factors = [1, 3],
        source_tol = 0.0,
        tol = 1e-12,
    )

    @test length(curves) == 2
    factor1, kz1, decay1 = curves[1]
    factor3, kz3, decay3 = curves[2]
    @test factor1 == 1
    @test factor3 == 3
    @test first(kz1) ≈ 0.0 atol = 1e-12
    @test first(kz3) ≈ 0.0 atol = 1e-12
    @test maximum(kz3) > maximum(kz1)
    @test first(decay1) ≈ 1.0 atol = 1e-12
    @test first(decay3) ≈ 1.0 atol = 1e-12
end

@testset "periodic spectral upsampling reproduces a cosine" begin
    n = 8
    factor = 4
    samples = [1.0 + 0.25 * cos(2π * j / n) for j in 0:(n - 1)]
    refined = spectral_upsample_periodic(samples; factor = factor)
    exact = [1.0 + 0.25 * cos(2π * j / (n * factor)) for j in 0:(n * factor - 1)]

    @test refined ≈ exact atol = 1e-10
end

@testset "periodic nufft upsampling reproduces a cosine" begin
    n = 8
    factor = 4
    samples = [1.0 + 0.25 * cos(2π * j / n) for j in 0:(n - 1)]
    refined = nufft_upsample_periodic(samples; factor = factor, eps = 1e-12)
    exact = [1.0 + 0.25 * cos(2π * j / (n * factor)) for j in 0:(n * factor - 1)]

    @test refined ≈ exact atol = 1e-10
end

@testset "global nufft upsampling preserves original samples" begin
    values = [0.0, 1.0, 2.0, 1.5, 0.5]
    factor = 4
    refined = ScreenedDensityAnalysis.nufft_upsample_global(values; factor = factor, pad_factor = 4, eps = 1e-12)

    @test length(refined) == length(values) * factor
    @test refined[1:factor:end] ≈ values atol = 1e-10
end

@testset "global nufft upsampling reduces graphene qz boundary ringing" begin
    _, vs = default_volume_source()
    _, qz, _ = plane_charge_density(vs, vs.density)
    periodic = nufft_upsample_periodic(qz; factor = 8, eps = 1e-12)
    global_refined = ScreenedDensityAnalysis.nufft_upsample_global(qz; factor = 8, pad_factor = 4, eps = 1e-12)

    @test minimum(global_refined) > minimum(periodic)
end

@testset "linear upsampling reproduces a linear profile" begin
    values = [1.0, 2.0, 3.0, 4.0]
    refined = linear_upsample_uniform(values; factor = 4)
    exact = collect(range(1.0, 4.0; length = (length(values) - 1) * 4 + 1))

    @test refined ≈ exact atol = 1e-12
end

@testset "refined qz spectrum matches coarse quadrature without screening" begin
    xs = [0.0, 0.5]
    ys = [0.0, 0.5]
    zs = [-0.5, 0.0, 0.5, 1.0]
    weights = ones(length(xs), length(ys), length(zs))
    density = Array{Float64}(undef, length(xs), length(ys), length(zs))
    for ix in eachindex(xs), iy in eachindex(ys), iz in eachindex(zs)
        density[ix, iy, iz] = 1 + 0.1 * zs[iz]
    end

    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    kz_values = [0.0, 0.5, 1.0]
    coarse = kz_zero_mode_spectrum(vs, vs.density, kz_values)
    refined = refined_kz_zero_mode_spectrum(
        vs,
        vs.density,
        ((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0)),
        1.0,
        1.0,
        BI.SharpScreening(),
        kz_values;
        upsample_factor = 1,
        tol = 0.0,
    )

    @test refined ≈ coarse atol = 1e-12
end

@testset "global nufft qz spectrum matches coarse quadrature without screening" begin
    xs = [0.0, 0.5]
    ys = [0.0, 0.5]
    zs = [-0.5, 0.0, 0.5, 1.0]
    weights = ones(length(xs), length(ys), length(zs))
    density = Array{Float64}(undef, length(xs), length(ys), length(zs))
    for ix in eachindex(xs), iy in eachindex(ys), iz in eachindex(zs)
        density[ix, iy, iz] = 1 + 0.1 * zs[iz]
    end

    vs = BI.VolumeSource((xs, ys, zs), weights, density)
    kz_values = [0.0, 0.5, 1.0]
    coarse = kz_zero_mode_spectrum(vs, vs.density, kz_values)
    refined = ScreenedDensityAnalysis.global_nufft_kz_zero_mode_spectrum(
        vs,
        vs.density,
        ((-10.0, -10.0, -10.0), (10.0, 10.0, 10.0)),
        1.0,
        1.0,
        BI.SharpScreening(),
        kz_values;
        upsample_factor = 1,
        pad_factor = 1,
        tol = 0.0,
    )

    @test refined ≈ coarse atol = 1e-12
end

@testset "log10 slice transform clamps nonpositive values" begin
    values = [1.0 0.1 NaN; 0.0 -2.0 1e-4]
    transformed = log10_clamped(values; floor_value = 1e-3)

    @test transformed[1, 1] ≈ 0.0 atol = 1e-12
    @test transformed[1, 2] ≈ -1.0 atol = 1e-12
    @test isnan(transformed[1, 3])
    @test transformed[2, 1] ≈ -3.0 atol = 1e-12
    @test transformed[2, 2] ≈ -3.0 atol = 1e-12
    @test transformed[2, 3] ≈ -3.0 atol = 1e-12
end

@testset "screened density script can be included without running main" begin
    mod = Module(:ScreenedDensityScriptIncludeTest)
    Core.eval(mod, :(include(path) = Base.include($mod, path)))
    Core.eval(mod, :(using BoundaryIntegral))
    Core.eval(mod, :(import BoundaryIntegral as BI))
    Core.eval(mod, :(using CairoMakie))
    Core.eval(mod, :(using LaTeXStrings))
    Base.include(mod, joinpath(@__DIR__, "..", "scripts", "screened_density_tkm.jl"))

    @test isdefined(mod, :main_all)
    @test isdefined(mod, :main_z_upsampled)
end
