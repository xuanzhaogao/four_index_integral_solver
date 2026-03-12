using Test

include(joinpath(@__DIR__, "..", "tkm3d_hubbard_utils.jl"))

@testset "tolerance sweep rows use tightest run as reference" begin
    samples = [
        (tol = 1.0e-2, U_raw = 8.0, U_ev = 80.0),
        (tol = 1.0e-4, U_raw = 9.0, U_ev = 90.0),
    ]

    rows = annotate_reference_errors("U_01", samples)

    @test length(rows) == 2
    @test rows[end].rel_err_raw == 0.0
    @test rows[end].rel_err_ev == 0.0
    @test rows[1].rel_err_raw == abs(8.0 - 9.0) / 9.0
    @test rows[1].rel_err_ev == abs(80.0 - 90.0) / 90.0
end

@testset "normalize_tolerances sorts and validates inputs" begin
    @test normalize_tolerances([1.0e-4, 1.0e-2, 1.0e-3]) == [1.0e-2, 1.0e-3, 1.0e-4]
    @test_throws ArgumentError normalize_tolerances(Float64[])
    @test_throws ArgumentError normalize_tolerances([1.0e-2, 1.0, 1.0e-4])
end
