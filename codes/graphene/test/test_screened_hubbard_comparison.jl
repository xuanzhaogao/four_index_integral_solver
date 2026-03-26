using Test
using DataFrames

include("../src/ScreenedHubbardComparison.jl")
using .ScreenedHubbardComparison

@testset "graphene cRPA reference values match Table I" begin
    @test GRAPHENE_CRPA_EV["U_00"] == 9.3
    @test GRAPHENE_CRPA_EV["U_01"] == 5.5
    @test GRAPHENE_CRPA_EV["U_02"] == 4.1
    @test GRAPHENE_CRPA_EV["U_03"] == 3.6
end

@testset "pair curve data is sorted by pair order and bandwidth" begin
    df = DataFrame(
        mode = ["SoftMixInversePermittivity", "SoftMixInversePermittivity", "SoftMixInversePermittivity", "SoftMixInversePermittivity"],
        bandwidth = [0.4, 0.1, 0.4, 0.1],
        pair = ["U_01", "U_01", "U_00", "U_00"],
        u_total_ev = [5.9, 5.8, 9.8, 9.6],
    )

    curves = pair_curve_data(df)

    @test [curve.pair for curve in curves] == ["U_00", "U_01"]
    @test curves[1].bandwidths == [0.1, 0.4]
    @test curves[1].values == [9.6, 9.8]
    @test curves[2].bandwidths == [0.1, 0.4]
    @test curves[2].values == [5.8, 5.9]
end
