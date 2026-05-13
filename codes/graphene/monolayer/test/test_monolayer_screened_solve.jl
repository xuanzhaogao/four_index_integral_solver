using Test

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const COQUI_OUT = joinpath(REF_DIR, "_coqui_crpa_loc.out")
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "parse_coqui_loc" begin
    @test isfile(COQUI_OUT)
    out = parse_coqui_loc(COQUI_OUT)
    @test sort(collect(keys(out))) == [:hund_ph, :hund_sf, :nn, :onsite]
    @test out[:onsite].v_ijkl ≈ 17.4029 atol = 1e-4
    @test out[:onsite].U_ijkl ≈ 9.7824  atol = 1e-4
    @test out[:nn].v_ijkl     ≈ 8.8319  atol = 1e-4
    @test out[:nn].U_ijkl     ≈ 5.1678  atol = 1e-4
    @test out[:hund_ph].v_ijkl ≈ 0.1316 atol = 1e-4
    @test out[:hund_ph].U_ijkl ≈ 0.0937 atol = 1e-4
    @test out[:hund_sf].v_ijkl ≈ 0.1316 atol = 1e-4
    @test out[:hund_sf].U_ijkl ≈ 0.0937 atol = 1e-4
end
