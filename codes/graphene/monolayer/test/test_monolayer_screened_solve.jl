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

@testset "monolayer_screened_sources — z-centering and product source" begin
    out = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = 1e-3, z_center = 1.675,
    )
    @test hasproperty(out, :vs1)
    @test hasproperty(out, :vs2)
    @test hasproperty(out, :vs_product)
    @test hasproperty(out, :Nphi1)
    @test hasproperty(out, :Nphi2)

    # vs1 density centroid should be at (0, 0, z_center) since the loader
    # zeros out the in-plane centroid and we add z_center on top.
    weights = out.vs1.weights .* out.vs1.density
    total = sum(weights)
    centroid = (
        sum(out.vs1.positions[1, :] .* weights) / total,
        sum(out.vs1.positions[2, :] .* weights) / total,
        sum(out.vs1.positions[3, :] .* weights) / total,
    )
    @test abs(centroid[1]) < 1e-3
    @test abs(centroid[2]) < 1e-3
    @test abs(centroid[3] - 1.675) < 1e-3

    # vs_product integrates to ~0 (orthogonal Wannier orbitals).
    n_prod = sum(out.vs_product.weights .* out.vs_product.density)
    @test abs(n_prod) < 0.1

    # Nphi1, Nphi2 are the orbital-norm integrals sum(w * |phi|^2) on the XSF grid.
    # For the k_323201 dataset this is ~78 (matches the Na column in
    # monolayer/data/bare_monolayer_kmesh_sweep.csv at k_323201).
    @test 50.0 < out.Nphi1 < 100.0
    @test 50.0 < out.Nphi2 < 100.0
end
