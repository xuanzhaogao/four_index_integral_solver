using Test
using LinearAlgebra

isdefined(Main, :MonolayerOrbitalLoader) || include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
isdefined(Main, :MonolayerBareIntegrals) || include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "channel_pair_sources — shape and norms" begin
    densities = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)
    signed    = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = false)

    for channel in (:onsite, :nn, :hund_sf, :hund_ph)
        pair = channel_pair_sources(densities, signed, channel)
        @test hasproperty(pair, :source)
        @test hasproperty(pair, :target)
        @test size(pair.source.positions, 1) == 3
        @test size(pair.target.positions, 1) == 3
    end

    # The Wannier90 XSF is not unit-normalized on the supercell grid; the `to_eV` conversion
    # in bare_channel_integral divides by the actual norms, so any finite positive norm is fine.
    # What we *do* check: onsite source/target are the same VolumeSource (same norm), and the
    # Hund's product integrates to near zero because phi1 and phi2 are orthogonal Wannier orbitals.
    onsite = channel_pair_sources(densities, signed, :onsite)
    n_onsite = sum(onsite.source.weights .* onsite.source.density)
    @test n_onsite > 0
    @test onsite.source === onsite.target  # onsite uses vs1 on both legs

    nn = channel_pair_sources(densities, signed, :nn)
    n_nn_src = sum(nn.source.weights .* nn.source.density)
    n_nn_tgt = sum(nn.target.weights .* nn.target.density)
    @test n_nn_src > 0
    @test n_nn_tgt > 0
    # Source is vs1 for both onsite and nn; their source norms must agree exactly.
    @test n_nn_src ≈ n_onsite

    # Hund's signed-product pair integrates to ≈ 0 because phi1 and phi2 are orthogonal Wannier orbitals.
    hund = channel_pair_sources(densities, signed, :hund_sf)
    n_hund = sum(hund.source.weights .* hund.source.density)
    @test abs(n_hund) / n_onsite < 5e-2  # Hund's product norm is ≤ 5% of the onsite density norm
end
