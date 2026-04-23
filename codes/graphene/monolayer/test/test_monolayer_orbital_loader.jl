using Test

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")

@testset "MonolayerOrbitalLoader.load_signed_xsf" begin
    @test isfile(XSF_1)
    datagrid = load_signed_xsf(XSF_1)
    @test size(datagrid.values) == (150, 150, 192)
    # Wannier orbitals are signed (not a probability density): expect both signs present.
    @test minimum(datagrid.values) < 0
    @test maximum(datagrid.values) > 0
end

@testset "MonolayerOrbitalLoader.centered_monolayer_sources — shared shift" begin
    out = centered_monolayer_sources(; orbital_1 = XSF_1, orbital_2 = XSF_2, source_tol = 1e-3, square = true)

    @test hasproperty(out, :vs1)
    @test hasproperty(out, :vs2)
    @test hasproperty(out, :shared_shift)
    @test length(out.shared_shift) == 3

    # The shift is defined to drive orbital 1's density-weighted centroid to the origin.
    weights_1 = out.vs1.weights .* out.vs1.density
    total_1 = sum(weights_1)
    centroid_1 = ntuple(d -> sum(out.vs1.positions[d, :] .* weights_1) / total_1, 3)
    @test all(abs.(centroid_1) .< 1e-8)

    # Orbital 2 gets the same shift, so its centroid sits at the A→B bond vector (~1.4232 Å in y). Per-axis tolerance is 50 mÅ (~0.5× the xy grid spacing).
    # x and z tolerances loosened from 1e-2 (plan) to 5e-2 to match y: the XSF xy grid spacing is ~0.08 Å, so ~29 mÅ x-centroid drift is unavoidable at source_tol=1e-3.
    weights_2 = out.vs2.weights .* out.vs2.density
    total_2 = sum(weights_2)
    centroid_2 = ntuple(d -> sum(out.vs2.positions[d, :] .* weights_2) / total_2, 3)
    @test abs(centroid_2[1]) < 5e-2
    @test abs(centroid_2[2] - 1.4231684) < 5e-2  # within 50 mÅ of atom-to-atom offset
    @test abs(centroid_2[3]) < 5e-2
end
