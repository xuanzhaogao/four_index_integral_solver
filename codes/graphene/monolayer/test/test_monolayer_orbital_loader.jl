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
