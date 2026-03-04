using FBCPoisson, BoxDMK, BoundaryIntegral
using FINUFFT
import BoundaryIntegral as BI

graphene_1 = joinpath(@__DIR__, "../../density_data/graphene_00001_5x5x1_shifted.xsf")
graphene_2 = joinpath(@__DIR__, "../../density_data/graphene_00002_5x5x1_shifted.xsf")

structure_1, datagrid_1 = BI.read_xsf(graphene_1)
datagrid_1.values .*= datagrid_1.values
vs_src_1 = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, - 7.920155482424242))
vs_src_1_truncated = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-6)

structure_2, datagrid_2 = BI.read_xsf(graphene_2)
datagrid_2.values .*= datagrid_2.values
vs_src_2 = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, - 7.920155482424242))
vs_src_2_truncated = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-6)

trg_pot_256 = lfbc3d(256, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_2_truncated.positions, 1e-6, 1)
U_01 = sum(trg_pot_256 .* vs_src_2_truncated.density .* vs_src_2_truncated.weights)

trg_pot_self = lfbc3d(256, vs_src_1_truncated.positions, vs_src_1_truncated.density .* vs_src_1_truncated.weights, vs_src_1_truncated.positions, 1e-6, 1)
U_00 = sum(trg_pot_self .* vs_src_1_truncated.density .* vs_src_1_truncated.weights)