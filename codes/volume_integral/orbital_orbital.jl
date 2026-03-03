using FBCPoisson, BoxDMK, BoundaryIntegral
using FINUFFT
import BoundaryIntegral as BI

graphene_1 = joinpath(@__DIR__, "../../density_data/graphene_00001_cropped.xsf")
structure_1, datagrid_1 = BI.read_xsf(graphene_1)
datagrid.values .*= datagrid.values
vs_src = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242))