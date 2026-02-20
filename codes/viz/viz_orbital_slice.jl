using BoundaryIntegral
import BoundaryIntegral as BI
using GLMakie

orbital_file = joinpath(@__DIR__, "../../density_data/graphene_00002_cropped.xsf")
structure, datagrid = BI.read_xsf(orbital_file)

z = 1.2
nz = Int(ceil((structure.primvec[3, 3] / 2 + 1.2) / structure.primvec[3, 3] * 100))
slice_val = datagrid.values[:, :, nz]

maximum(slice_val)