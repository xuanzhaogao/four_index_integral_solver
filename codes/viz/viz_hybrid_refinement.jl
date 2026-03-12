using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using GLMakie

orbital_file = joinpath(@__DIR__, "../../density_data/graphene_00001_5x5x1_shifted.xsf")
structure, datagrid = BI.read_xsf(orbital_file)

datagrid.values .*= datagrid.values

vs = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242), tol = 1e-4)


L = 90.0
Lx = L
Ly = L
Lz = 2.4

l_panel = 1.0
p = 6
eps_out = 1.0
eps_in = 6.0

max_order = 128
l_ec = 10.0 / 2^4 * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64, tkm_kmax = 25.18050978)

fig = BI.viz_3d(; interfaces = [interface], sources = [vs], show_points = false, highlight_edges = true)

save(joinpath(@__DIR__, "figs/hybrid_refinement.png"), fig)

fig = BI.viz_3d(; interfaces = [interface], sources = [], show_points = false, highlight_edges = true)

save(joinpath(@__DIR__, "figs/hybrid_refinement_2.png"), fig)