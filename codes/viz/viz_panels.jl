using BoundaryIntegral
import BoundaryIntegral as BI
using GLMakie

L = 20.0
Lx = L
Ly = L
Lz = 1.0

l_panel = 1.0
p = 4
r = 3
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, PointSource((6.0, 7.0, 1.5), 1.0), 1.0, l_ec, 1e-6, eps_in, eps_out, Float64, max_depth = 12)

fig = viz_3d(interface, show_points=false, highlight_edges = true)

neighbor_list = BI.build_neighbor_list(interface, max_order, 1e-6, false, false, range_factor = 5.0, distance_only = false)
fig = viz_3d(interface, highlight_panel = 100, show_points=false, highlight_edges = true, neighbor_list=neighbor_list)


interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, PointSource((1.0, 1.0, 1.0), 1.0), 1.0, l_ec, 1e-5, eps_in, eps_out, Float64)

fig = viz_3d(interface, show_points=true, highlight_edges = true)