using GLMakie
using BoundaryIntegral
import BoundaryIntegral as BI

L = 20.0
Lx = L
Ly = L
Lz = 1.0

l_panel = 1.0
p = 6
r = 3
eps_out = 1.0
eps_in = 4.0

max_order = 128
l_ec = 1 / 2^r * 1.01

# Gaussian volume source on a tensor grid
center = (1.0, 2.0, 1.1)
σ = 0.1

# tol controls truncation of the Gaussian support
source = BI.GaussianVolumeSource(center, σ, 1e-8)

interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, source, 1.0, l_ec, 1e-6, eps_in, eps_out, Float64)

fig = BI.viz_3d(; interfaces = [interface], sources = [], min_density = 1e-6, algorithm = :mip, show_points = false, highlight_edges = true)
