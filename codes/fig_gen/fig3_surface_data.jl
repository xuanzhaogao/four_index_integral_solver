#=
One-off solve for the +z-face σ used in fig3 panel (a). Independent of
fig3_data.jl (which runs the full edge-local / uniform sweep for panel (b)).
We solve a single edge-local configuration and serialize the panel mesh +
σ so the plotter can resample σ at arbitrary (x, y) without re-solving.

Run once; the artifact is reused by fig3_plot.jl.
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const p_quad = 4
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol = 1e-6
const gmres_tol = 1e-6

# Deep enough edge refinement that the smallest panels at the edge reach
# ~10^-3, giving the log-log resample several decades of clean range.
const l_ec = 1.01 / 2.0^7

iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
npan = length(iface.panels)
@info "Surface artifact mesh" l_ec npan N=npan * p_quad^2

Lhs = lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
rhs = rhs_dielectric_box3d(iface, ps, eps_src)
sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
gmres_res = norm(Lhs * sigma - rhs) / norm(rhs)
@info "Solve done" gmres_res

out = (
    Lx = Lx, Ly = Ly, Lz = Lz,
    p_quad = p_quad,
    l_ec = l_ec,
    eps_d = eps_d, eps_0 = eps_0,
    source = (point = ps.point, charge = ps.charge),
    interface = iface,
    sigma = sigma,
)

path = joinpath(@__DIR__, "fig3_surface_data.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path bytes=stat(path).size
