#=
One-off deep solve for fig3 panel (b) convergence reference.

Solves the same problem as fig3_data_pscan.jl at the single deepest
configuration (p = 6, k = 13, l_ec = 1.01/2^13 ≈ 1.23e-4), evaluates the
far-field single-layer potential at the same target sphere, and serializes
just the u_far vector + metadata. The plotter loads this artifact and uses
it as the reference for all relative errors.

Keeping this solve separate from fig3_data_pscan.jls (which is already
populated) avoids re-running the 51 baseline solves.
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol  = 1e-6
const gmres_tol = 1e-6

const p_ref = 4
const k_ref = 13
const l_ec  = 1.01 / 2.0^k_ref

const R_far = 5.0
const n_theta = 12
const n_phi   = 24
far_targets = Matrix{Float64}(undef, 3, n_theta * n_phi)
let kk = 1
    for i in 1:n_theta
        θ = π * (i - 0.5) / n_theta
        for j in 1:n_phi
            φ = 2π * (j - 1) / n_phi
            far_targets[1, kk] = R_far * sin(θ) * cos(φ)
            far_targets[2, kk] = R_far * sin(θ) * sin(φ)
            far_targets[3, kk] = R_far * cos(θ)
            kk += 1
        end
    end
end

@info "Reference solve" p_ref k_ref l_ec
t0 = time()
iface = single_dielectric_box3d(Lx, Ly, Lz, p_ref, l_ec, eps_d, eps_0)
npan = length(iface.panels)
N = npan * p_ref^2
@info "  mesh" npan N

Lhs = lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
rhs = rhs_dielectric_box3d(iface, ps, eps_src)
sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
gmres_res = norm(Lhs * sigma - rhs) / norm(rhs)
@info "  GMRES done" gmres_res

S = BI.laplace3d_pottrg(iface, far_targets)
u_far = S * sigma
dt = time() - t0
@info "  far-field evaluated" dt

out = (
    p_ref = p_ref,
    k_ref = k_ref,
    l_ec  = l_ec,
    npanels = npan,
    N = N,
    far_R = R_far,
    n_theta = n_theta,
    n_phi   = n_phi,
    u_far = u_far,
    gmres_res = gmres_res,
)

path = joinpath(@__DIR__, "fig3_data_pref.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path bytes=stat(path).size
