#=
Surface-density sweep for Figure 3 panel (a).

For contrasts alpha = (eps_d - 1) / (eps_d + 1) in {10/33, 20/33, 30/33},
solve the corrected dielectric-cube BIE on the same edge-local mesh used for
the panel (a) profile. The plotter resamples the saved σ fields along a fixed
top-face slice.
=#

using LinearAlgebra
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

include("fig3_eps_utils.jl")

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_0 = 1.0
const p_quad = 4
const ps = BI.PointSource((0.1, 0.2, 0.3), 1.0)
const fmm_tol = 1e-6
const gmres_tol = 1e-6
const up_tol = 1e-6
const max_order = 64
const contrast_specs = [(num = 10, den = 33), (num = 20, den = 33), (num = 30, den = 33)]

# Same edge refinement as the existing single-contrast panel (a) artifact.
const l_ec = 1.01 / 2.0^7

surface_results = NamedTuple[]
for spec in contrast_specs
    alpha = spec.num / spec.den
    eps_d = Float64(epsilon_from_contrast(alpha; eps_0 = eps_0))
    @info "Surface contrast solve" alpha eps_d l_ec

    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    npan = length(iface.panels)
    @info "  mesh" npan N = npan * p_quad^2

    Lhs = lhs_dielectric_box3d_fmm3d_corrected(iface, fmm_tol, up_tol, max_order;
                                               correct_edges = true)
    rhs = rhs_dielectric_box3d(iface, ps, eps_d)
    sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    gmres_res = norm(Lhs * sigma - rhs) / norm(rhs)
    charge = charge_integral(BI.all_weights(iface), sigma)
    charge_theory = charge_theory_interior_source(eps_d; eps_0 = eps_0, source_charge = ps.charge)
    @info "  solved" gmres_res charge charge_theory charge_error_abs=abs(charge - charge_theory)

    push!(surface_results, (
        contrast_num = spec.num,
        contrast_den = spec.den,
        alpha = alpha,
        eps_d = eps_d,
        interface = iface,
        sigma = sigma,
        npanels = npan,
        N = npan * p_quad^2,
        l_min = panel_min_edge_length(iface),
        charge = charge,
        charge_theory = charge_theory,
        charge_error = charge - charge_theory,
        charge_error_abs = abs(charge - charge_theory),
        gmres_res = gmres_res,
    ))
end

out = (
    Lx = Lx, Ly = Ly, Lz = Lz,
    eps_0 = eps_0,
    eps_src = "eps_d",
    p_quad = p_quad,
    l_ec = l_ec,
    source = (point = ps.point, charge = ps.charge),
    fmm_tol = fmm_tol,
    gmres_tol = gmres_tol,
    up_tol = up_tol,
    max_order = max_order,
    correct_edges = true,
    surface_results = surface_results,
)

path = joinpath(@__DIR__, "fig3_surface_eps_sweep.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path bytes = stat(path).size
