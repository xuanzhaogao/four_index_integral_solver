#=
Data generation for Figure 7: end-to-end convergence of the screened
four-index integral on the stacked-box dielectric model.

Geometry (matches slab_model.md):
    Omega_2 = [-5,5] x [-5,5] x [ 0, 10]   (substrate, ε_2 = 10)
    Omega_1 = [-5,5] x [-5,5] x [10, 11]   (thin layer, ε_1 = 4)
    shared face Γ_12 = [-5,5]^2 × {10}
    exterior vacuum, ε_0 = 1

Densities (smooth normalized Gaussians, width s = 0.06):
    x_s = (-0.35, 0, 11.3),  x_t = ( 0.35, 0, 11.3)

Quantity of interest:
    V = V_vac + ΔV,   V_vac = ∫∫ ρ_tar(x) ρ_src(y) / (4π|x-y|) dx dy
                            = erf(R/(s_eff √2)) / (4π R), s_eff = s√2
    ΔV = ∫ ρ_tar(x) φ_scat(x) dx
       ≈ Σ_i ρ_tar(x_i) * φ_scat(x_i) * w_i,
    where {x_i, w_i, ρ_tar(x_i)} are the GaussianVolumeSource quadrature
    points / weights / sampled density for the target, and φ_scat(x_i) is
    the BIE-reconstructed single-layer potential evaluated at x_i.

Sweep: p ∈ p_list, r ∈ r_list ∪ {r_ref}. Reference at (p_ref, r_ref).
Records V, ΔV, GMRES n_iter, residual, panel count, DOF per (p, r).
=#

using LinearAlgebra
using SpecialFunctions
using Serialization
using Krylov
using BoundaryIntegral
const BI = BoundaryIntegral

# ---- Geometry (multi_dielectric_box3d_rhs_adaptive expects centered boxes) -
const SUB_CENTER   = (0.0, 0.0, 5.0)
const SUB_SIZE     = (10.0, 10.0, 10.0)
const LAYER_CENTER = (0.0, 0.0, 10.5)
const LAYER_SIZE   = (10.0, 10.0, 1.0)

const BOXES = [
    (center = SUB_CENTER,   Lx = SUB_SIZE[1],   Ly = SUB_SIZE[2],   Lz = SUB_SIZE[3]),
    (center = LAYER_CENTER, Lx = LAYER_SIZE[1], Ly = LAYER_SIZE[2], Lz = LAYER_SIZE[3]),
]
const EPSES   = [10.0, 4.0]   # matches BOXES order: substrate, layer
const EPS_OUT = 1.0           # exterior vacuum; the screened-source overloads
                              #   pick the per-sample-point permittivity from
                              #   (BOXES, EPSES, EPS_OUT) automatically.

# ---- Source / target Gaussian densities ------------------------------------
const S_GAUSS  = 0.06
const X_S      = (-0.35, 0.0, 11.3)
const X_T      = ( 0.35, 0.0, 11.3)
const N_VOL    = 12
const VOL_TOL  = 1e-6

const vs_src = BI.GaussianVolumeSource(X_S, S_GAUSS, N_VOL, VOL_TOL)
const vs_tar = BI.GaussianVolumeSource(X_T, S_GAUSS, N_VOL, VOL_TOL)

# Analytic vacuum reference V_vac for two normalized 3D Gaussians.
function v_vac_two_gaussians(x_s, x_t, s_s, s_t)
    R = sqrt(sum((x_s .- x_t).^2))
    s_eff = sqrt(s_s^2 + s_t^2)
    return erf(R / (s_eff * sqrt(2.0))) / (4π * R)
end
const V_VAC = v_vac_two_gaussians(X_S, X_T, S_GAUSS, S_GAUSS)
@info "V_vac (analytic)" V_VAC

# ---- Sweep and tolerances --------------------------------------------------
const FMM_TOL    = 1.0e-6
const GMRES_TOL  = 1.0e-6
const RHS_ATOL   = 1.0e-5   # looser than fmm_tol: avoid runaway RHS-adaptive
const RHS_MAXDEPTH = 8       # safety cap on adaptive depth
const MAX_ORDER  = 128
const LEC_BASE   = 1.0
const RANGE_FAC  = 5.0   # hcubature near-field correction range factor

const p_list = [2, 3, 4]
const r_list = [1, 2, 3, 4]
const p_ref  = 6
const r_ref  = 6   # r_ref = 10 is very expensive but ensures a well-converged reference solution

# ---- Per-configuration solve ----------------------------------------------
function delta_v_from_sigma(interface, sigma, vs_tar)
    targets = vs_tar.positions   # already a 3 × N_tar matrix
    pot_op  = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
        interface, targets, FMM_TOL, FMM_TOL, RANGE_FAC;
        include_edges_src = false)
    phi_scat = vec(pot_op * sigma)
    return sum(@. vs_tar.density * phi_scat * vs_tar.weights)
end

function solve_one(p::Int, r::Int)
    l_ec = LEC_BASE / 2.0^r
    @info "configuration" p r l_ec
    t0 = time()

    interface = BI.multi_dielectric_box3d_rhs_adaptive(
        p, l_ec, BOXES, EPSES, vs_src, RHS_ATOL;
        eps_out = EPS_OUT, max_depth = RHS_MAXDEPTH)
    n_pts   = BI.num_points(interface)
    n_panel = length(interface.panels)
    N_dof   = n_panel * p^2
    @info "  mesh" panels = n_panel n_pts N_dof

    @info "  assembling BIE operators with FMM" tol = FMM_TOL max_order = MAX_ORDER
    lhs = BI.lhs_dielectric_box3d_fmm3d_corrected(
        interface, FMM_TOL, FMM_TOL, MAX_ORDER;
        include_edges_src = false, include_edges_trg = false)
    
    @info "  assembling RHS with TKM/FMM hybrid (multi-box screened)" fmm_tol = FMM_TOL
    rhs = BI.rhs_dielectric_box3d_hybrid(interface, BOXES, EPSES, EPS_OUT, vs_src, FMM_TOL)

    @info "  solving linear system with GMRES" tol = GMRES_TOL

    sigma, stats = Krylov.gmres(lhs, rhs;
                                rtol = GMRES_TOL, atol = GMRES_TOL,
                                verbose = 1)
    n_iter   = stats.niter
    residual = norm(lhs * vec(sigma) - rhs) / max(norm(rhs), eps(Float64))
    @info "  GMRES" n_iter residual

    ΔV = delta_v_from_sigma(interface, vec(sigma), vs_tar)
    V  = V_VAC + ΔV
    dt = time() - t0
    @info "  result" ΔV V dt

    return (p = p, r = r, l_ec = l_ec, n_panels = n_panel,
            n_pts = n_pts, N = N_dof,
            n_iter = n_iter, residual = residual,
            ΔV = ΔV, V = V, dt = dt)
end

# ---- Reference solve first -------------------------------------------------
@info "=== Reference solve ===" p_ref r_ref
ref_result = solve_one(p_ref, r_ref)
V_REF = ref_result.V
@info "Reference V" V_REF V_VAC ΔV_ref = ref_result.ΔV

# ---- Convergence sweep -----------------------------------------------------
results = NamedTuple[]
for p in p_list
    for r in r_list
        push!(results, solve_one(p, r))
    end
end

# Append reference for downstream consumption (last entry).
push!(results, merge(ref_result, (; is_ref = true)))

out = (
    geometry = (boxes = BOXES, epses = EPSES, eps_out = EPS_OUT),
    density  = (s_gauss = S_GAUSS, x_s = X_S, x_t = X_T,
                n_vol = N_VOL, vol_tol = VOL_TOL),
    sweep    = (p_list = p_list, r_list = r_list, p_ref = p_ref, r_ref = r_ref),
    tols     = (fmm_tol = FMM_TOL, gmres_tol = GMRES_TOL,
                lec_base = LEC_BASE, range_fac = RANGE_FAC),
    V_vac    = V_VAC,
    V_ref    = V_REF,
    results  = results,
)
path = joinpath(@__DIR__, "fig7_data.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path bytes=stat(path).size

println("\n=== Summary ===")
println("V_vac (analytic)     = $(round(V_VAC; sigdigits=8))")
println("V_ref (p=$p_ref, r=$r_ref) = $(round(V_REF; sigdigits=8))")
println("\n  p    r    N        n_iter   E_V              ΔV")
for res in results
    e_v = abs(res.V - V_REF) / abs(V_REF)
    println("  $(res.p)    $(res.r)    $(res.N)    $(res.n_iter)    $(round(e_v; sigdigits=4))    $(round(res.ΔV; sigdigits=4))")
end
