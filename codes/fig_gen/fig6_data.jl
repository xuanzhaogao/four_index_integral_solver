#=
Data generation for Figure 6 (post-refinement for repeated layer-potential
evaluation). Sweeps the near-evaluation cutoff multiplier c (= range_factor
in r_near = c · h_P / n_quad) and the post-refinement threshold h₀.

For one fixed Gaussian volume target straddling the +z face of the cube:
  - solve BIE once, freeze σ on original Nyström grid;
  - for each (c, h₀), build a post-refined source mesh, evaluate the
    layer potential at the 64³ target grid via plain FMM (no near correction),
    and integrate against the Gaussian to get V_FMM;
  - reference V_ref uses FMM + hcubature with a generous range_factor;
  - report N_FMM (refined source DOF), N_HCub (number of near (target, panel)
    pairs), and E_V = |V_FMM − V_ref| / |V_ref|.
=#

using LinearAlgebra
using SparseArrays
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const p_quad = 4
const l_ec   = 0.125
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol     = 1e-10
const hcub_atol   = 1e-10
const gmres_tol   = 1e-10

# Target grid (64³ uniform) and Gaussian density centered on the +z face
const n_grid    = 50
const box_half  = 0.4
const gauss_center = (0.0, 0.0, Lz / 2)
const gauss_sigma  = 0.15

# Sweep
const c_list  = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
const h0_list = [0.5, 0.25, 0.125, 0.0625, 0.03125]
const c_ref   = 10.0      # generous c for the reference evaluation

# ---------------------------------------------------------------------------
# Solve BIE
# ---------------------------------------------------------------------------
@info "Building geometry and solving BIE ..."
const iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
@info "  panels=$(length(iface.panels))  N_src=$(BI.num_points(iface))"

Lhs = lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
rhs = rhs_dielectric_box3d(iface, ps, eps_src)
const sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
@info "  GMRES residual" rel = norm(Lhs * sigma - rhs) / norm(rhs)

# ---------------------------------------------------------------------------
# Build target grid and Gaussian density once
# ---------------------------------------------------------------------------
function build_targets_and_density()
    cx, cy, cz = gauss_center
    xs = collect(range(cx - box_half, cx + box_half; length = n_grid))
    ys = collect(range(cy - box_half, cy + box_half; length = n_grid))
    zs = collect(range(cz - box_half, cz + box_half; length = n_grid))
    targets = Matrix{Float64}(undef, 3, n_grid^3)
    rho = Vector{Float64}(undef, n_grid^3)
    dx = xs[2] - xs[1]
    dy = ys[2] - ys[1]
    dz = zs[2] - zs[1]
    norm_g = inv((sqrt(2π) * gauss_sigma)^3)
    kk = 1
    for i in 1:n_grid, j in 1:n_grid, k in 1:n_grid
        x, y, z = xs[i], ys[j], zs[k]
        targets[1, kk] = x; targets[2, kk] = y; targets[3, kk] = z
        r2 = (x - cx)^2 + (y - cy)^2 + (z - cz)^2
        rho[kk] = norm_g * exp(-r2 / (2 * gauss_sigma^2))
        kk += 1
    end
    return targets, rho, dx * dy * dz
end

const targets, rho_vals, dV = build_targets_and_density()
@info "  target grid" n = size(targets, 2) dV

# ---------------------------------------------------------------------------
# Reference V_ref
# ---------------------------------------------------------------------------
@info "Computing reference V_ref (c_ref = $c_ref) ..."
const op_ref = laplace3d_pottrg_fmm3d_corrected_hcubature(
    iface, targets, fmm_tol, hcub_atol, c_ref)
const u_ref = op_ref * sigma
const V_ref = sum(rho_vals .* u_ref) * dV
@info "  V_ref" V_ref

# ---------------------------------------------------------------------------
# Sweep (c, h₀)
# ---------------------------------------------------------------------------
function count_near_pairs(interface, c::Float64)
    nl = BI.build_target_neighbor_list(interface, targets, false; range_factor = c)
    return isempty(nl) ? 0 : sum(length(v) for v in values(nl))
end

results = Dict{Float64, Vector{NamedTuple}}()
for c in c_list
    @info "===== c = $c ====="
    rows = NamedTuple[]
    for h0 in h0_list
        @info "  h0 = $h0"
        # Post-refine source mesh
        ri, pids, splits = BI._refine_interface_for_targets(
            iface, targets, h0; range_factor = c)
        P = BI._refined_interface_prolongation(iface, ri, pids, splits)
        sig_r = P * sigma

        # FMM evaluation
        pot_base = BI.laplace3d_pottrg_fmm3d(ri, targets, fmm_tol)
        u_far = pot_base * sig_r

        # Hcubature near correction with the same cutoff c
        tnl = BI.build_target_neighbor_list(ri, targets, false; range_factor = c)
        N_HCub = isempty(tnl) ? 0 : sum(length(v) for v in values(tnl))
        if N_HCub > 0
            corr = BI.laplace3d_pottrg_corrections_hcubature(ri, targets, tnl, hcub_atol)
            u_near = corr * sig_r
        else
            u_near = zeros(Float64, length(u_far))
        end
        u = u_far .+ u_near
        V = sum(rho_vals .* u) * dV

        N_FMM = BI.num_points(ri)
        E_V   = abs(V - V_ref) / abs(V_ref)
        push!(rows, (h0 = h0, N_FMM = N_FMM, N_HCub = N_HCub, V = V, E_V = E_V))
        @info "    " N_FMM N_HCub E_V
    end
    results[c] = rows
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    geometry = (Lx = Lx, Ly = Ly, Lz = Lz),
    p_quad = p_quad,
    l_ec = l_ec,
    n_grid = n_grid,
    box_half = box_half,
    gauss_center = gauss_center,
    gauss_sigma = gauss_sigma,
    c_list = c_list,
    h0_list = h0_list,
    c_ref = c_ref,
    V_ref = V_ref,
    results = results,
)
datapath = joinpath(@__DIR__, "fig6_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes = stat(datapath).size

println("\n=== Sweep summary ===")
for c in c_list
    println("c = $c:")
    for r in results[c]
        println("  h0=$(r.h0)  N_FMM=$(r.N_FMM)  N_HCub=$(r.N_HCub)  E_V=$(round(r.E_V; sigdigits=4))")
    end
end
