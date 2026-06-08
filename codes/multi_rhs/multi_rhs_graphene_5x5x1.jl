# Multi-RHS demo on real graphene Wannier functions (5x5x1 supercell).
#
# Builds three pair densities from two Wannier orbitals,
#     rho_11 = phi_1 * phi_1,  rho_12 = phi_1 * phi_2,  rho_22 = phi_2 * phi_2,
# truncates them to their shared overlap support, centers them at the dielectric-slab
# midplane (graphene embedded in the slab, matching the monolayer reference geometry),
# builds ONE shared interface refined to resolve all three sources, and block-GMRES solves
# A Σ = F for the three layer densities at once.
#
# Data: the production cRPA monolayer run (k_323201_nb_144_c_15). The sheet is at z=7.5
# (mid-cell, c=15), so each Wannier function is a single contiguous blob — no periodic-wrap
# issue, naive centroid is correct (z=7.5), and the density is fine for a free-space solve.
#
# Run:  julia --project scripts/multi_rhs_graphene_5x5x1.jl

using BoundaryIntegral
using Krylov
using LinearAlgebra
const BI = BoundaryIntegral

# ----------------------------------------------------------------------------- parameters
const DATA_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1    = joinpath(DATA_DIR, "graphene_00001.xsf")
const XSF_2    = joinpath(DATA_DIR, "graphene_00002.xsf")

# Geometry / physics: matched to ~/work/four_index_integral_solver/codes/graphene/monolayer
# (screened_monolayer_hund_calibrated.jl). Graphene sits at the slab MIDPLANE z_center=Lz/2,
# embedded in the dielectric (eps_in); the pz lobes poke out of the thin slab into vacuum.
const L        = 90.0    # lateral slab size (Lx = Ly = L)
const LZ       = 3.35    # slab thickness
const Z_CENTER = LZ / 2  # graphene plane = slab midplane
const EPS_IN   = 3.5     # dielectric permittivity inside the slab
const EPS_OUT  = 1.0     # vacuum outside

support_rtol = 1e-4      # keep grid points where envelope >= rtol * max
# edge refinement target size, same formula as the reference solve_screened_mode:
#   l_ec = Lz / 2^edge_refine_level * 1.01.  Reference (production) uses level 4 (l_ec≈0.21,
#   cluster-scale on the L=90 rim); a smaller level keeps this demo interactive.
edge_refine_level = 2
l_ec         = LZ / 2.0^edge_refine_level * 1.01
n_quad       = 6         # Gauss-Legendre order per panel (reference value)
rhs_atol     = 1e-5      # RHS-driven refinement tolerance (reference rhs_tol)
fmm_tol      = 1e-5      # reference lhs_tol
up_tol       = 1e-5
max_order    = 64
gmres_rtol   = 1e-5      # reference gmres_rtol
gmres_itmax  = 200

# ------------------------------------------------------------------ 1. load orbitals
println("Loading orbitals ...")
_, dg1 = BI.read_xsf(XSF_1)
_, dg2 = BI.read_xsf(XSF_2)
@assert BI.datagrids_compatible(dg1, dg2) "orbitals must share the same grid"
println("  grid = ", (dg1.nx, dg1.ny, dg1.nz), "  (", dg1.nx * dg1.ny * dg1.nz, " points)")
println("  centroid phi_1 = ", round.(BI.density_centroid(dg1), digits = 3))
println("  centroid phi_2 = ", round.(BI.density_centroid(dg2), digits = 3))

# --------------------------------------------------- 2. pair densities + union truncation
println("Building pair densities rho_11, rho_12, rho_22 ...")
vs11 = VolumeSource(dg1, dg1.values .* dg1.values)
vs12 = VolumeSource(dg1, dg1.values .* dg2.values)
vs22 = VolumeSource(dg1, dg2.values .* dg2.values)

env  = sqrt.(vs11.density .^ 2 .+ vs12.density .^ 2 .+ vs22.density .^ 2)
keep = findall(>=(support_rtol * maximum(env)), env)
pos  = vs11.positions[:, keep]
wts  = vs11.weights[keep]

# Center the source like the reference: shift so the phi_1^2 centroid lands at (0,0,z_center)
# (orbital core at the slab midplane). density_centroid(dg1) = centroid of phi_1^2.
c1 = BI.density_centroid(dg1)
shift = (0.0 - c1[1], 0.0 - c1[2], Z_CENTER - c1[3])
pos[1, :] .+= shift[1]; pos[2, :] .+= shift[2]; pos[3, :] .+= shift[3]
println("  centering shift = ", round.(shift, digits = 3))

rho11 = VolumeSource(copy(pos), copy(wts), vs11.density[keep])
rho12 = VolumeSource(copy(pos), copy(wts), vs12.density[keep])
rho22 = VolumeSource(copy(pos), copy(wts), vs22.density[keep])
sources = [rho11, rho12, rho22]
println("  kept ", length(keep), " / ", length(env), " points (support_rtol = ", support_rtol, ")")

xmin, xmax = extrema(pos[1, :]); ymin, ymax = extrema(pos[2, :]); zmin, zmax = extrema(pos[3, :])
println("  centered source AABB: x[", round(xmin, digits = 2), ",", round(xmax, digits = 2),
        "] y[", round(ymin, digits = 2), ",", round(ymax, digits = 2),
        "] z[", round(zmin, digits = 2), ",", round(zmax, digits = 2), "]")

# ------------------------------------------------------------- 3. dielectric slab geometry
boxes = [(center = (0.0, 0.0, Z_CENTER), Lx = L, Ly = L, Lz = LZ)]
epses = [EPS_IN]
eps_out = EPS_OUT
println("Dielectric slab (reference geometry): center = (0, 0, ", Z_CENTER,
        ")  L = ", L, "  Lz = ", LZ, "  eps_in = ", EPS_IN, "  eps_out = ", EPS_OUT,
        "  (l_ec = ", round(l_ec, digits = 3), ", edge_refine_level = ", edge_refine_level, ")")

# --------------------------------------------------------- 4. build the shared interface
println("Refining interface to resolve all 3 sources ...")
interface = multi_dielectric_box3d_rhs_adaptive(
    n_quad, l_ec, boxes, epses, sources, rhs_atol; eps_out = eps_out)
println("  interface: ", length(interface.panels), " panels, ",
        BI.num_points(interface), " quadrature points")

# ----------------------------------------- 5. build operator + RHS ONCE (shared by both solves)
# A depends only on the interface (geometry + contrast), not on the source, so the block solve
# and the per-density single solves all use the SAME operator and the SAME right-hand sides.
op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, fmm_tol, up_tol, max_order)
F  = rhs_dielectric_box3d_fmm3d(interface, sources, fmm_tol)     # N x K
K  = size(F, 2)
labels = ["rho_11", "rho_12", "rho_22"]

# ----------------------------------------------------------------- 6. K=3 block GMRES solve
println("Block-GMRES solving A Σ = F  (K = ", K, " right-hand sides) ...")
Σ_block, bstats = Krylov.block_gmres(op, F; rtol = gmres_rtol, atol = 0.0, itmax = gmres_itmax)
println("  block: ", bstats.niter, " iters, solved = ", bstats.solved)

# ------------------------------------------------- 7. solve each density separately (same op/F)
println("Single-RHS GMRES solving each column separately (same interface/operator) ...")
Σ_single = similar(Σ_block)
single_iters = Int[]
for k in 1:K
    σk, sstats = Krylov.gmres(op, F[:, k]; rtol = gmres_rtol, atol = 0.0, itmax = gmres_itmax)
    Σ_single[:, k] .= σk
    push!(single_iters, sstats.niter)
end
println("  single-RHS iters per column = ", single_iters)

# ------------------------------------------------------------------- 8. compare Σ column by column
println("\nBlock vs per-density solve (same interface), per column:")
println(rpad("col", 8), rpad("||σ_block||", 16), rpad("||σ_single||", 16),
        rpad("max|Δ|", 14), "rel ||Δ||")
for k in 1:K
    b = @view Σ_block[:, k]
    s = @view Σ_single[:, k]
    nb = sqrt(sum(abs2, b)); ns = sqrt(sum(abs2, s))
    absd = maximum(abs.(b .- s))
    reld = sqrt(sum(abs2, b .- s)) / ns
    println(rpad(labels[k], 8), rpad(round(nb, sigdigits = 7), 16),
            rpad(round(ns, sigdigits = 7), 16),
            rpad(round(absd, sigdigits = 4), 14), round(reld, sigdigits = 4))
end

# ------------------------- 9. four-index integrals V_ab = ∫ rho_a (u_inc[rho_b] + u[sigma_b]) dx
# Same recipe as the monolayer reference (ScreenedOrbitalSolve): incident potential from the
# SCREENED source (TKM volume solve), scattered layer potential from sigma (FMM + hcubature near
# correction), contracted with the RAW target density. All sources share the grid, so the
# layer-potential operator and the per-source u_inc are built once and reused for block & single.
println("\nComputing V_ab = ∫ rho_a (u_inc[rho_b] + u[sigma_b]) dx  (", K, "x", K, ") ...")
targets = sources[1].positions
pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, targets, fmm_tol, fmm_tol, 5.0)

u_inc = Vector{Vector{Float64}}(undef, K)   # incident potential of each screened source
for b in 1:K
    sb = BI.screened_volume_source(interface, sources[b], BI.SharpScreening())
    kmax = BI._estimate_tkm3dc_kmax(sb)
    vals = BI.TKM3D.ltkm3dc(fmm_tol, sb.positions; charges = sb.weights .* sb.density,
                            targets = targets, pgt = 1, kmax = kmax)
    vals.ier == 0 || error("TKM3D.ltkm3dc failed, ier = $(vals.ier)")
    u_inc[b] = real.(vals.pottarg)
end
tw = [sources[a].weights .* sources[a].density for a in 1:K]   # raw target contraction weights

function Vmatrix(Σ)
    V = Matrix{Float64}(undef, K, K)
    for b in 1:K
        φb = u_inc[b] .+ (pottrg * Σ[:, b])
        for a in 1:K
            V[a, b] = dot(tw[a], φb)
        end
    end
    return V
end
V_block  = Vmatrix(Σ_block)
V_single = Vmatrix(Σ_single)

println("V_ab (raw, from block sigma):  [rows/cols = ", join(labels, ", "), "]")
for a in 1:K
    println("   ", join([rpad(round(V_block[a, b], sigdigits = 7), 16) for b in 1:K], ""))
end
println("V from block vs per-density single solve:")
println("  max |ΔV|     = ", maximum(abs.(V_block .- V_single)))
println("  max rel |ΔV| = ", maximum(abs.(V_block .- V_single) ./ abs.(V_single)))
println("done.")
