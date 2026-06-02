#=
p=2 self-convergence probe for fig3 panel (b).

Fix the polynomial order at p=2 and increase the edge-refinement depth
k = 1, 2, 3, ... (l_ec = 1.01/2^k). At each level solve the dielectric-box
BIE and evaluate the single-layer potential u_k at a FIXED set of exterior
test points lying on the diagonal plane y = x of the box. The points wrap
around the box cross-section (corners, edge-bands, face-bands) at offsets
log-spaced from 1e-2 (near field) to 5 (far field), so the metric is
sensitive to BOTH near- and far-field accuracy.

Near targets make the plain smooth-quadrature evaluator inaccurate, so the
potential is evaluated with the near-field-corrected operator
laplace3d_pottrg_fmm3d_corrected_hcubature (FMM base + hcubature corrections
for near target-panel pairs). The same evaluator is used at every level, so
the successive change reflects sigma convergence, not evaluation noise.

Self-convergence indicator:
    delta_k = ||u_k - u_{k-1}||_2 / ||u_{k-1}||_2 .

CHECKPOINTING: after every completed level the full results are re-serialized
to fig3_p2_selfconv.jls and a flushed one-line summary is appended to
fig3_p2_selfconv_progress.txt, so an OOM-killed deep level loses nothing.

Stops when delta_k < floor_tol or stops decreasing, capped at k_max = 16.
=#

using LinearAlgebra
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol  = 1e-6
const gmres_tol = 1e-6

# near-field-corrected potential evaluator settings
const hcub_atol     = 1e-8
const range_factor  = 5.0

const p_quad = 2
const k_max  = 16          # hard cap on refinement depth
const floor_tol = 2e-6     # delta below this ⇒ at tolerance floor, stop
const patience  = 1        # stop after this many non-decreasing deltas (once below 1e-3)

const jls_path  = joinpath(@__DIR__, "fig3_p2_selfconv.jls")
const prog_path = joinpath(@__DIR__, "fig3_p2_selfconv_progress.txt")

# ---------------------------------------------------------------------------
# Fixed test points on the diagonal plane y = x, exterior to the box.
# A point on the plane is (a, a, z); the box is {|a|<=1, |z|<=1}. Distances:
#   edge band  (|a|>1, |z|<=1): dist = sqrt(2)*(|a|-1)   ⇒ a = 1 + off/sqrt(2)
#   face band  (|a|<=1, |z|>1): dist = |z|-1             ⇒ z = 1 + off
#   corner     (|a|>1, |z|>1) : sqrt(2(|a|-1)^2 + (|z|-1)^2)
# Offsets are log-spaced from 1e-2 (near) to 5 (far). A few interior
# coordinates are added so points sit directly off faces/edges; the interior
# rectangle is then dropped, leaving an exterior grid wrapping the box.
# ---------------------------------------------------------------------------
const n_off   = 9
const offs    = 10 .^ range(log10(1e-2), log10(5.0); length = n_off)
const a_inner = collect(range(-0.9, 0.9; length = 3))
const z_inner = collect(range(-0.9, 0.9; length = 3))

const a_list = sort(unique(vcat(-(1 .+ offs ./ sqrt(2)), a_inner, 1 .+ offs ./ sqrt(2))))
const z_list = sort(unique(vcat(-(1 .+ offs), z_inner, 1 .+ offs)))

function build_plane_targets()
    pts = NTuple{3, Float64}[]
    for a in a_list, z in z_list
        if abs(a) > 1 + 1e-12 || abs(z) > 1 + 1e-12     # exterior only
            push!(pts, (a, a, z))
        end
    end
    T = Matrix{Float64}(undef, 3, length(pts))
    for (j, p) in enumerate(pts)
        T[1, j], T[2, j], T[3, j] = p
    end
    return T
end

const targets = build_plane_targets()

function solve_and_far(l_ec::Float64)
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    Lhs = lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
    rhs = rhs_dielectric_box3d(iface, ps, eps_src)
    sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    gres = norm(Lhs * sigma - rhs) / norm(rhs)
    M = laplace3d_pottrg_fmm3d_corrected_hcubature(iface, targets, fmm_tol, hcub_atol, range_factor)
    return length(iface.panels), M * sigma, gres
end

# checkpoint helper: rewrite the .jls with everything collected so far
function checkpoint(levels, stop_reason)
    out = (p_quad = p_quad, targets = targets, n_targets = size(targets, 2),
           offs = offs, far_R = 5.0, floor_tol = floor_tol,
           stop_reason = stop_reason, levels = levels)
    open(io -> serialize(io, out), jls_path, "w")
end

open(prog_path, "w") do io
    println(io, "# p=2 self-convergence on diagonal plane y=x  ($(size(targets,2)) targets, off 1e-2..5)")
    println(io, "# fmm_tol=$fmm_tol gmres_tol=$gmres_tol hcub_atol=$hcub_atol range_factor=$range_factor k_max=$k_max")
    println(io, "# k    l_ec          npanels       N            delta_k        gmres_res     dt[s]")
    flush(io)
end

levels = NamedTuple[]
u_prev = nothing
delta_prev = Inf
n_nondecr = 0
stop_reason = "k_max"

for k in 1:k_max
    l_ec = 1.01 / 2.0^k
    @info "level" p_quad k l_ec
    t0 = time()
    npan, u_far, gres = solve_and_far(l_ec)
    dt = time() - t0
    N = npan * p_quad^2
    delta = u_prev === nothing ? NaN : norm(u_far .- u_prev) / norm(u_prev)
    push!(levels, (k = k, l_ec = l_ec, npanels = npan, N = N,
                   u_far = u_far, delta = delta, gmres_res = gres, dt = dt))
    @info "  done" k npan N delta gmres_res=gres dt

    # --- checkpoint to disk + live progress line (flushed) ---
    checkpoint(levels, stop_reason)
    open(prog_path, "a") do io
        ds = isnan(delta) ? "—" : string(round(delta; sigdigits = 5))
        println(io, "  $k    $(round(l_ec; sigdigits=4))   $npan   $N   $ds   $(round(gres; sigdigits=3))   $(round(dt; digits=1))")
        flush(io)
    end

    global u_prev = u_far
    if !isnan(delta)
        if delta < floor_tol
            global stop_reason = "below floor_tol ($floor_tol)"
            checkpoint(levels, stop_reason)
            break
        end
        if delta >= delta_prev && delta < 1e-3
            global n_nondecr += 1
            if n_nondecr >= patience
                global stop_reason = "delta stopped decreasing"
                checkpoint(levels, stop_reason)
                break
            end
        else
            global n_nondecr = 0
        end
        global delta_prev = delta
    end
end

checkpoint(levels, stop_reason)
@info "stopped" stop_reason last_k=levels[end].k
open(prog_path, "a") do io
    println(io, "# stopped: $stop_reason  (last k = $(levels[end].k))")
    flush(io)
end
@info "Saved" jls_path bytes=stat(jls_path).size
