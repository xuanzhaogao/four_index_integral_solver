# 6.1.5 / 6.1.6 / 6.1.7 — System II correctness checks + qualitative figure data.
# One solve at eps = 1e-9, p = 8, r = 4 shared by all checks.
#
# 6.1.7  shared-face census: the Omega1-Omega2 face must enter Gamma exactly once
#        with gamma12 = (eps1-eps2)/(eps1+eps2); explicit duplicate-panel scan.
# 6.1.6  transmission conditions [phi] = 0 and [eps dn phi] = 0 across the
#        Omega1-Omega2 face and the Omega_m top face; sum(sigma dS) per body;
#        Gauss's law on an enclosing box.
# 6.1.5  panelization (corners + level + sigma) and sigma along the top surface
#        approaching the triple-junction edge (x = 0, z = 0.5).

include(joinpath(@__DIR__, "..", "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf, Statistics
using FastGaussQuadrature: gausslegendre

const SYS = system2()
const DATA = joinpath(@__DIR__, "..", "data")
const EPS = parse(Float64, get(ENV, "CHK_EPS", "1e-9"))
const P = parse(Int, get(ENV, "CHK_P", "8"))
const R = parse(Int, get(ENV, "CHK_R", "4"))

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

println(">>> solve System II  eps=$(EPS) p=$(P) r=$(R)")
res = solve_system(SYS; eps = EPS, p = P, r = R)
@printf("N=%d niter=%d residual=%.2e\n", res.N, res.niter, res.residual)

iface = res.interface
panels = iface.panels
centers = [(p.corners[1] .+ p.corners[2] .+ p.corners[3] .+ p.corners[4]) ./ 4 for p in panels]
lens = [BI._panel_max_length(p) for p in panels]

# ---------------------------------------------------------------- 6.1.7 census
println(">>> 6.1.7 geometry sanity")
shared_ids = [i for i in eachindex(panels)
              if abs(centers[i][1]) < 1e-12 && abs(panels[i].normal[1]) > 0.999 &&
                 0.0 < centers[i][3] < 0.5]
gammas = unique([round((iface.eps_in[i] - iface.eps_out[i]) / (iface.eps_in[i] + iface.eps_out[i]); digits = 8)
                 for i in shared_ids])
# duplicate-panel scan: identical centers anywhere in Gamma
dup_count = let seen = Dict{NTuple{3, Float64}, Int}(), cnt = 0
    for c in centers
        key = (round(c[1]; digits = 10), round(c[2]; digits = 10), round(c[3]; digits = 10))
        cnt += (seen[key] = get(seen, key, 0) + 1) > 1 ? 1 : 0
    end
    cnt
end
@printf("shared-face panels: %d   gamma values: %s (expected |gamma|=0.5)\n", length(shared_ids), gammas)
@printf("duplicate-center panels anywhere: %d (expected 0)\n", dup_count)

# ------------------------------------------------------- 6.1.6 transmission
# Generic one-sided FD weights via Vandermonde: sum_i w_i f(t_i) ~ f^(deriv)(0).
function fd_weights(nodes::Vector{Float64}, deriv::Int)
    n = length(nodes)
    Vm = [nodes[j]^k for k in 0:(n - 1), j in 1:n]
    rhs = zeros(n); rhs[deriv + 1] = factorial(deriv)
    return Vm \ rhs
end

"Transmission residuals across a planar interface at x0 along unit normal nrm.
4-node one-sided stencils (O(d^3) derivative); the check itself floors around
max(d^3 scale, eval_tol/d) — reported alongside."
function transmission(res, pts0::Vector{NTuple{3, Float64}}, nrm::NTuple{3, Float64},
                      eps_minus::Float64, eps_plus::Float64; d = 0.005)
    nodes = [d, 2d, 3d, 4d]
    w0 = fd_weights(nodes, 0)    # extrapolate value to 0
    w1 = fd_weights(nodes, 1)    # derivative at 0
    nn = length(nodes)
    K = length(pts0)
    X = Matrix{Float64}(undef, 3, 2nn * K)
    for (k, x0) in enumerate(pts0)
        for (j, t) in enumerate(nodes)
            X[:, 2nn * (k - 1) + j] .= collect(x0) .+ t .* collect(nrm)
            X[:, 2nn * (k - 1) + nn + j] .= collect(x0) .- t .* collect(nrm)
        end
    end
    phi = eval_phi(res, X)
    jumps = Float64[]; fluxes = Float64[]; phivals = Float64[]; fluxvals = Float64[]
    for k in 1:K
        fp = phi[(2nn * (k - 1) + 1):(2nn * (k - 1) + nn)]
        fm = phi[(2nn * (k - 1) + nn + 1):(2nn * k)]
        phi_p = sum(w0 .* fp); phi_m = sum(w0 .* fm)
        dn_p = sum(w1 .* fp)                 # d/dn from + side
        dn_m = -sum(w1 .* fm)                # nodes mirrored -> flip sign
        push!(jumps, abs(phi_p - phi_m)); push!(phivals, abs(phi_m))
        push!(fluxes, abs(eps_plus * dn_p - eps_minus * dn_m)); push!(fluxvals, abs(eps_minus * dn_m))
    end
    # normalize by face-wide flux scale (pointwise |eps dn phi| can vanish by symmetry)
    return (; jump_max = maximum(jumps ./ phivals),
            flux_max = maximum(fluxes) / maximum(fluxvals),
            d = d, floor_est = max(d^3, res.eps / d))
end

println(">>> 6.1.6 transmission checks")
pts12 = [(0.0, y, z) for (y, z) in zip(range(-0.35, 0.35; length = 8), range(0.08, 0.42; length = 8))]
t12 = transmission(res, pts12, (1.0, 0.0, 0.0), 4.0, 12.0)          # Omega1 (x<0) -> Omega2 (x>0)
ptsm = [(x, y, 0.8) for (x, y) in zip(range(-0.22, 0.22; length = 8), range(-0.2, 0.2; length = 8))]
tm = transmission(res, ptsm, (0.0, 0.0, 1.0), 2.0, 1.0)             # Omega_m top: inside -> vacuum
@printf("Omega1-Omega2: max|[phi]|/|phi| = %.3e   max|[eps dnphi]|/|.| = %.3e  (check floor ~%.0e)\n",
        t12.jump_max, t12.flux_max, t12.floor_est)
@printf("dOmega_m top : max|[phi]|/|phi| = %.3e   max|[eps dnphi]|/|.| = %.3e  (check floor ~%.0e)\n",
        tm.jump_max, tm.flux_max, tm.floor_est)

# --------------------------------------------- 6.1.6 sigma neutrality + Gauss
println(">>> sigma body sums + Gauss law")
offsets = cumsum(vcat(0, [BI.num_points(p) for p in panels]))
panel_charge = [sum(panels[i].weights .* res.sigma[(offsets[i] + 1):(offsets[i + 1])]) for i in eachindex(panels)]

function body_charge(bi::Int)
    b = SYS.boxes[bi]; half = (b.Lx / 2, b.Ly / 2, b.Lz / 2)
    Q = 0.0; npan = 0
    for i in eachindex(panels)
        c = centers[i]; rel = c .- b.center
        on = false
        for ax in 1:3
            if abs(abs(rel[ax]) - half[ax]) < 1e-10 && abs(panels[i].normal[ax]) > 0.999 &&
               all(abs(rel[a]) <= half[a] + 1e-10 for a in 1:3 if a != ax)
                on = true
            end
        end
        on || continue
        outward = sum(panels[i].normal .* rel) > 0 ? 1.0 : -1.0
        Q += outward * panel_charge[i]; npan += 1
    end
    return Q, npan
end
# Expected: bodies without enclosed source -> 0; Omega_m (contains source) ->
# bound charge q_b = 1 - int rho/eps(x) dx; total far-field monopole
# int rho_scr + sum sigma = q_free = 1.
S_scr = sum(res.screened_vs.weights .* res.screened_vs.density)
expected = (0.0, 0.0, 1.0 - S_scr)
for (bi, lbl) in zip(1:3, ("Omega1", "Omega2", "Omega_m"))
    Q, npan = body_charge(bi)
    @printf("  body %-7s : sum sigma dS = %+.4e  expected %+.4e  dev %.2e  (%d panels)\n",
            lbl, Q, expected[bi], abs(Q - expected[bi]), npan)
end
@printf("  monopole int rho_scr + sum sigma = %.6f (expected 1)\n", S_scr + sum(panel_charge))

# Gauss law: flux of grad phi through a box of side 4 (vacuum), GL grid/face, central FD
function gauss_flux(res; L = 4.0, npp = 10, d = 0.01)
    gl_x, gl_w = gausslegendre(npp)
    pts = NTuple{3, Float64}[]; nrms = NTuple{3, Float64}[]; ws = Float64[]
    for ax in 1:3, sgn in (-1.0, 1.0)
        n = ntuple(a -> a == ax ? sgn : 0.0, 3)
        for i in 1:npp, j in 1:npp
            u, v = gl_x[i] * L / 2, gl_x[j] * L / 2
            others = setdiff(1:3, ax)
            pt = zeros(3); pt[ax] = sgn * L / 2; pt[others[1]] = u; pt[others[2]] = v
            push!(pts, (pt[1], pt[2], pt[3])); push!(nrms, n)
            push!(ws, gl_w[i] * gl_w[j] * (L / 2)^2)
        end
    end
    K = length(pts)
    X = Matrix{Float64}(undef, 3, 2K)
    for k in 1:K
        X[:, 2k - 1] .= collect(pts[k]) .+ d .* collect(nrms[k])
        X[:, 2k] .= collect(pts[k]) .- d .* collect(nrms[k])
    end
    phi = eval_phi(res, X)
    return sum(ws[k] * (phi[2k - 1] - phi[2k]) / (2d) for k in 1:K)
end
flux = gauss_flux(res)
@printf("  Gauss flux (box side 4): %.6f  expected -1 (free charge)  dev %.2e\n", flux, abs(flux + 1))

# ----------------------------------------------------------- 6.1.5 qualitative
println(">>> 6.1.5 panelization + sigma line data")
lvl = [round(Int, log2(maximum(lens) / l)) for l in lens]
sigma_mean = [mean(res.sigma[(offsets[i] + 1):(offsets[i + 1])]) for i in eachindex(panels)]
pan_rows = [(corners = panels[i].corners, normal = panels[i].normal, level = lvl[i],
             sigma = sigma_mean[i], eps_in = iface.eps_in[i], eps_out = iface.eps_out[i])
            for i in eachindex(panels)]

# sigma along the substrate top surface approaching the junction edge x=0,z=0.5
line = NamedTuple[]
for i in eachindex(panels)
    abs(panels[i].normal[3]) > 0.999 && abs(centers[i][3] - 0.5) < 1e-10 || continue
    for (j, q) in enumerate(panels[i].points)
        if abs(q[2]) < 0.05
            push!(line, (x = q[1], y = q[2], sigma = res.sigma[offsets[i] + j]))
        end
    end
end
save_ref(joinpath(DATA, "raw", "checks_qual.jls"), (;
    eps = EPS, p = P, r = R, N = res.N, niter = res.niter,
    shared_count = length(shared_ids), gammas, dup_count,
    t12, tm, flux, S_scr,
    body_charges = [body_charge(b)[1] for b in 1:3], total_charge = sum(panel_charge),
    panels = pan_rows, sigma_line = line))
println("CHECKS+QUAL DONE")
