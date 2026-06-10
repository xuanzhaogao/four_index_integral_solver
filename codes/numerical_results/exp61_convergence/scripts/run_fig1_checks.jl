# Correctness checks for the Fig.-1 system (slab on two cubes), 6.1.6/6.1.7 style.
# One solve (default eps = 1e-4, p = 6, r = 4; override via CHK_EPS/CHK_P/CHK_R),
# then:
#   - shared-face census + duplicate-panel scan (geometry sanity)
#   - transmission [phi] and [eps dn phi] across: Omega1|Omega2 (x = 0),
#     slab|Omega1 contact (z = 0, x < 0), slab|Omega2 contact (z = 0, x > 0),
#     slab top (z = 1)
#   - sum(sigma dS) per body (expected: 0, 0, 1 - 1/eps_slab) and Gauss's law
#     on an enclosing box (expected flux -1)
# Results: console + data/raw/fig1_checks.jls
#
# Run:  julia --project=. exp61_convergence/scripts/run_fig1_checks.jl

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf, Statistics
using FastGaussQuadrature: gausslegendre

const SYS = system_fig1()
const DATA = joinpath(@__DIR__, "..", "data")
const EPS = parse(Float64, get(ENV, "CHK_EPS", "1e-4"))
const P = parse(Int, get(ENV, "CHK_P", "6"))
const R = parse(Int, get(ENV, "CHK_R", "4"))

println(">>> warm-up"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1)

println(">>> solve Fig.1 system  eps=$(EPS) p=$(P) r=$(R)"); flush(stdout)
res = solve_system(SYS; eps = EPS, p = P, r = R)
@printf("N=%d niter=%d residual=%.2e adaptive_pairs=%d\n",
        res.N, res.niter, res.residual, res.n_adaptive_pairs)

iface = res.interface
panels = iface.panels
centers = [(p.corners[1] .+ p.corners[2] .+ p.corners[3] .+ p.corners[4]) ./ 4 for p in panels]

# ---------------------------------------------------------------- geometry sanity
println(">>> geometry sanity (census + duplicates)")
pairkey(i) = (round(min(iface.eps_in[i], iface.eps_out[i]); digits = 6),
              round(max(iface.eps_in[i], iface.eps_out[i]); digits = 6))
census = Dict{Tuple{Float64, Float64}, Int}()
for i in eachindex(panels)
    census[pairkey(i)] = get(census, pairkey(i), 0) + 1
end
for (k, v) in sort(collect(census))
    @printf("  eps pair %-14s %6d panels\n", string(k), v)
end
dup_count = let seen = Dict{NTuple{3, Float64}, Int}(), cnt = 0
    for c in centers
        key = (round(c[1]; digits = 9), round(c[2]; digits = 9), round(c[3]; digits = 9))
        cnt += (seen[key] = get(seen, key, 0) + 1) > 1 ? 1 : 0
    end
    cnt
end
println("  duplicate-center panels: ", dup_count, " (expected 0)")

# ---------------------------------------------------------------- transmission
function fd_weights(nodes::Vector{Float64}, deriv::Int)
    n = length(nodes)
    Vm = [nodes[j]^k for k in 0:(n - 1), j in 1:n]
    rhs = zeros(n); rhs[deriv + 1] = factorial(deriv)
    return Vm \ rhs
end

"4-node one-sided FD transmission check across a planar interface."
function transmission(res, pts0, nrm::NTuple{3, Float64},
                      eps_minus::Float64, eps_plus::Float64; d = 0.005)
    nodes = [d, 2d, 3d, 4d]
    w0 = fd_weights(nodes, 0); w1 = fd_weights(nodes, 1)
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
        push!(jumps, abs(sum(w0 .* fp) - sum(w0 .* fm)))
        push!(phivals, abs(sum(w0 .* fm)))
        dn_p = sum(w1 .* fp); dn_m = -sum(w1 .* fm)
        push!(fluxes, abs(eps_plus * dn_p - eps_minus * dn_m))
        push!(fluxvals, abs(eps_minus * dn_m))
    end
    return (; jump_max = maximum(jumps ./ phivals),
            flux_max = maximum(fluxes) / maximum(fluxvals),
            d = d, floor_est = max(d^3, res.eps / d))
end

println(">>> transmission checks (floor ~", max(0.005^3, EPS / 0.005), ")")
checks = [
    ("Omega1|Omega2 (x=0)", [(0.0, y, z) for (y, z) in zip(range(-4.0, 4.0; length = 8),
                                                           range(-8.5, -1.0; length = 8))],
     (1.0, 0.0, 0.0), 4.0, 12.0),
    ("slab|Omega1  (z=0)", [(x, y, 0.0) for (x, y) in zip(range(-4.2, -0.8; length = 8),
                                                          range(-3.5, 3.5; length = 8))],
     (0.0, 0.0, 1.0), 4.0, 10.0),
    ("slab|Omega2  (z=0)", [(x, y, 0.0) for (x, y) in zip(range(0.8, 4.2; length = 8),
                                                          range(-3.5, 3.5; length = 8))],
     (0.0, 0.0, 1.0), 12.0, 10.0),
    ("slab|vac     (z=1)", [(x, y, 1.0) for (x, y) in zip(range(-4.2, 4.2; length = 8),
                                                          range(-3.5, 3.5; length = 8))],
     (0.0, 0.0, 1.0), 10.0, 1.0),
]
tres = NamedTuple[]
for (lbl, pts, nrm, em, ep) in checks
    t = transmission(res, pts, nrm, em, ep)
    push!(tres, (; lbl, t.jump_max, t.flux_max))
    @printf("  %-20s max|[phi]|/|phi| = %.3e   max|[eps dnphi]|/scale = %.3e\n",
            lbl, t.jump_max, t.flux_max)
    flush(stdout)
end

# ------------------------------------------------- sigma body sums + Gauss law
println(">>> sigma body sums + Gauss law")
offsets = cumsum(vcat(0, [BI.num_points(p) for p in panels]))
panel_charge = [sum(panels[i].weights .* res.sigma[(offsets[i] + 1):(offsets[i + 1])])
                for i in eachindex(panels)]

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
# NOTE: bodies share faces here, so bound charge on a shared interface is not
# attributable to a single body -> no per-body invariant (unlike an isolated
# body in vacuum). Per-body sums are reported as diagnostics only; the
# validated conservation targets are the global monopole and Gauss's law.
S_scr = sum(res.screened_vs.weights .* res.screened_vs.density)
body_Q = Float64[]
for (bi, lbl) in zip(1:3, ("Omega1", "Omega2", "slab"))
    Q, npan = body_charge(bi)
    push!(body_Q, Q)
    @printf("  body %-7s : sum sigma dS = %+.4e  (diagnostic; %d panels)\n", lbl, Q, npan)
end
monopole = S_scr + sum(panel_charge)
@printf("  monopole int rho_scr + sum sigma = %.6f  expected 1  dev %.2e\n",
        monopole, abs(monopole - 1))

"Gauss flux through a box (center, side L), GL grid per face, central FD."
function gauss_flux(res; center = (0.0, 0.0, -4.5), L = 23.0, npp = 10, d = 0.01)
    gl_x, gl_w = gausslegendre(npp)
    pts = NTuple{3, Float64}[]; nrms = NTuple{3, Float64}[]; ws = Float64[]
    for ax in 1:3, sgn in (-1.0, 1.0)
        n = ntuple(a -> a == ax ? sgn : 0.0, 3)
        for i in 1:npp, j in 1:npp
            u, v = gl_x[i] * L / 2, gl_x[j] * L / 2
            others = setdiff(1:3, ax)
            pt = collect(center); pt[ax] += sgn * L / 2
            pt[others[1]] += u; pt[others[2]] += v
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
@printf("  Gauss flux (box side 23 around system): %.6f  expected -1  dev %.2e\n",
        flux, abs(flux + 1))

save_ref(joinpath(DATA, "raw", "fig1_checks.jls"), (;
    eps = EPS, p = P, r = R, N = res.N, niter = res.niter,
    census = collect(census), dup_count, transmission = tres,
    body_Q, S_scr, total_charge = sum(panel_charge), monopole, flux))
println("FIG1 CHECKS DONE")
