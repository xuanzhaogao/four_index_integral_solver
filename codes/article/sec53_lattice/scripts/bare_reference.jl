# Bare (vacuum, unscreened) density-density reference U^bare_ij = V^bare_iijj for the §5.3
# lattice, the denominator of the effective-screening ratio ε_eff(i,j) = U^bare_ij / U_ij
# plotted by scripts/plot_junction_analysis.jl.
#
# U^bare depends only on the pair's TRANSLATION CLASS (template_i, template_j, steps_j - steps_i):
# every orbital is a template translated by integer grid steps (centers.tsv), so the bare integral
# is the same for every pair in a class. Per template, the four extreme orbitals (min/max of x±y)
# as sources, each against every orbital as target, cover every class in the window: a pair's
# displacement fits from the corner on the opposite side of its quadrant. The central orbital is
# added as a cross-check (its classes are also reached from a corner; duplicates must agree).
#
# Same quadrature as codes/graphene/monolayer/src/MonolayerBareIntegrals.jl: source charges
# w·φ_i² on the template grid, TKM3D.ltkm3dc (1/(4πr) kernel) at the grid points of φ_j², eV via
# raw · 4π · E2 / (‖φ_i‖² ‖φ_j‖²) -- the same orbital normalization assemble_v uses for V_eV.
#
# Run: julia --project=<env with BoundaryIntegral> scripts/bare_reference.jl [campaign.toml]
# Writes figs/bare_Uij_<campaign>.tsv: ti tj dsx dsy dsz r_A U_bare_eV
using BoundaryIntegral, Printf, LinearAlgebra
const BI = BoundaryIntegral
using BoundaryIntegral.TKM3D: ltkm3dc

const E2 = 14.3996                      # e^2/(4πε0), eV·Å
const RHO_RTOL = 1e-6                   # keep φ² ≥ RHO_RTOL·max(φ²): tighter than the campaign's 1e-4
const VOLUME_TOL = 1e-6

toml = length(ARGS) >= 1 ? ARGS[1] :
       joinpath(@__DIR__, "..", "campaigns", "lattice_conv_l3_eps2.4_k46.toml")
c = BI.load_campaign(toml)
temps = BI.load_templates!(c)
grids = [t[2] for t in temps]
dg = grids[1]
centers = BI.read_centers(joinpath(c.root, "centers.tsv"))
norb = length(centers)

# Truncated density support of one template: local grid indices and φ² values.
function template_support(g)
    v2 = g.values .^ 2
    m = maximum(v2)
    idx = [I for I in CartesianIndices(v2) if v2[I] >= RHO_RTOL * m]
    return idx, [v2[I] for I in idx]
end
supp = [template_support(g) for g in grids]
At, Bt, Ct = BI.true_cell_vectors(dg)
w = abs(det(hcat(At, Bt, Ct))) / (dg.nx * dg.ny * dg.nz)
h = cbrt(w)
@printf("template grid %d×%d×%d, w = %.3e Å^3, support sizes %s\n",
        dg.nx, dg.ny, dg.nz, w, string([length(s[1]) for s in supp]))

# orbital norm ‖φ‖² = Σ w φ² over the FULL template (not the truncated support)
nrm2 = [w * sum(abs2, g.values) for g in grids]

function positions_of(ct)
    idx, _ = supp[ct.template_id]; s = ct.steps
    P = Matrix{Float64}(undef, 3, length(idx))
    for (k, I) in enumerate(idx)
        p = BI.grid_point(dg, I[1] + s[1], I[2] + s[2], I[3] + s[3])
        P[1, k] = p[1]; P[2, k] = p[2]; P[3, k] = p[3]
    end
    return P
end

cxy = [(ct.center[1], ct.center[2]) for ct in centers]
mx = sum(first.(cxy)) / norb; my = sum(last.(cxy)) / norb
rows = Dict{Tuple{Int,Int,NTuple{3,Int}},Tuple{Float64,Float64}}()
maxdup = 0.0
for t in 1:length(grids), pick in (:center, :sw, :se, :nw, :ne)
    cand = [k for k in 1:norb if centers[k].template_id == t]
    score(k) = pick === :center ? hypot(cxy[k][1] - mx, cxy[k][2] - my) :
               pick === :sw ?  cxy[k][1] + cxy[k][2] : pick === :ne ? -cxy[k][1] - cxy[k][2] :
               pick === :se ? -cxy[k][1] + cxy[k][2] : cxy[k][1] - cxy[k][2]
    i0 = cand[argmin(score.(cand))]
    ci = centers[i0]
    src = positions_of(ci)
    q = w .* supp[t][2]
    js = collect(1:norb)
    tgt = [positions_of(centers[j]) for j in js]
    off = cumsum([0; [size(P, 2) for P in tgt]])
    T = reduce(hcat, tgt)
    t0 = time()
    res = ltkm3dc(VOLUME_TOL, src; charges = q, targets = T, pgt = 1,
                  kmax = BI._estimate_tkm3dc_kmax(h))
    res.ier == 0 || error("ltkm3dc ier=$(res.ier)")
    pot = real.(res.pottarg)
    @printf("source %d (template %d, %s): %d src pts, %d targets, %.1f s\n",
            i0, t, pick, size(src, 2), size(T, 2), time() - t0)
    for (k, j) in enumerate(js)
        cj = centers[j]
        raw = dot(view(pot, off[k]+1:off[k+1]), w .* supp[cj.template_id][2])
        U = raw * 4π * E2 / (nrm2[t] * nrm2[cj.template_id])
        r = hypot(cj.center[1] - ci.center[1], cj.center[2] - ci.center[2])
        key = (t, cj.template_id, cj.steps .- ci.steps)
        if haskey(rows, key)
            global maxdup = max(maxdup, abs(rows[key][2] - U) / U)
        else
            rows[key] = (r, U)
        end
    end
end
@printf("%d classes; max rel. disagreement between duplicate classes: %.2e\n", length(rows), maxdup)
rows = [(k[1], k[2], k[3], v[1], v[2]) for (k, v) in rows]

out = joinpath(@__DIR__, "..", "figs", "bare_Uij_$(c.name).tsv")
open(out, "w") do io
    println(io, "ti\ttj\tdsx\tdsy\tdsz\tr_A\tU_bare_eV")
    for (ti, tj, ds, r, U) in sort(rows; by = x -> (x[1], x[4]))
        @printf(io, "%d\t%d\t%d\t%d\t%d\t%.6f\t%.8f\n", ti, tj, ds..., r, U)
    end
end
println("wrote ", out)
for (ti, tj, ds, r, U) in sort(rows; by = x -> x[4])[1:12]
    @printf("  t%d->t%d r=%6.3f  U_bare=%.4f eV   E2/r=%s\n", ti, tj, r, U,
            r > 0 ? @sprintf("%.4f", E2 / r) : "-")
end
