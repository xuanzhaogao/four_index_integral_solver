# §5.3 junction analysis (1×2): how crossing the buried Si|SiO₂ junction changes the on-site,
# short-range and long-range density-density interactions U_ij = V_iijj, quantitatively rather
# than through the colour maps of fig_63. Both panels plot the deviation from a homogeneous
# ε_slab medium with the true orbital shapes, U - U^bare/ε_slab (zero = no substrate effect):
#   (a) along x: on-site U_ii vs Δx_i and NN / NNN / 3NN U_ij vs the pair midpoint;
#   (b) vs r_ij for three source orbitals -- mid-Si (Δx ≈ -5.25), the central junction orbital,
#       mid-SiO₂ (Δx ≈ +5.25), "mid" = middle of that half of the window -- each against all
#       its partners and itself (r = 0).
# (a) shows per-lattice-column averages as scatterlines; (b) shows every matrix element as a small
# marker with the average (per neighbour shell below 6 Å, per 1-Å bin beyond) as a dashed line.
# U^bare is the orbital-resolved vacuum integral from scripts/bare_reference.jl (not E2/r, which
# is 15% off at the NN distance).
# Δx is measured from the junction plane x = XJ; Si (ε=11.9) at Δx < 0, SiO₂ (ε=3.9) at Δx > 0.
#
# Run: julia --project=<env with CairoMakie> scripts/plot_junction_analysis.jl
# Writes figs/fig_junction_<camp>.{pdf,png} and figs/junction_<camp>_{onsite,shortrange,longrange}.tsv
# (the TSVs and printed numbers also cover the ε_eff = U^bare/U analysis of earlier versions).
using Serialization, Printf, Statistics, DelimitedFiles, CairoMakie, LaTeXStrings
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

# The resolved eps=2.4 tensor (K = 46, max_order = 64), i.e. fig_63's "eps2.4_mo64" variant.
const CAMP = get(ENV, "JUNCTION_CAMP", "lattice_conv_l3_eps2.4_k46")
const E2 = 14.3996                      # e^2/(4πε0), eV·Å
const XJ = 5.547                        # junction plane (Si | SiO₂ cube interface)
const EPS_SLAB = 2.4
const EPS_SI, EPS_SIO2 = 11.9, 3.9
const DEEP = 8.0                        # "deep" = |Δx| ≥ DEEP (the window ends at |Δx| ≈ 10.5)

root = joinpath(CEPH, CAMP)
d = open(deserialize, joinpath(root, "V_full_eV.jls"))
C = readdlm(joinpath(root, "centers.tsv"), '\t'; skipstart = 1)
norb = size(C, 1)
tid = Int.(C[:, 2]); st = [Tuple(Int.(C[k, 3:5])) for k in 1:norb]
px = Float64.(C[:, 6]); py = Float64.(C[:, 7]); dx = px .- XJ

drow = zeros(Int, norb)
for (a, p) in enumerate(d.pair_ids); p[1] == p[2] && (drow[p[1]] = a); end
all(>(0), drow) || error("every orbital needs its diagonal pair (i,i)")
U = d.V[drow, drow]                                    # U[i,j] = V_iijj (eV)
r = [hypot(px[i] - px[j], py[i] - py[j]) for i in 1:norb, j in 1:norb]

B = readdlm(joinpath(@__DIR__, "..", "figs", "bare_Uij_$(CAMP).tsv"), '\t'; skipstart = 1)
bare = Dict((Int(B[k, 1]), Int(B[k, 2]), (Int(B[k, 3]), Int(B[k, 4]), Int(B[k, 5]))) => Float64(B[k, 7])
            for k in 1:size(B, 1))
# U^bare_ij = U^bare_ji, so class (ti, tj, d) is also stored as (tj, ti, -d)
bare_of(i, j) = get(bare, (tid[i], tid[j], st[j] .- st[i]),
                    get(bare, (tid[j], tid[i], st[i] .- st[j]), NaN))
Ub = [bare_of(i, j) for i in 1:norb, j in 1:norb]
# The corner sources of bare_reference.jl miss a few classes spanning the whole window (r > 20 Å,
# not in any plotted series). There the exact tail has converged to E2/r: r·U^bare/E2 ≥ 0.998
# already at r = 12-14 Å, so E2/r is good to < 0.2%. Anything shorter must be exact.
let miss = isnan.(Ub)
    any(miss) && minimum(r[miss]) < 20 && error("bare reference misses classes at r < 20 Å; rerun bare_reference.jl")
    Ub[miss] .= E2 ./ r[miss]
    @printf("bare: %d of %d pairs exact, %d (r = %.1f..%.1f Å) from E2/r\n", count(!, miss), norb^2,
            count(miss), any(miss) ? minimum(r[miss]) : 0.0, any(miss) ? maximum(r[miss]) : 0.0)
end
epsf = Ub ./ U                                          # ε_eff(i,j)

Uon = [U[i, i] for i in 1:norb]

# ---- (a) a zigzag chain crossing the junction: the two adjacent rows nearest mid-window ----
ys = sort(unique(round.(py; digits = 2)))
ymid = (minimum(py) + maximum(py)) / 2
k0 = argmin(abs.(ys .- ymid))
pairrow = abs(ys[k0+1] - ys[k0]) < abs(ys[k0] - ys[k0-1]) ? (ys[k0], ys[k0+1]) : (ys[k0-1], ys[k0])
chain = [k for k in 1:norb if round(py[k]; digits = 2) in pairrow]
chain = chain[sortperm(dx[chain])]

# ---- (b) shells: pairs (i<j) at the NN, NNN and 3NN distances, binned by midpoint Δx ----
shells = (("NN", 1.4225), ("NNN", 2.465), ("3NN", 2.845))
binx(v) = round(v / 0.25) * 0.25      # midpoints sit on a ~0.62 Å lattice; 0.25 Å bins separate them
function binned(xs, vs)
    bx = binx.(xs); ub = sort(unique(bx))
    return ub, [mean(vs[bx .== b]) for b in ub], [length(vs[bx .== b]) for b in ub]
end
shell_data = Dict{String,Any}()
for (lab, s) in shells
    ps = [(i, j) for i in 1:norb for j in i+1:norb if abs(r[i, j] - s) < 0.05]
    xm = [(dx[i] + dx[j]) / 2 for (i, j) in ps]; v = [U[i, j] for (i, j) in ps]
    ub = [Ub[i, j] for (i, j) in ps]
    shell_data[lab] = (; ps, xm, v, ub)
end

# ---- (c),(d) representative sources ----
pick(x0) = argmin([hypot(dx[k] - x0, py[k] - ymid) for k in 1:norb])
iSi = pick(-maximum(abs.(dx))); iOx = pick(maximum(abs.(dx))); iJ = pick(0.0)
side_si(j) = dx[j] < 0; side_ox(j) = dx[j] > 0
# (c): the middle of each half of the window, and the junction orbital
xhalf = maximum(abs.(dx)) / 2
csrc = [(lab = "mid-Si source",   i = pick(-xhalf), col = QUAL.blue,   mk = QUAL_MK.blue),
        (lab = "junction source", i = iJ,           col = QUAL.purple, mk = QUAL_MK.purple),
        (lab = "mid-SiO₂ source", i = pick(xhalf),  col = QUAL.red,    mk = QUAL_MK.red)]
series = [
    (lab = "Si source, Si partners",       i = iSi, js = [j for j in 1:norb if j != iSi && side_si(j)], col = QUAL.blue,   mk = QUAL_MK.blue),
    (lab = "SiO₂ source, SiO₂ partners",   i = iOx, js = [j for j in 1:norb if j != iOx && side_ox(j)], col = QUAL.red,    mk = QUAL_MK.red),
    (lab = "junction source → Si side",    i = iJ,  js = [j for j in 1:norb if j != iJ && side_si(j)],  col = QUAL.green,  mk = QUAL_MK.green),
    (lab = "junction source → SiO₂ side",  i = iJ,  js = [j for j in 1:norb if j != iJ && side_ox(j)],  col = QUAL.orange, mk = QUAL_MK.orange),
]

# ---------------------------------------------------------------------------- TSV outputs
figs = joinpath(@__DIR__, "..", "figs")
open(joinpath(figs, "junction_$(CAMP)_onsite.tsv"), "w") do io
    println(io, "orbital\ttemplate\tx_A\ty_A\tdx_A\tU_eV\tU_bare_eV\teps_eff\tin_chain")
    for k in sortperm(dx)
        @printf(io, "%d\t%d\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.5f\t%d\n",
                k, tid[k], px[k], py[k], dx[k], Uon[k], Ub[k, k], epsf[k, k], k in chain)
    end
end
open(joinpath(figs, "junction_$(CAMP)_shortrange.tsv"), "w") do io
    println(io, "shell\ti\tj\tr_A\tdx_i_A\tdx_j_A\tdx_mid_A\tU_ij_eV\tU_bare_eV\teps_eff")
    for (lab, _) in shells, (k, (i, j)) in enumerate(shell_data[lab].ps)
        @printf(io, "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.5f\n", lab, i, j, r[i, j],
                dx[i], dx[j], shell_data[lab].xm[k], U[i, j], Ub[i, j], epsf[i, j])
    end
end
open(joinpath(figs, "junction_$(CAMP)_longrange.tsv"), "w") do io
    println(io, "series\ti\tj\tdx_i_A\tdx_j_A\tr_A\tU_ij_eV\tU_bare_eV\teps_eff\tr_U_over_E2")
    for s in series, j in s.js
        @printf(io, "%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.6f\t%.6f\t%.5f\t%.5f\n", s.lab, s.i, j,
                dx[s.i], dx[j], r[s.i, j], U[s.i, j], Ub[s.i, j], epsf[s.i, j], r[s.i, j] * U[s.i, j] / E2)
    end
end

# ---------------------------------------------------------------------------- figure
# Paper style (codes/fig_gen/fig_style.jl, cf. fig5_plot.jl): FIG_W × FIG_H, 1×2, panel tags inside
# the top-left corner, colour + marker locked per series. (a) draws the per-column averages as
# scatterlines; (b) draws every matrix element as a small opaque marker and the average as a
# dashed line.
const MS_PT = 5                 # (b) raw-value marker size
const A_PT  = 0.5               # (b) raw-value marker opacity
lab_tex(s) = latexstring(s)

fig = Figure(size = (FIG_W, FIG_H))

# (a) U - U^bare/ε_slab along x: on-site vs Δx_i, NN / NNN / 3NN vs the pair midpoint;
#     average = mean over each lattice column (0.25 Å bin of Δx)
ax_a = Axis(fig[1, 1];
            xlabel = L"\Delta x\ (\text{\AA})",
            ylabel = L"U - U^{\mathrm{bare}}/\varepsilon_{\mathrm{slab}}\ (\mathrm{eV})",
            xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
            yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
hlines!(ax_a, [0.0]; color = (:black, 0.6), linestyle = :dot, linewidth = LW_GUIDE)
vlines!(ax_a, [0.0]; color = (:black, 0.6), linestyle = :dash, linewidth = LW_GUIDE)
bseries = vcat([(L"\text{on‐site}",dx, Uon, [Ub[k, k] for k in 1:norb], LINE_COLORS[1], LINE_MARKERS[1], "on-site")],
               [(lab_tex("\\text{$lab}"),
                 shell_data[lab].xm, shell_data[lab].v, shell_data[lab].ub,
                 LINE_COLORS[k+1], LINE_MARKERS[k+1], lab) for (k, (lab, d0)) in enumerate(shells)])
bdev = Dict{String,Tuple{Float64,Float64}}()
for (lab, xs, vs, ubs, col, mk, key) in bseries
    du = vs .- ubs ./ EPS_SLAB
    bx, bv, _ = binned(xs, du)
    scatterlines!(ax_a, bx, bv; color = col, marker = mk, markersize = MS, linewidth = LW_DATA, label = lab)
    bdev[key] = (mean(du[xs .<= -DEEP]), mean(du[xs .>= DEEP]))
end
axislegend(ax_a; position = :rb)
xlims!(ax_a, -11.5, 11.5); ylims!(ax_a, -0.225, 0.1)

# (b) U_ij - U^bare_ij/ε_slab vs r_ij for three source orbitals and all their partners (r = 0 is
#     the on-site term); average = mean over each neighbour shell below 6 Å, 1-Å bins beyond
ax_b = Axis(fig[1, 2];
            xlabel = L"r_{ij}\ (\text{\AA})",
            ylabel = L"U_{ij} - U^{\mathrm{bare}}_{ij}/\varepsilon_{\mathrm{slab}}\ (\mathrm{eV})",
            xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
            yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5))
hlines!(ax_b, [0.0]; color = (:black, 0.6), linestyle = :dot, linewidth = LW_GUIDE)
rbin(rr) = [x < 6 ? round(x; digits = 1) : round(x) for x in rr]
function shellmean(rr, vv)
    b = rbin(rr); ub = sort(unique(b))
    return [mean(rr[b .== x]) for x in ub], [mean(vv[b .== x]) for x in ub]
end
# legend: Δx of each source (centre) orbital
csty = [(QUAL.blue, QUAL_MK.blue, :solid), (QUAL.purple, QUAL_MK.purple, :dash), (QUAL.red, QUAL_MK.red, :dashdot)]
for (s, (col, mk, ls)) in zip(csrc, csty)
    lab = latexstring(@sprintf("\\Delta x = {%s}%.1f\\ \\text{\\AA}", dx[s.i] < 0 ? "-" : "+", abs(dx[s.i])))
    rr = r[s.i, :]; du = U[s.i, :] .- Ub[s.i, :] ./ EPS_SLAB      # includes j = i (r = 0)
    scatter!(ax_b, rr, du; color = (col, A_PT), marker = mk, markersize = MS_PT, label = lab)
    mr, mu = shellmean(rr, du)
    lines!(ax_b, mr, mu; color = col, linestyle = ls, linewidth = LW_DATA, label = lab)
end
axislegend(ax_b; position = :rb, merge = true)
xlims!(ax_b, -0.8, 20.5); ylims!(ax_b, -0.215, 0.03)

colgap!(fig.layout, 1, 30)
for (ax, lab) in ((ax_a, "(a)"), (ax_b, "(b)"))
    text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
          offset = (6, -6), font = :bold, fontsize = FS_BASE)
end

out = joinpath(figs, "fig_junction_$(CAMP).pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
save(replace(out, ".pdf" => ".png"), fig; px_per_unit = 2)
println("wrote ", out)

# ---------------------------------------------------------------------------- numbers for the text
deep(x) = (x .<= -DEEP, x .>= DEEP)
let (si, ox) = deep(dx)
    @printf("on-site U: deep Si %.3f eV, deep SiO2 %.3f eV (%+.1f%%); range %.3f .. %.3f; chain y = %s\n",
            mean(Uon[si]), mean(Uon[ox]), 100 * (mean(Uon[ox]) / mean(Uon[si]) - 1),
            minimum(Uon), maximum(Uon), string(pairrow))
    @printf("           eps_eff on-site: Si %.3f, SiO2 %.3f (U_bare %.3f eV)\n",
            mean([epsf[k, k] for k in findall(si)]),mean([epsf[k, k] for k in findall(ox)]),
            mean([Ub[k, k] for k in 1:norb]))
end
for (lab, s) in shells
    sd = shell_data[lab]; si, ox = deep(sd.xm)
    @printf("%-3s r=%.2f: n=%d  deep Si %.3f eV, deep SiO2 %.3f eV (%+.1f%%); bare %.3f eV; eps_eff Si %.2f SiO2 %.2f\n",
            lab, s, length(sd.ps), mean(sd.v[si]), mean(sd.v[ox]), 100 * (mean(sd.v[ox]) / mean(sd.v[si]) - 1),
            mean(sd.ub), mean(sd.ub[si] ./ sd.v[si]), mean(sd.ub[ox] ./ sd.v[ox]))
end
for s in series
    rr = r[s.i, s.js]; ee = epsf[s.i, s.js]; uu = U[s.i, s.js]
    @printf("%-28s i=%d dx=%+.2f: eps_eff(0)=%.2f", s.lab, s.i, dx[s.i], epsf[s.i, s.i])
    for rb in (1.42, 2.46, 5.0, 8.0, 12.0, 15.0)
        k = abs.(rr .- rb) .< 0.6
        any(k) && @printf("  r~%g: eps %.2f rU/E2 %.3f", rb, mean(ee[k]), mean(rr[k] .* uu[k]) / E2)
    end
    println()
end
for s in csrc
    js = [j for j in 1:norb if j != s.i]; rr = r[s.i, js]; uu = U[s.i, js]
    @printf("(c) %-16s i=%d dx=%+.2f y=%.2f: U_ii=%.3f", s.lab, s.i, dx[s.i], py[s.i], U[s.i, s.i])
    @printf("  dU(0) %+.3f", U[s.i, s.i] - Ub[s.i, s.i] / EPS_SLAB)
    ub = Ub[s.i, js]
    for rb in (1.42, 2.46, 5.0, 10.0, 15.0, 20.0)
        k = abs.(rr .- rb) .< 0.6
        any(k) && @printf("  r~%g: dU %+.3f", rb, mean(uu[k] .- ub[k] ./ EPS_SLAB))
    end
    @printf("  (max r %.1f)\n", maximum(rr))
end
for (k, v) in bdev
    @printf("(b) %-18s U - U^bare/eps_slab: deep Si %+.3f, deep SiO2 %+.3f eV\n", k, v...)
end
