# §6.3 figure (2×2), style matching numerical_results/fig_gen (fig_style.jl):
#   (a) |V[ρ_ij,ρ_kl]| heatmap (eV, log), l_ec = 1.14 (level 3);
#   (b) on-site U top-down map (Δx to the buried junction), l_ec = 1.14;
#   (c) single-orbital on-site U convergence vs number of unknowns N (l_ec sweep), with the
#       Richardson-extrapolated U∞ as a dashed guide;
#   (d) non-local density-density interaction U_ij = V_iijj vs |r_i - r_j|, for a representative
#       i on each side of the buried junction (Si, ε=11.9 vs SiO₂, ε=3.9).  The two sides separate
#       because the environment screening differs across the junction; the dashed 1/r guide shows
#       that these terms stay long-ranged (cf. the local Hund's / pair-hopping terms V_ijij).
# Panel tags (a)/(b)/(c)/(d) sit OUTSIDE, at the top-left corner of each panel.
# Run: julia --project=codes/lattice_scale scripts/plot_fig_63.jl
using BoundaryIntegral, Serialization, Printf, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))   # sets theme: fontsize=FS_BASE
const CEPH = "/mnt/ceph/users/xgao1/four_index"

# Data variants, selected with FIG63_VARIANT (default "v2", the published figure).
# The three inputs can no longer share a single suffix: the eps=2.4 rerun writes V to its OWN
# campaign root rather than to a lattice_conv_l3* tree, and it has no l_ec convergence sweep of
# its own, so each input is named explicitly.
const VARIANT = get(ENV, "FIG63_VARIANT", "v2")
const SRC = Dict(
    ""       => (camp = "lattice_10x10_het3x",        vroot = "lattice_conv_l3",
                 tsv  = "lec_single_conv.tsv",        out   = "fig_63.pdf",       conv_eps = 10.0),
    "v2"     => (camp = "lattice_10x10_het3x_v2",     vroot = "lattice_conv_l3_v2",
                 tsv  = "lec_single_conv_v2.tsv",     out   = "fig_63_v2.pdf",    conv_eps = 10.0),
    # level-3, all solver tolerances 1e-5 -- the like-for-like counterpart of lattice_conv_l3_v2.
    # (The earlier lattice_10x10_het3x_eps2.4 tree is level 2 with rhs_tol=1e-3; superseded.)
    "eps2.4" => (camp = "lattice_conv_l3_eps2.4",     vroot = "lattice_conv_l3_eps2.4",
                 tsv  = "lec_single_conv_eps2.4.tsv", out   = "fig_63_eps2.4.pdf", conv_eps = 2.4),
)[VARIANT]
const SLAB_EPS = VARIANT == "eps2.4" ? 2.4 : 10.0
if SLAB_EPS != SRC.conv_eps
    @warn("panel (c) convergence series is from a DIFFERENT slab permittivity",
          figure_eps = SLAB_EPS, convergence_eps = SRC.conv_eps, tsv = SRC.tsv)
end
c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", SRC.camp * ".toml"))
d3 = open(deserialize, joinpath(CEPH, SRC.vroot, "V_full_eV.jls"))   # l_ec = 1.14 (ref)
p3 = sortperm(d3.pair_ids); V3 = d3.V[p3, p3]; n = size(V3, 1)

# on-site U (diagonal pairs) from level 3
px = [o.pos[1] for o in c.orbitals]; py = [o.pos[2] for o in c.orbitals]; norb = length(c.orbitals)
U = fill(NaN, norb); for (a, p) in enumerate(d3.pair_ids); p[1] == p[2] && (U[p[1]] = d3.V[a, a]); end
XJ = 5.547

# single-orbital l_ec convergence: U vs dof
Uc = Float64[]; Nc = Float64[]
# Panel (d) uses refinement levels >= 2 (an open upper end, so levels 6/7 appear as soon as
# they are computed). Level 1 (l_ec = 4.55 A) sits outside the asymptotic
# range: at eps=2.4 it is non-monotonic with level 2 (U rises before it falls), so including it
# fits a power law through a point the trend does not govern.
const CONV_MIN_LEVEL = 2
for ln in eachline(joinpath(@__DIR__, "..", "figs", SRC.tsv))
    (isempty(ln) || startswith(ln, "level")) && continue
    f = split(ln, '\t')
    parse(Int, f[1]) >= CONV_MIN_LEVEL || continue
    push!(Uc, parse(Float64, f[3])); push!(Nc, parse(Float64, f[4]))
end
q = sortperm(Nc); UU = Uc[q]; NN = Nc[q]
dd = [UU[k] - UU[k+1] for k in 1:length(UU)-1]; r = dd[end] / dd[end-1]
Uinf = UU[end] - dd[end] * r / (1 - r)     # Richardson (first order in l_ec)

begin
    # Uses the shared FIG_W with the usual 0.99\textwidth include. On-page label size is
    # 6900 * f / canvas_units pt, so it is the RATIO canvas/f that must be held at ~1000 to match
    # the other figures, not the product: 1000 at 0.99 gives 6.83 pt, as do the FIG_W figures.
    # (800 at 0.8 also works, at a smaller printed figure; 640 at 0.99 would give 10.7 pt.)
    #
    # Height, not gap size, is what closes the axis->colourbar distance: every panel is square
    # (DataAspect, or aspect=1 for (d)), so too short a figure leaves the axis boxes narrower
    # than their columns and the slack shows up between axis and colourbar.
    fig = Figure(size = (FIG_W, 820), fontsize = FS_BASE)

    # ---- (a) |V| heatmap ----
    ga = fig[1, 1] = GridLayout()
    Aflo = 1e-6; A = abs.(V3); LA = log10.(clamp.(A, Aflo, Inf))
    axa = Axis(ga[1, 1]; aspect = DataAspect(),
            yreversed = true, xticks = [], yticks = [], xlabel = L"(i, j)", ylabel = L"(k, l)")
    hma = heatmap!(axa, 1:n, 1:n, LA; colormap = Reverse(:viridis), colorrange = (log10(Aflo), maximum(LA)), rasterize = 4)
    Colorbar(ga[1, 2], hma; label = L"\log_{10}|V|\ (\mathrm{eV})", width = 8)

    # ---- (b) top-down on-site U ----
    gb = fig[1, 2] = GridLayout()
    axb = Axis(gb[1, 1]; xlabel = L"\Delta x\ (\text{\AA})", ylabel = L"y\ (\text{\AA})", aspect = DataAspect())
    vlines!(axb, [0.0]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    sc = scatter!(axb, px .- XJ, py; color = U, colormap = :viridis, markersize = 11,
                colorrange = (minimum(U), maximum(U)))
    Colorbar(gb[1, 2], sc; label = L"U\ (\mathrm{eV})", width = 8)

    # ---- (c) non-local U_ij = V_iijj from a central atom, mapped like (b) ----
    gc = fig[2, 1] = GridLayout()
    # row of the diagonal pair (i,i) for each orbital i; 0 if that pair was not retained
    drow = zeros(Int, norb)
    for (a, p) in enumerate(d3.pair_ids); p[1] == p[2] && (drow[p[1]] = a); end
    # i is fixed at the atom nearest the centre of the lattice
    cx = sum(px) / norb; cy = sum(py) / norb
    i0 = argmin([drow[k] == 0 ? Inf : hypot(px[k] - cx, py[k] - cy) for k in 1:norb])
    js = [j for j in 1:norb if j != i0 && drow[j] > 0]
    Uij = [abs(d3.V[drow[i0], drow[j]]) for j in js]
    axc = Axis(gc[1, 1]; xlabel = L"\Delta x\ (\text{\AA})", ylabel = L"y\ (\text{\AA})",
            aspect = DataAspect())
    vlines!(axc, [0.0]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    scm = scatter!(axc, px[js] .- XJ, py[js]; color = Uij, colormap = :viridis,
            colorrange = (minimum(Uij), maximum(Uij)), markersize = 11)
    scatter!(axc, [px[i0] - XJ], [py[i0]]; color = :red, markersize = 11)   # the fixed orbital i
    Colorbar(gc[1, 2], scm; label = L"U_{ij}\ (\mathrm{eV})", width = 8)

    # ---- (d) on-site U convergence vs DOF (relative error, log-log) ----
    gd = fig[2, 2] = GridLayout()
    Er = abs.(UU .- Uinf) ./ abs(Uinf)     # relative error vs Richardson U∞
    lx = log10.(NN); ly = log10.(Er)       # least-squares power-law fit E_r ~ N^b
    b = (length(lx) * sum(lx .* ly) - sum(lx) * sum(ly)) / (length(lx) * sum(lx .^ 2) - sum(lx)^2)
    a = (sum(ly) - b * sum(lx)) / length(lx)
    # Limits data-derived, NOT hard-coded: E_r and N shift by orders of magnitude between
    # variants (edge correction, permittivity, level set), and a fixed window renders (d) empty
    # or crowds the points into a corner. Margins are in decades.
    xlo = 10^(log10(minimum(NN)) - 0.22); xhi = 10^(log10(maximum(NN)) + 0.22)
    ylo = 10^(log10(minimum(Er)) - 0.25); yhi = 10^(log10(maximum(Er)) + 0.25)
    # The fit spans the full axis rather than just the fitted points, so the guide reads as a
    # trend line across the panel.
    xf = [xlo, xhi]; yf = 10 .^ (a .+ b .* log10.(xf))
    # Round decade/half-decade ticks chosen from the actual range: Makie's default on a
    # sub-decade log window falls back to 10^5.75, 10^6.00, ... which is unreadable at this size.
    # One formatter for both axes, valid for any magnitude (an earlier version returned an empty
    # label for v < 1, which would have silently blanked the E_r axis).
    # Built with rich()/superscript(), NOT LaTeXStrings: L"..." ticks render through the math
    # engine in a serif math font and do not match the theme font used by the tick labels of
    # panels (b) and (c).
    function ticklab(v)
        e = Int(floor(log10(v) + 1e-9)); m = Int(round(v / 10.0^e))
        m == 1 ? rich("10", superscript(string(e))) : rich("$(m)×10", superscript(string(e)))
    end
    # Prefer whole decades; only fall back to the 1/2/5 sequence when the window is too narrow
    # to show at least two of them (which is what forces the 5x10^5 / 2x10^6 labels).
    function tickvals(lo, hi)
        dec = [10.0^e for e in ceil(Int, log10(lo)):floor(Int, log10(hi))]
        length(dec) >= 2 && return dec
        return [m * 10.0^e for e in floor(Int, log10(lo)):ceil(Int, log10(hi))
                for m in (1, 2, 5) if lo <= m * 10.0^e <= hi]
    end
    xtv = tickvals(xlo, xhi)
    ytv = tickvals(ylo, yhi)
    axd = Axis(gd[1, 1]; xscale = log10, yscale = log10, xlabel = L"N", ylabel = L"\mathcal{E}_r",
            aspect = 1,                      # square box, matching the DataAspect panels
            xticks = (xtv, ticklab.(xtv)), yticks = (ytv, ticklab.(ytv)),
            xminorticksvisible = true, yminorticksvisible = true,
            xminorticks = IntervalsBetween(6), yminorticks = IntervalsBetween(6),
            xminorgridvisible = true, yminorgridvisible = true)
    lines!(axd, xf, yf; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    scatter!(axd, NN, Er; color = QUAL.blue, marker = :circle, markersize = MS)
    xlims!(axd, xlo, xhi); ylims!(axd, ylo, yhi)

    # ---- panel tags OUTSIDE, top-left corner of each panel ----
    for (g, lab) in ((ga, "(a)"), (gb, "(b)"), (gc, "(c)"), (gd, "(d)"))
        Label(g[1, 1, TopLeft()], lab; font = :bold, fontsize = FS_BASE,
            halign = :right, padding = (0, 8, 4, 0))
    end
    # tighten each panel against its own colourbar, then the panels against each other.
    # gd has a single column (no colourbar), so it has no column gap to set.
    for g in (ga, gb, gc); colgap!(g, 5); end
    colgap!(fig.layout, 6); rowgap!(fig.layout, 6)

    fig
end

out = joinpath(@__DIR__, "..", "figs", SRC.out)
save(out, fig; px_per_unit = PX_PER_UNIT)
@printf("wrote %s   max|V|=%.3f eV, U∞=%.4f eV (r=%.3f)\n", out, maximum(abs.(V3)), Uinf, r)
# panel (c) numbers for the text: the fixed central orbital and the U_ij range it sees
let rr = [hypot(px[j] - px[i0], py[j] - py[i0]) for j in js],
    lx = log10.(rr), ly = log10.(Uij), m = length(rr)
    sd = (m * sum(lx .* ly) - sum(lx) * sum(ly)) / (m * sum(lx .^ 2) - sum(lx)^2)
    @printf("  (c) i=%d at (x=%.2f, y=%.2f) A, dx_from_junction=%+.2f A\n",
            i0, px[i0], py[i0], px[i0] - XJ)
    @printf("      U_ii=%.3f eV;  U_ij over %d neighbours: %.3f .. %.3f eV;  U_ij ~ r^%.2f\n",
            d3.V[drow[i0], drow[i0]], length(js), minimum(Uij), maximum(Uij), sd)
end
