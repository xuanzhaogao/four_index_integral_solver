# §6.3 figure (1×3), style matching numerical_results/fig_gen (fig_style.jl):
#   (a) |V[ρ_ij,ρ_kl]| heatmap (eV, log), l_ec = 1.14 (level 3);
#   (b) on-site U top-down map (Δx to the buried junction), l_ec = 1.14;
#   (c) single-orbital on-site U convergence vs number of unknowns N (l_ec sweep), with the
#       Richardson-extrapolated U∞ as a dashed guide.
# Panel tags (a)/(b)/(c) sit OUTSIDE, at the top-left corner of each panel.
# Run: julia --project=codes/lattice_scale scripts/plot_fig_63.jl
using BoundaryIntegral, Serialization, Printf, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))   # sets theme: fontsize=FS_BASE
const CEPH = "/mnt/ceph/users/xgao1/four_index"

# RERUN_TAG: plot the rerun tree (data_v2/, *_v2 campaigns) instead of the published one.
const TAG = get(ENV, "RERUN_TAG", "")
c = BI.load_campaign(joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x$(TAG).toml"))
d3 = open(deserialize, joinpath(CEPH, "lattice_conv_l3$(TAG)", "V_full_eV.jls"))   # l_ec = 1.14 (ref)
p3 = sortperm(d3.pair_ids); V3 = d3.V[p3, p3]; n = size(V3, 1)

# on-site U (diagonal pairs) from level 3
px = [o.pos[1] for o in c.orbitals]; py = [o.pos[2] for o in c.orbitals]; norb = length(c.orbitals)
U = fill(NaN, norb); for (a, p) in enumerate(d3.pair_ids); p[1] == p[2] && (U[p[1]] = d3.V[a, a]); end
XJ = 5.547

# single-orbital l_ec convergence: U vs dof
Uc = Float64[]; Nc = Float64[]
for ln in eachline(joinpath(@__DIR__, "..", "figs", "lec_single_conv$(TAG).tsv"))
    (isempty(ln) || startswith(ln, "level")) && continue
    f = split(ln, '\t'); push!(Uc, parse(Float64, f[3])); push!(Nc, parse(Float64, f[4]))
end
q = sortperm(Nc); UU = Uc[q]; NN = Nc[q]
dd = [UU[k] - UU[k+1] for k in 1:length(UU)-1]; r = dd[end] / dd[end-1]
Uinf = UU[end] - dd[end] * r / (1 - r)     # Richardson (first order in l_ec)

begin
    fig = Figure(size = (FIG_W, 300), fontsize = FS_BASE)

    # ---- (a) |V| heatmap ----
    ga = fig[1, 1] = GridLayout()
    Aflo = 1e-6; A = abs.(V3); LA = log10.(clamp.(A, Aflo, Inf))
    axa = Axis(ga[1, 1]; aspect = DataAspect(),
            yreversed = true, xticks = [], yticks = [], xlabel = L"(i, j)", ylabel = L"(k, l)")
    hma = heatmap!(axa, 1:n, 1:n, LA; colormap = :viridis, colorrange = (log10(Aflo), maximum(LA)), rasterize = 4)
    Colorbar(ga[1, 2], hma; label = L"\log_{10}|V|\ (\mathrm{eV})", width = 8)

    # ---- (b) top-down on-site U ----
    gb = fig[1, 2] = GridLayout()
    axb = Axis(gb[1, 1]; xlabel = L"\Delta x\ (\text{\AA})", ylabel = L"y\ (\text{\AA})", aspect = DataAspect())
    vlines!(axb, [0.0]; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    sc = scatter!(axb, px .- XJ, py; color = U, colormap = :viridis, markersize = 11,
                colorrange = (minimum(U), maximum(U)))
    Colorbar(gb[1, 2], sc; label = L"U\ (\mathrm{eV})", width = 8)

    # ---- (c) on-site U convergence vs DOF (relative error, log-log) ----
    gc = fig[1, 3] = GridLayout()
    Er = abs.(UU .- Uinf) ./ abs(Uinf)     # relative error vs Richardson U∞
    lx = log10.(NN); ly = log10.(Er)       # least-squares power-law fit E_r ~ N^b
    b = (length(lx) * sum(lx .* ly) - sum(lx) * sum(ly)) / (length(lx) * sum(lx .^ 2) - sum(lx)^2)
    a = (sum(ly) - b * sum(lx)) / length(lx)
    xf = [1e5, 1e7]; yf = 10 .^ (a .+ b .* log10.(xf))
    axc = Axis(gc[1, 1]; xscale = log10, yscale = log10, xlabel = L"N", ylabel = L"\mathcal{E}_r",
            xminorticksvisible = true, yminorticksvisible = true,
            xminorticks = IntervalsBetween(6), yminorticks = IntervalsBetween(6),
            xminorgridvisible = true, yminorgridvisible = true)
    lines!(axc, xf, yf; color = :black, linestyle = :dash, linewidth = LW_GUIDE)
    scatter!(axc, NN, Er; color = QUAL.blue, marker = :circle, markersize = MS)
    # Data-derived, NOT hard-coded: with the edge correction on, E_r drops ~50x
    # and a fixed window sized for the uncorrected errors renders panel (c) empty.
    ylims!(axc, 10^(log10(minimum(Er)) - 0.4), 10^(log10(maximum(Er)) + 0.4));
    # xlims!(axc, 1e5, 1e7)

    # ---- panel tags OUTSIDE, top-left corner of each panel ----
    for (g, lab) in ((ga, "(a)"), (gb, "(b)"), (gc, "(c)"))
        Label(g[1, 1, TopLeft()], lab; font = :bold, fontsize = FS_BASE,
            halign = :right, padding = (0, 8, 4, 0))
    end
    colgap!(fig.layout, 18)

    fig
end

out = joinpath(@__DIR__, "..", "figs", "fig_63$(TAG).pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
@printf("wrote %s   max|V|=%.3f eV, U∞=%.4f eV (r=%.3f)\n", out, maximum(abs.(V3)), Uinf, r)
