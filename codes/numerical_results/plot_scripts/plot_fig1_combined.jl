# Combined Fig.-1 figure, merging the system geometry and the convergence study
# into one 3-panel figure:
#   (a) system geometry (left, spanning both rows; 3D axis hidden)
#   (b) relative V-error vs N        (top right)
#   (c) GMRES iteration count vs N   (bottom right)
# Supersedes visualize_fig1_system.jl + plot_fig1.jl for the paper.
# Output: figs/fig61_fig1_combined.pdf
#
# The geometry is the RHS-adaptive interface at the sweep parameters, an
# expensive BoundaryIntegral build; it is cached to data/cache/fig1_viz_iface.jls
# so re-rendering the figure (layout tweaks) does not rebuild it. Delete that
# file to force a rebuild.

include(joinpath(@__DIR__, "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using CairoMakie, LaTeXStrings, Printf, Serialization
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

# RERUN_TAG: plot the rerun tree (data_v2/, *_v2 campaigns) instead of the published one.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "exp61_convergence", "data" * TAG)
const FIGS = joinpath(@__DIR__, "..", "figs" * TAG); mkpath(FIGS)
mkpath(FIGS)

# ---- geometry: RHS-adaptive interface at the sweep parameters --------------
const SYS = system_fig1()
const VIZ_EPS, VIZ_P, VIZ_R = 1e-4, 6, 2
const GEOM_COLOR = RGBAf(0.40, 0.52, 0.68, 1.0)
const CACHE = joinpath(DATA, "cache"); mkpath(CACHE)
const IFACE_CACHE = joinpath(CACHE, "fig1_viz_iface.jls")

iface = if isfile(IFACE_CACHE)
    @info "loading cached interface" IFACE_CACHE
    deserialize(IFACE_CACHE)
else
    @info "building RHS-adaptive interface (heavy; will be cached)"
    vs  = gaussian_source(SYS.src_center, SYS.src_sigma, VIZ_EPS)
    svs = Harness.screened_source(SYS, vs)
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
    ifc = BI.multi_dielectric_box3d_rhs_adaptive(VIZ_P, Harness.l_ec_of(SYS, VIZ_R),
        SYS.boxes, SYS.epses, svs, 1.0, VIZ_EPS, SYS.eps_out;
        max_depth = 128, tkm_kmax = kmax)
    try
        serialize(IFACE_CACHE, ifc)
    catch e
        @warn "interface cache write failed; will rebuild next time" exception = e
    end
    ifc
end
println("geometry panels: ", length(iface.panels))

# ---- convergence data ------------------------------------------------------
ref = deserialize(joinpath(DATA, "raw", "fig1_ref.jls"))
load_series(p) = [deserialize(joinpath(DATA, "raw", "fig1_p$(p)_r$(r).jls")) for r in 1:5]
const NOEDGE = joinpath(DATA, "raw", "fig1_noedges")
load_noedge(p) = [deserialize(joinpath(NOEDGE, "fig1_p$(p)_r$(r).jls"))
                  for r in 1:5 if isfile(joinpath(NOEDGE, "fig1_p$(p)_r$(r).jls"))]

const _SC = sweep_colors(3)
const _SM = sweep_markers(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])
const MK  = Dict(2 => _SM[1], 4 => _SM[2], 6 => _SM[3])

const MExt = Base.get_extension(BI, :MakieExt)

begin
    fig = Figure(size = (FIG_W, 380))

    # ---- (a) geometry: left panel, axis removed ---------------------------
    ax_g = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.72π, elevation = 0.16π,
        protrusions = 0)
    MExt.viz_3d!(ax_g, iface; show_points = false, highlight_edges = false,
        base_color = GEOM_COLOR)
    scatter!(ax_g, [Point3f(SYS.src_center...)]; color = :red, markersize = 14)
    hidedecorations!(ax_g)
    hidespines!(ax_g)

    # ---- (b) relative V-error vs N: middle --------------------------------
    ax1 = Axis(fig[1, 2]; xscale = log10, yscale = log10,
        xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
        yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
        xlabel = L"N", ylabel = L"\mathcal{E}_{r}")
    # ---- (c) GMRES iterations vs N: right ---------------------------------
    ax2 = Axis(fig[1, 3]; xscale = log10,
        xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
        yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
        xlabel = L"N", ylabel = L"N_\mathrm{iter}")

    for p in (2, 4, 6)
        s = load_series(p)
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        scatterlines!(ax1, Ns, ev; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
        scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS)
    end
    # no-edge-correction contrast: same color + marker per p, dashed (faded)
    for p in (2, 4)
        s = load_noedge(p)
        isempty(s) && continue
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = MK[p],
            markersize = MS, linewidth = LW_GUIDE, linestyle = :dash,
            label = L"p = %$p \text{ (no edge corr.)}")
        scatterlines!(ax2, Ns, [d.niter for d in s]; color = (COL[p], 0.55),
            marker = MK[p], markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
    end
    axislegend(ax1; position = :lb, labelsize = FS_LEGEND - 3, nbanks = 1,
        rowgap = 0, padding = (5, 5, 3, 3))
    ylims!(ax1, 10^(-4.5), 10^(-1.5))
    ylims!(ax2, 15, 30)

    # panel tags inside each panel (in-axis text), matching the other figures
    for (ax, lab) in ((ax_g, "(a)"), (ax1, "(b)"), (ax2, "(c)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
            offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    colgap!(fig.layout, 10)

    save(joinpath(FIGS, "fig61_fig1_combined.pdf"), fig; px_per_unit = PX_PER_UNIT)
    println("wrote ", joinpath(FIGS, "fig61_fig1_combined.pdf"))
    fig
end
