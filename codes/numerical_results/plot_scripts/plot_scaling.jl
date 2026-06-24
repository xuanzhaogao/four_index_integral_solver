# 6.6 figure: multicube multi-RHS runtime + peak RAM vs K. GATHERS the per-cutoff records
# written by run_multicube.jl (data/raw/multicube_K*_cut*.jls — one per Slurm-array task)
# into data/multicube.csv, then renders figs/fig66_multicube_scaling.{pdf,png}.
# Mirrors exp65's fig65_multirhs_scaling, but WITHOUT the pottrg curve (it is ~0.2 s here —
# the thick slab keeps the orbital clear of the interface, so no near-field corrections; see NOTES).
#
# Run AFTER the array finishes:
#   julia --project=codes/numerical_results exp66_multicube/scripts/plot_scaling.jl
#
# Records are plain NamedTuples of scalars/vectors (no BoundaryIntegral dependency).

using CairoMakie, LaTeXStrings, Serialization
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const SMOKE = get(ENV, "MULTICUBE_SMOKE", "0") == "1"
const DATA = joinpath(@__DIR__, "..", "exp66_multicube", "data", SMOKE ? "smoke" : "")
const RAW  = joinpath(DATA, "raw")
const FIGS = joinpath(@__DIR__, "..", "figs")
const FIGNAME = SMOKE ? "fig66_multicube_scaling_smoke" : "fig66_multicube_scaling"
mkpath(FIGS)

files = sort(filter(f -> occursin(r"^multicube_K\d+_cut.*\.jls$", f), readdir(RAW)))
isempty(files) && error("no per-cutoff records in $RAW — run the Slurm array first")

# one record per K; if a K repeats, keep the newest file
byK = Dict{Int, Any}()
for f in files
    r = deserialize(joinpath(RAW, f))
    if !haskey(byK, r.K) || mtime(joinpath(RAW, f)) > byK[r.K][1]
        byK[r.K] = (mtime(joinpath(RAW, f)), r)
    end
end
recs = [byK[k][2] for k in sort(collect(keys(byK)))]

K     = [r.K for r in recs]
t_pre = [r.t_precompute for r in recs]
t_pot = [r.t_pottrg for r in recs]
t_blk = [r.t_solve_block for r in recs]
t_evl = [r.t_eval for r in recs]
rss   = [r.rss_peak_gb for r in recs]
niter = [r.niter for r in recs]
t_tot = t_pre .+ t_pot .+ t_blk .+ t_evl     # pottrg included in the total but not shown separately

open(joinpath(DATA, "multicube.csv"), "w") do io
    println(io, join(["hostname", "nthreads", "smoke", "correct_edges", "cutoff", "K",
        "n_points", "n_src", "niter", "block_resid", "v11_raw", "v11_vac", "screen_ratio",
        "t_precompute", "t_pottrg", "t_solve_block", "t_eval", "t_total", "rss_peak_gb"], ","))
    for r in recs
        println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges, r.cutoff,
            r.K, r.n_points, r.n_src, r.niter, r.block_resid, r.v11_raw, r.v11_vac, r.screen_ratio,
            r.t_precompute, r.t_pottrg, r.t_solve_block, r.t_eval, r.t_total, r.rss_peak_gb]), ","))
    end
end

xlab = L"K"
# multi-RHS speedup vs repeating a single-RHS solve K times (precompute is ~3% of total)
naive = K .* t_tot[1]
speedup = naive ./ t_tot

begin
    fig = Figure(size = (FIG_W, FIG_H))

    # ----- Panel (a): runtime decomposition + multi-RHS speedup -----------
    ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = K,
        xscale = log10, yscale = log10)
    # naive baseline: repeat a single-RHS solve K times (the gap to total = speedup)
    lines!(ax1, K, naive; color = (:gray50, 0.9), linestyle = :dash,
        linewidth = LW_GUIDE, label = L"K \times \mathrm{single\text{-}RHS}")
    # fixed cost as its average level (pottrg omitted — negligible here)
    hlines!(ax1, [sum(t_pre) / length(t_pre)]; color = QUAL.blue, linestyle = :dash,
        linewidth = LW_GUIDE, label = "precompute")
    scatterlines!(ax1, K, t_blk; color = QUAL.orange, marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = "block solve")
    scatterlines!(ax1, K, t_evl; color = QUAL.green, marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = "eval")
    scatterlines!(ax1, K, t_tot; color = :black, marker = :rect,
        linewidth = LW_DATA, markersize = MS, label = "total")
    # speedup labels at each K >= 4 (skip the K=1 baseline)
    for i in eachindex(K)
        K[i] == 1 && continue
        text!(ax1, K[i], t_tot[i]; text = string(round(speedup[i]; digits = 1), "×"),
            align = (:center, :top), offset = (0, -6), fontsize = FS_ANNOT - 4,
            color = :black)
    end
    axislegend(ax1; position = :rb)
    ylims!(ax1, 1.0, 10^(3.7))

    # ----- Panel (b): GMRES iterations vs K -------------------------------
    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "GMRES iterations", xticks = K)
    scatterlines!(ax2, K, niter; color = QUAL.purple, marker = :circle,
        linewidth = LW_DATA, markersize = MS)
    ylims!(ax2, 0, 40)

    for (ax, lab) in ((ax1, "(a)"), (ax2, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
              offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    fig
end

save(joinpath(FIGS, FIGNAME * ".pdf"), fig; px_per_unit = PX_PER_UNIT)
println("gathered $(length(K)) cutoffs (K = $K) -> data/multicube.csv + figs/$(FIGNAME).pdf")
