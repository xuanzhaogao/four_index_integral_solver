# 6.6 multi-RHS figure, redesigned as a 2-panel strong-scaling / amortization study:
#   (a) K=1 thread scaling: total (precompute + block solve) runtime vs nthreads in
#       {1,2,4,8,16,...,96}, with an ideal-linear-speedup reference and per-point
#       speedup-vs-1-thread labels.
#   (b) nthreads=96 multi-RHS scaling: PER-RHS runtime (total/K) vs K (REUSES the existing
#       K-sweep records unchanged), with a flat single-RHS-cost reference and per-K speedup
#       labels (identical ratios to a total-vs-K view, since both series are divided by K).
# "total" in both panels is precompute + block solve only — the eval stage is a small,
# roughly K-independent tail (a few seconds to ~1 minute) that isn't the point of either
# scaling story, so it is dropped from the headline metric (still on disk if needed later).
# Renders figs/fig66_multicube_scaling.pdf.
#
# Panel (a) reconstructs the full block-solve time at low thread counts from a TIMING-ONLY
# itmax=1 run (data/raw_threads/multicube_nt<N>.jls, written by run_multicube.jl with
# MULTICUBE_GMRES_ITMAX=1) as t_RHS_assembly + t_1iter * niter_ref: RHS assembly is a
# one-time cost, and only the per-iteration block-GMRES cost (t_1iter) is scaled by
# niter_ref, the GMRES iteration count from the real converged K=1 run at 96 threads
# (data/raw/multicube_K1_cut0.0.jls). GMRES's iteration count is a property of the operator,
# not the thread count, so one reference transfers across the whole sweep. This avoids paying
# for an hours-long converged solve at nthreads=1. (Scaling the LUMPED t_solve_block, which
# folds RHS assembly into the single iteration, by niter_ref would over-count RHS assembly
# ~niter_ref times and depress the high-thread points — a bug fixed here.)
#
# Run AFTER both Slurm arrays finish (run_multicube_array.sbatch for (b),
# run_multicube_threads_array.sbatch for (a)):
#   julia --project=codes/numerical_results plot_scripts/plot_scaling.jl
#
# Records are plain NamedTuples of scalars/vectors (no BoundaryIntegral dependency).

using CairoMakie, LaTeXStrings, Serialization
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const SMOKE = get(ENV, "MULTICUBE_SMOKE", "0") == "1"
# RERUN_TAG: plot the rerun tree (data_v2/, *_v2 campaigns) instead of the published one.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "exp66_multicube", "data" * TAG, SMOKE ? "smoke" : "")
const RAW  = joinpath(DATA, "raw")
const RAW_THREADS = joinpath(DATA, "raw_threads")
const FIGS = joinpath(@__DIR__, "..", "figs" * TAG); mkpath(FIGS)
const FIGNAME = SMOKE ? "fig66_multicube_scaling_smoke" : "fig66_multicube_scaling"
mkpath(FIGS)

# ============================================================================
# Panel (a) data: K=1 thread sweep (timing-only itmax=1 records)
# ============================================================================
# Prefer the pinned single-node sweep (multicube_pin_nt<N>.jls, run_multicube_threads_pinned
# .sbatch) — placement-clean; fall back to the older unpinned per-node array records if the
# pinned set is absent.
thread_files = let files = isdir(RAW_THREADS) ? readdir(RAW_THREADS) : String[]
    pinned = sort(filter(f -> occursin(r"^multicube_pin_nt\d+\.jls$", f), files))
    isempty(pinned) ? sort(filter(f -> occursin(r"^multicube_nt\d+\.jls$", f), files)) : pinned
end
isempty(thread_files) && error("no thread-sweep records in $RAW_THREADS — run " *
    "run_multicube_threads_pinned.sbatch first")
@info "thread sweep source" nfiles=length(thread_files) pinned=any(occursin("pin_", f) for f in thread_files)

thread_recs = [deserialize(joinpath(RAW_THREADS, f)) for f in thread_files]
sort!(thread_recs; by = r -> r.nthreads)
any(r -> r.gmres_itmax > 1, thread_recs) &&
    @warn "a thread-sweep record has gmres_itmax > 1 (not timing-only) — check RUNTAG files"

# reference GMRES iteration count from the real converged K=1 run at 96 threads
k1_ref_file = joinpath(RAW, "multicube_K1_cut0.0.jls")
isfile(k1_ref_file) || error("missing reference $k1_ref_file — run run_multicube_array.sbatch first")
k1_ref = deserialize(k1_ref_file)
const NITER_REF = k1_ref.niter
@info "thread-sweep reconstruction reference" niter_ref=NITER_REF k1_ref.nthreads k1_ref.t_total

stage_t(r, lbl) = first(s.t for s in r.stages if s.label == lbl)

nthreads  = [r.nthreads for r in thread_recs]
pre_th    = [r.t_precompute for r in thread_recs]
# Reconstruct the full block solve as (one-time RHS assembly) + niter_ref * (per-iteration
# block-GMRES cost). NB: r.t_solve_block LUMPS the RHS assembly with the single GMRES
# iteration, so scaling that whole number by niter_ref would count the RHS assembly niter_ref
# times — instead pull the two stages apart and multiply only the per-iteration GMRES cost.
rhs_th    = [stage_t(r, "RHS assembly (per-source, multi-region)") for r in thread_recs]
gmres1_th = [stage_t(r, "block GMRES") for r in thread_recs]   # itmax=1 -> one iteration
blk_th    = rhs_th .+ gmres1_th .* NITER_REF                    # reconstructed full block-solve time
tot_th    = pre_th .+ blk_th                                    # reconstructed total runtime (eval excluded)

# ============================================================================
# Panel (b) data: nthreads=96 K sweep (unchanged from the previous panel (a))
# ============================================================================
files = sort(filter(f -> occursin(r"^multicube_K\d+_cut.*\.jls$", f), readdir(RAW)))
isempty(files) && error("no per-cutoff records in $RAW — run the Slurm array first")

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
t_tot = t_pre .+ t_blk                    # total (eval excluded, see header note)

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

# per-RHS (amortized) runtime and its naive baseline (paying the K=1 cost every time,
# i.e. no batching benefit at all — a flat line at t_tot[1])
per_rhs = t_tot ./ K
setup_per_rhs = t_pre ./ K          # system setup amortized over K RHS (built once -> ~1/K)
block_per_rhs = t_blk ./ K          # block solve amortized over K RHS (batched-FMM -> ~4x floor)
naive_per_rhs = fill(t_tot[1], length(K))
speedup_K = t_tot[1] ./ per_rhs

# ============================================================================
# Two-panel scaling summary, in the shared fig style (fig_gen/fig_style.jl):
#   (a) strong scaling : total runtime vs thread count  (log2 x, log10 y)
#   (b) amortization   : runtime per RHS vs K at 96 threads (linear y)
# Light-gray secondary references, one annotation per panel, minimal legend.
# ============================================================================
const REF_CLR  = (:gray50, 0.9)      # light gray, visually secondary references
const BAND_CLR = (:steelblue, 0.10)  # very light optimal-K band

sp96  = round(tot_th[1] / tot_th[end]; digits = 1)   # strong-scaling speedup at 96 threads
spmax = round(maximum(speedup_K); digits = 1)        # best per-RHS speedup (near K = 20-25)

begin
    fig = Figure(size = (FIG_W, FIG_H))

    # ----- Panel (a): strong scaling -------------------------------------
    ax1 = Axis(fig[1, 1];
        xlabel = L"N_\mathrm{threads}", ylabel = "total runtime (s)",
        xscale = log2, yscale = log10,
        xticks = (nthreads, string.(nthreads)),
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(5))

    # ideal linear speedup is a straight line on log-log, anchored at the 1-thread total
    lines!(ax1, [nthreads[1], 132], tot_th[1] ./ [nthreads[1], 132];
        color = :black, linestyle = :dot, linewidth = 2, label = "ideal")
    hlines!(ax1, [tot_th[1]]; color = REF_CLR, linestyle = :solid,
        linewidth = LW_GUIDE, label = "baseline")
    scatterlines!(ax1, nthreads, tot_th; color = QUAL.green, marker = QUAL_MK.green,
        linewidth = LW_DATA, markersize = MS, label = "measured")

    xlims!(ax1, 0.85, 132)
    ylims!(ax1, 20, 6000)

    # single annotation: strong-scaling speedup at the final point.
    # Derived from the data, NOT hard-coded: this number changes whenever
    # niter_ref changes (e.g. when the preconditioner is on).
    text!(ax1, 96, 5e2; text = latexstring(string(round(tot_th[1] / tot_th[end]; digits = 1), " \\times")),
        align = (:right, :bottom), offset = (-4, 8), fontsize = FS_ANNOT + 2, color = :black)

    annotation!(ax1, 96, tot_th[1], 96, tot_th[end];
        # text = L"20.3 \times",
        # path = Ann.Paths.Arc(0.3),
        # style = Ann.Styles.LineArrow(),
        style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(), tail = Ann.Arrows.Head()),
        fontsize = FS_ANNOT + 2,
        color = :black,
        labelspace = :data
    )

    # annotation!(ax1, 96, 100, 0, 900, style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(), tail = Ann.Arrows.Head()))

    # annotation!(ax1, 96, 96, 100, 1000;
            # text = "minimum ≈ $(spmax)×\nK ≈ 20–25", fontsize = FS_ANNOT - 4)

    axislegend(ax1; position = :lb, labelsize = FS_LEGEND - 2, rowgap = 1)

    # ----- Panel (b): multi-RHS amortization -----------------------------
    ax2 = Axis(fig[1, 2];
        xlabel = L"K", ylabel = "runtime per RHS (s)",
        xticks = (K, string.(K)),
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(5))

    # vspan!(ax2, 19, 25; color = BAND_CLR)   # optimal-K band, drawn behind the data
    hlines!(ax2, [naive_per_rhs[1]]; color = REF_CLR, linestyle = :solid, linewidth = LW_GUIDE)
    # text!(ax2, K[end], naive_per_rhs[1]; text = "single-RHS baseline",
    #     align = (:right, :top), offset = (-4, -6), fontsize = FS_ANNOT - 4, color = (:gray, 0.9))
    scatterlines!(ax2, K, per_rhs; color = QUAL.green, marker = QUAL_MK.green,
        linewidth = LW_DATA, markersize = MS)

    # Derived, not hard-coded: the K sweep was extended to K = 46 (cutoff 6.5) and a fixed
    # upper limit of 33 clipped every point past K = 31 -- including the amortization turnover.
    xlims!(ax2, 0, maximum(K) + 2)
    # Derived, not pinned. The old fixed 0-250 window was sized when the sweep stopped at
    # K = 31 and the data topped out at 91 s; it left the panel less than half full, which
    # mattered once the sweep reached K = 46: the amortization TURNOVER (30.3 -> 36.6 s per
    # RHS beyond K = 37) is a 6.3 s feature and was rendering at 2.5% of the panel height.
    # Anchored at 0 because this axis is a runtime, so the origin is meaningful.
    ylims!(ax2, 0, 1.08 * max(maximum(per_rhs), naive_per_rhs[1]))

    # amortization maximum, derived from the data rather than hard-coded. Positioned off the
    # arrow's midpoint rather than at fixed coordinates -- the old (19.5, 120) with a +16px
    # offset sat above the new axis top and clipped.
    text!(ax2, 22, (minimum(per_rhs) + naive_per_rhs[1]) / 2;
        text = latexstring(string(round(naive_per_rhs[1] / minimum(per_rhs); digits = 1), " \\times")),
        align = (:right, :center), offset = (-8, 0), fontsize = FS_ANNOT + 2, color = :black)

    annotation!(ax2, 22, per_rhs[argmin(per_rhs)], 22, naive_per_rhs[1];
        # text = L"3.6 \times",
        # path = Ann.Paths.Arc(0.3),
        # style = Ann.Styles.LineArrow(),
        style = Ann.Styles.LineArrow(head = Ann.Arrows.Head(), tail = Ann.Arrows.Head()),
        fontsize = FS_ANNOT + 2,
        color = :black,
        labelspace = :data
    )

    # ----- panel tags -----------------------------------------------------
    for (ax, lab) in ((ax1, "(a)"), (ax2, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
            offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    colgap!(fig.layout, 1, 30)
    fig
end

save(joinpath(FIGS, FIGNAME * ".pdf"), fig; px_per_unit = PX_PER_UNIT)
println("thread sweep: nthreads = $nthreads")
println("K sweep: K = $K")
println("wrote $(joinpath(FIGS, FIGNAME * ".pdf"))")
