# 6.5 figure: multi-RHS runtime + peak RAM vs K (number of right-hand sides / pair
# densities), production env (96 threads) with the in-function FINUFFT nthreads cap (16).
# (a) per-phase runtime (precompute / block solve / eval) + production total vs K, with
#     the OMP=96 pathological total as a dashed reference; (b) peak RSS vs K.
# Output: figs/fig65_multirhs_scaling.{pdf,png}
#
# Run: julia --project=codes/numerical_results exp65_orbital_bench/scripts/plot_scaling.jl

using CairoMakie, LaTeXStrings

const DATA = joinpath(@__DIR__, "..", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

# minimal CSV reader -> NamedTuple of columns (data are plain numbers; no deps on BI)
function readcsv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = String.(split(lines[1], ","))
    idx = Dict(h => i for (i, h) in enumerate(header))
    rows = [String.(split(l, ",")) for l in lines[2:end]]
    col(name) = [parse(Float64, r[idx[name]]) for r in rows]
    p = sortperm(col("K"))
    return (; K = Int.(col("K"))[p], t_pre = col("t_precompute")[p],
            t_blk = col("t_solve_block")[p], t_seq = col("t_solve_seq")[p],
            t_eval = col("t_eval")[p], rss = col("rss_peak_gb")[p])
end

fix = readcsv(joinpath(DATA, "multi_rhs.csv"))   # in-function FINUFFT cap (16)
t_total_fix = fix.t_pre .+ fix.t_blk .+ fix.t_eval

const C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]  # Okabe-Ito
xlab = L"K\;\;(\text{right-hand sides / pair densities})"

fig = Figure(size = (880, 350))

ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = fix.K,
    title = "(a) runtime vs K  (96 threads, FINUFFT→16 in-function)")
scatterlines!(ax1, fix.K, fix.t_pre;     color = C[1], marker = :circle, label = "precompute")
scatterlines!(ax1, fix.K, fix.t_blk;     color = C[2], marker = :circle, label = "block solve")
scatterlines!(ax1, fix.K, fix.t_eval;    color = C[3], marker = :circle, label = "eval (4-index)")
scatterlines!(ax1, fix.K, t_total_fix;   color = :black, marker = :rect, linewidth = 2.5,
              label = "total (precompute+block+eval)")
axislegend(ax1; position = :lt, framevisible = false, labelsize = 9)

ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "peak RSS (GB)", xticks = fix.K,
    title = "(b) peak memory vs K")
scatterlines!(ax2, fix.K, fix.rss; color = C[4], marker = :circle)
hlines!(ax2, [fix.rss[1]]; color = :gray70, linestyle = :dot)
text!(ax2, fix.K[2], fix.rss[1] + 4; text = "≈18 GB baseline + ~7.5 GB/RHS\n(block-GMRES nd=K FMM scratch)",
      fontsize = 9, color = :gray40, align = (:left, :bottom))

save(joinpath(FIGS, "fig65_multirhs_scaling.pdf"), fig)
save(joinpath(FIGS, "fig65_multirhs_scaling.png"), fig; px_per_unit = 2)
println("saved -> ", joinpath(FIGS, "fig65_multirhs_scaling.{pdf,png}"))
