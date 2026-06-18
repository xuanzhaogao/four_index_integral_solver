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

const SMOKE = get(ENV, "MULTICUBE_SMOKE", "0") == "1"
const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
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

const C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]  # Okabe-Ito
xlab = L"K"

begin
    fig = Figure(size = (1400, 400), fontsize = 18)

    ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = K, xscale = log10, yscale = log10)
    # reference slopes (faint guides)
    lines!(ax1, K, fill(sum(t_pre) / length(t_pre), length(K)); color = :gray80)              # slope 0
    lines!(ax1, K, t_evl[1] .* K;                                color = :gray80, linestyle = :dot)  # slope 1 (∝K)
    # fixed cost as its average level (pottrg omitted — negligible here)
    hlines!(ax1, [sum(t_pre) / length(t_pre)]; color = C[1], linestyle = :dash, linewidth = 2, label = "precompute")
    # scaling costs vs K
    scatterlines!(ax1, K, t_blk; color = C[2], marker = :circle, label = "block solve")
    scatterlines!(ax1, K, t_evl; color = C[3], marker = :circle, label = "eval")
    scatterlines!(ax1, K, t_tot; color = :black, marker = :rect, linewidth = 2, label = "total")
    axislegend(ax1; position = :rb, labelsize = 16)
    ylims!(ax1, 1.0, 10^(3.3))

    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "peak RAM (GB)", xticks = K)
    scatterlines!(ax2, K, rss; color = C[4], marker = :circle)
    length(K) > 1 && hlines!(ax2, [rss[1]]; color = :gray70, linestyle = :dot)
    ylims!(ax2, 0, 350)

    ax3 = Axis(fig[1, 3]; xlabel = xlab, ylabel = "GMRES iterations", xticks = K)
    scatterlines!(ax3, K, niter; color = C[5], marker = :circle)
    ylims!(ax3, 0, 40)

    fig
end

save(joinpath(FIGS, FIGNAME * ".pdf"), fig)
save(joinpath(FIGS, FIGNAME * ".png"), fig)
println("gathered $(length(K)) cutoffs (K = $K) -> data/multicube.csv + figs/$(FIGNAME).{pdf,png}")
