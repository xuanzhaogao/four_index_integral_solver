# 6.5 figure: multi-RHS runtime + peak RAM vs K. GATHERS the per-cutoff records written by
# run_multi_rhs.jl (data/raw/multi_rhs_K*_cut*.jls — one per Slurm-array task) into
# data/multi_rhs.csv, then renders figs/fig65_multirhs_scaling.{pdf,png}.
#
# Run AFTER all per-cutoff jobs finish:
#   julia --project=codes/numerical_results exp65_orbital_bench/scripts/plot_scaling.jl
#
# Records are plain NamedTuples of scalars/vectors (no BoundaryIntegral dependency).

using CairoMakie, LaTeXStrings, Serialization

const SMOKE = get(ENV, "ORBBENCH_SMOKE", "0") == "1"
const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
const RAW  = joinpath(DATA, "raw")
const FIGS = joinpath(@__DIR__, "..", "figs")
const FIGNAME = SMOKE ? "fig65_multirhs_scaling_smoke" : "fig65_multirhs_scaling"
mkpath(FIGS)

files = sort(filter(f -> occursin(r"^multi_rhs_K\d+_cut.*\.jls$", f), readdir(RAW)))
isempty(files) && error("no per-cutoff records in $RAW — run the Slurm array (or run_multi_rhs.jl) first")

# one record per K; if a K repeats, keep the newest file (defensive against stale runs)
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
t_tot = t_pre .+ t_pot .+ t_blk .+ t_evl

open(joinpath(DATA, "multi_rhs.csv"), "w") do io
    println(io, join(["hostname", "nthreads", "smoke", "correct_edges", "cutoff", "K",
        "n_points", "n_src", "niter", "block_resid", "t_precompute", "t_pottrg",
        "t_solve_block", "t_eval", "t_total", "rss_peak_gb", "u_onsite_ev"], ","))
    for r in recs
        println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges, r.cutoff,
            r.K, r.n_points, r.n_src, r.niter, r.block_resid, r.t_precompute, r.t_pottrg,
            r.t_solve_block, r.t_eval, r.t_total, r.rss_peak_gb, r.u_onsite_ev]), ","))
    end
end

const C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]  # Okabe-Ito
xlab = L"K"

begin
    fig = Figure(size = (1000, 400), fontsize = 18)

      ax1 = Axis(fig[1,1]; xlabel=xlab, ylabel="runtime (s)", xticks=K, xscale=log10, yscale=log10)
    # reference slopes (faint guides)
    lines!(ax1, K, fill(sum(t_pre)/length(t_pre), length(K)); color=:gray80)          # slope 0
    lines!(ax1, K, t_evl[1] .* K;                              color=:gray80, linestyle=:dot) # slope 1 (∝K)
    # fixed costs as their average level
    hlines!(ax1, [sum(t_pre)/length(t_pre)]; color=C[1], linestyle=:dash, linewidth=2, label="precompute")
    hlines!(ax1, [sum(t_pot)/length(t_pot)]; color=C[5], linestyle=:dash, linewidth=2, label="pottrg build")
    # scaling costs vs K
    scatterlines!(ax1, K, t_blk; color=C[2], marker=:circle, label="block solve")
    scatterlines!(ax1, K, t_evl; color=C[3], marker=:circle, label="eval")
    scatterlines!(ax1, K, t_tot; color=:black, marker=:rect, linewidth=2, label="total")
    axislegend(ax1; position=:rb, labelsize=16)
    ylims!(ax1, 1.0, 10^(2.5))

    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "peak RAM (GB)", xticks = K)
    scatterlines!(ax2, K, rss; color = C[4], marker = :circle)
    length(K) > 1 && hlines!(ax2, [rss[1]]; color = :gray70, linestyle = :dot)
    ylims!(ax2, 0, 180)

    fig
end

save(joinpath(FIGS, FIGNAME * ".pdf"), fig)
