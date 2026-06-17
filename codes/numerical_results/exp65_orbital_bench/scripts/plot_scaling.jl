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
xlab = L"K\;\;(\text{right-hand sides / pair densities})"

begin
    fig = Figure(size = (880, 350))

    ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = K, title = "(a) runtime vs K", xscale = log10, yscale = log10)
    scatterlines!(ax1, K, t_pre ./ K; color = C[1], marker = :circle, label = "precompute")
    scatterlines!(ax1, K, t_blk ./ K; color = C[2], marker = :circle, label = "block solve")
    scatterlines!(ax1, K, t_evl ./ K; color = C[3], marker = :circle, label = "eval (onsite row)")
    scatterlines!(ax1, K, t_pot ./ K; color = C[5], marker = :circle, label = "eval pottrg build")
    # scatterlines!(ax1, K, t_tot; color = :black, marker = :rect, linewidth = 2.5,
    #             label = "total (precompute+block+eval)")
    axislegend(ax1; position = :lt, framevisible = false, labelsize = 9)

    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "peak RSS (GB)", xticks = K, title = "(b) peak memory vs K")
    scatterlines!(ax2, K, rss; color = C[4], marker = :circle)
    length(K) > 1 && hlines!(ax2, [rss[1]]; color = :gray70, linestyle = :dot)

    fig
end

save(joinpath(FIGS, FIGNAME * ".pdf"), fig)
save(joinpath(FIGS, FIGNAME * ".png"), fig; px_per_unit = 2)
println("gathered ", length(K), " cutoffs (K = ", K, ") -> ", joinpath(DATA, "multi_rhs.csv"),
        " + figs/", FIGNAME, ".{pdf,png}")
