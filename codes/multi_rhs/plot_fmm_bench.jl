# Plot per-RHS FMM3D kernel time from the bench_fmm.jl sweep (data: figs/ + data//fmm_bench.csv).
# Two panels: (a) per-RHS FMM time vs OMP_NUM_THREADS (nd=36) — thread scaling/saturation;
#             (b) per-RHS FMM time vs nd (at 16 threads) — amortization over RHS.
# julia --project plot_fmm_bench.jl  ->  figs/ + data//fmm_bench.png

using CairoMakie, DelimitedFiles

data, _ = readdlm(joinpath(@__DIR__, "data", "fmm_bench.csv"), ','; header = true)
omp   = Int.(data[:, 1]); level = Int.(data[:, 2]); N = Int.(data[:, 3])
nd    = Int.(data[:, 4]); tfmm  = Float64.(data[:, 5])
Ns = sort(unique(N))
cols = [:dodgerblue, :seagreen, :darkorange, :crimson]

sel(pred) = [(N[i], nd[i], omp[i], tfmm[i]) for i in eachindex(N) if pred(i)]

fig = Figure(size = (1120, 480))

# (a) per-RHS FMM vs threads, nd = 36
axa = Axis(fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = "OMP_NUM_THREADS", ylabel = "per-RHS FMM time (s)",
    title = "(a) thread scaling  (nd = 36, fmm_tol = 1e-4)")
for (c, Nv) in enumerate(Ns)
    rows = sort(sel(i -> N[i] == Nv && nd[i] == 36), by = x -> x[3])
    ts = [r[3] for r in rows]; pr = [r[4] / 36 for r in rows]
    scatterlines!(axa, ts, pr; color = cols[c], markersize = 9, label = "N=$Nv")
end
# ideal linear-scaling reference from the 1-thread, largest-N point
let rows = sort(sel(i -> N[i] == Ns[end] && nd[i] == 36), by = x -> x[3])
    t1 = rows[1][3]; p1 = rows[1][4] / 36
    tt = [r[3] for r in rows]
    lines!(axa, tt, p1 .* t1 ./ tt; color = :gray, linestyle = :dot, label = "ideal ∝ 1/threads")
end
axislegend(axa; position = :lb)

# (b) per-RHS FMM vs nd, at 16 threads
axb = Axis(fig[1, 2]; xscale = log10, yscale = log10,
    xlabel = "number of RHS  (nd)", ylabel = "per-RHS FMM time (s)",
    xticks = ([1, 4, 16, 36], ["1", "4", "16", "36"]),
    title = "(b) amortization over RHS  (16 threads)")
for (c, Nv) in enumerate(Ns)
    rows = sort(sel(i -> N[i] == Nv && omp[i] == 16), by = x -> x[2])
    nds = [r[2] for r in rows]; pr = [r[4] / r[2] for r in rows]
    scatterlines!(axb, nds, pr; color = cols[c], markersize = 9, label = "N=$Nv")
end
axislegend(axb; position = :lb)

out = joinpath(@__DIR__, "figs", "fmm_bench.png")
save(out, fig)
println("saved ", out)
