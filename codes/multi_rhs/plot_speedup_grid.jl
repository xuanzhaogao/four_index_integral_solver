# Slab FMM speedup at fixed N=461,376, as a 3-panel figure:
#   (a) per-RHS nd-amortization speedup vs nd, one line per thread count
#       S_a(nd,T) = nd * t(1,T) / t(nd,T)          [vs nd=1 at the same thread count]
#   (b) thread speedup vs OMP_NUM_THREADS, one line per nd
#       S_b(T,nd) = t(1,nd) / t(T,nd)              [vs 1 thread at the same nd]
#   (c) combined per-RHS speedup heatmap over nd x threads
#       S_c(nd,T) = nd * t(1,1) / t(nd,T)          [vs (nd=1, 1 thread)]
# Data: data/fmm_bench.csv (slab kernel).
#
# julia --project=. scripts? -> here: julia --project plot_speedup_grid.jl -> figs/fmm_speedup_grid.png

using CairoMakie, DelimitedFiles

datadir = joinpath(@__DIR__, "data")
figdir  = joinpath(@__DIR__, "figs")
d, _ = readdlm(joinpath(datadir, "fmm_bench.csv"), ','; header = true)
omp = Int.(d[:, 1]); N = Int.(d[:, 3]); nd = Int.(d[:, 4]); t = Float64.(d[:, 5])

const N4 = 461376
T = Dict((omp[i], nd[i]) => t[i] for i in eachindex(t) if N[i] == N4)
omps = sort(unique(o for (o, _) in keys(T)))
nds  = sort(unique(n for (_, n) in keys(T)))

begin

    fig = Figure(size = (1480, 450), fontsize = 18)

    # (a) nd amortization, per thread count
    axa = Axis(fig[1, 1]; xscale = log10, xlabel = "number of RHS  nd",
        ylabel = "per-RHS speedup vs nd=1", xticks = (nds, string.(nds)),
        title = "(a) nd speedup")
    thr_a = [1, 4, 16, 64, 96]
    ca = cgrad(:viridis, length(thr_a); categorical = true)
    for (k, Tt) in enumerate(thr_a)
        ys = [nds[j] * T[(Tt, 1)] / T[(Tt, nds[j])] for j in eachindex(nds)]
        scatterlines!(axa, nds, ys; color = ca[k], markersize = 9, label = "$Tt thr")
    end
    axislegend(axa; position = :lt, framevisible = false)
    ylims!(axa, 0.9, 4.1)

    # (b) thread scaling, per nd  (log-log so ideal y=x is a diagonal)
    axb = Axis(fig[1, 2]; xscale = log10, yscale = log10,
        xlabel = "OMP_NUM_THREADS", ylabel = "speedup vs 1 thread",
        xticks = (omps, string.(omps)), yticks = ([1, 2, 5, 10, 20, 50], string.([1, 2, 5, 10, 20, 50])),
        title = "(b) thread speedup")
    cb = cgrad(:plasma, length(nds); categorical = true)
    for (k, n) in enumerate(nds)
        ys = [T[(1, n)] / T[(o, n)] for o in omps]
        scatterlines!(axb, omps, ys; color = cb[k], markersize = 9, label = "nd=$n")
    end
    lines!(axb, omps, Float64.(omps); color = :gray, linestyle = :dot, label = "ideal ∝ threads")
    axislegend(axb; position = :lt, framevisible = false)

    # (c) combined per-RHS speedup heatmap, nd x threads
    axc = Axis(fig[1, 3]; xlabel = "OMP_NUM_THREADS", ylabel = "nd",
        xticks = (1:length(omps), string.(omps)), yticks = (1:length(nds), string.(nds)),
        title = "(c) combined speedup vs (nd=1, 1 thr)")
    M = [nds[i] * T[(1, 1)] / T[(omps[j], nds[i])] for i in 1:length(nds), j in 1:length(omps)]
    hm = heatmap!(axc, 1:length(omps), 1:length(nds), permutedims(M); colormap = :viridis)
    hi = maximum(M)
    for i in 1:length(nds), j in 1:length(omps)
        text!(axc, j, i; text = string(round(M[i, j], digits = 1)) * "x",
            align = (:center, :center), fontsize = 10,
            color = M[i, j] > 0.55hi ? :black : :white)
    end
    Colorbar(fig[1, 4], hm; label = "per-RHS speedup ×")

    fig
end

out = joinpath(figdir, "fmm_speedup_grid.pdf")
save(out, fig)
println("saved ", out, "   peak combined = ", round(hi, digits = 1), "x")
