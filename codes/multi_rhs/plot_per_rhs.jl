# Plot per-RHS solve time vs number of right-hand sides, from a PINNED-thread bench_per_rhs.jl run
# (data: figs//per_rhs.csv, columns omp,K,t_block,iters; the K=0 row carries the single-RHS
# baseline as t_block = per-RHS baseline). The earlier version hardcoded numbers from an UNPINNED
# run whose small-K/large-K points ran at different thread counts, producing a spurious ~18x that
# "tracked 1/K". With threads pinned the per-RHS solve time *plateaus* at the FMM amortization
# ceiling (~4x), it does NOT keep falling like 1/K.
#
# julia --project=scripts plot_per_rhs.jl  ->  figs//per_rhs_runtime.png

using CairoMakie, DelimitedFiles

figdir  = joinpath(@__DIR__, "figs")
datadir = joinpath(@__DIR__, "data")
data, _ = readdlm(joinpath(datadir, "per_rhs.csv"), ','; header = true)
omp  = Int.(data[:, 1]); K = Int.(data[:, 2]); tblk = Float64.(data[:, 3])
ompv = maximum(omp)

base_i = findfirst(==(0), K)                          # K=0 row = single-RHS baseline (s/RHS)
single = base_i === nothing ? NaN : tblk[base_i]
ord = sortperm(K); keep = [i for i in ord if K[i] >= 1]
Ks = K[keep]; tb = tblk[keep]; per = tb ./ Ks

amort = per[1] / per[end]                             # K=1 -> Kmax per-RHS amortization

fig = Figure(size = (820, 560))
ax = Axis(fig[1, 1];
    xlabel = "number of right-hand sides  K",
    ylabel = "per-RHS solve time  (s)",
    xscale = log10, yscale = log10,
    xticks = (Ks, string.(Ks)),
    title = "Multi-RHS amortization — graphene center 1, pinned $(ompv) threads")

# misleading 'ideal 1/K' line, shown dotted ONLY to contrast with reality
Kd = range(1, maximum(Ks); length = 200)
lines!(ax, collect(Kd), per[1] ./ Kd;
    color = :gray, linestyle = :dot, label = "ideal 1/K (NOT achieved)")
# the real floor: FMM amortization ceiling ~4x below the K=1 point
hlines!(ax, [per[1] / 4]; color = :seagreen, linestyle = :dashdot,
    label = "FMM amortization floor (≈4× below K=1)")
isfinite(single) && hlines!(ax, [single];
    color = :firebrick, linestyle = :dash, label = "single-RHS baseline ($(round(single,digits=1)) s/RHS)")
scatterlines!(ax, Ks, per;
    color = :dodgerblue, markersize = 12, linewidth = 2,
    label = "block GMRES  t_block(K)/K")

text!(ax, Ks[end], per[end] * 1.18; text = "$(round(amort,digits=1))× vs K=1\n(plateau, not 1/K)",
    align = (:right, :bottom), fontsize = 13, color = :dodgerblue)

axislegend(ax; position = :lb)
out = joinpath(figdir, "per_rhs_runtime.png")
save(out, fig)
println("saved ", out, "   (K=1→$(Ks[end]) per-RHS amortization = $(round(amort,digits=2))×)")
