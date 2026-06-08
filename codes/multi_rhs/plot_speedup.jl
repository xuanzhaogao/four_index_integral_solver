# Per-RHS FMM speedup as a 2D function of (nd, OMP threads), relative to the nd=1 / 1-thread
# baseline, at N=461,376:  speedup(nd, thr) = nd * t_fmm(nd=1, 1 thr) / t_fmm(nd, thr).
# Along the threads axis this is pure thread scaling; along the nd axis it is the batched-RHS
# amortization (~4x ceiling); the top-right corner folds both. slab vs 3D-uniform side by side.
#
# julia --project=scripts plot_speedup.jl  ->  figs//fmm_speedup.png

using CairoMakie, DelimitedFiles

figdir  = joinpath(@__DIR__, "figs")
datadir = joinpath(@__DIR__, "data")
load(path) = begin
    d, _ = readdlm(path, ','; header = true)
    Dict((Int(d[i,1]), Int(d[i,3]), Int(d[i,4])) => Float64(d[i,5]) for i in axes(d,1))
end
slab = load(joinpath(datadir, "fmm_bench.csv"))
uni  = load(joinpath(datadir, "fmm_bench_uniform.csv"))

omps = sort(unique(k[1] for k in keys(slab)))
nds  = sort(unique(k[3] for k in keys(slab)))
NBIG = maximum(k[2] for k in keys(slab))

function speedmat(D)
    base = D[(1, NBIG, 1)]                       # nd=1, 1 thread
    M = fill(NaN, length(nds), length(omps))
    for i in eachindex(nds), j in eachindex(omps)
        k = (omps[j], NBIG, nds[i])
        haskey(D, k) && (M[i,j] = nds[i] * base / D[k])
    end
    M
end
Ms, Mu = speedmat(slab), speedmat(uni)
hi = maximum(filter(isfinite, vcat(vec(Ms), vec(Mu))))

fig = Figure(size = (1050, 460))
Label(fig[0, 1:4], "Per-RHS FMM speedup vs (nd=1, 1 thread)  —  N=$NBIG";
      fontsize = 16, font = :bold)

function heat!(col, M, title)
    ax = Axis(fig[1, col]; title = title, xlabel = "OMP_NUM_THREADS", ylabel = "nd (RHS per block)",
        xticks = (1:length(omps), string.(omps)), yticks = (1:length(nds), string.(nds)))
    hm = heatmap!(ax, 1:length(omps), 1:length(nds), permutedims(M);
        colormap = :viridis, colorrange = (1, hi))
    for i in eachindex(nds), j in eachindex(omps)
        isnan(M[i,j]) && continue
        text!(ax, j, i; text = string(round(M[i,j], digits = 1)) * "×",
            align = (:center, :center), fontsize = 12,
            color = M[i,j] > 0.55hi ? :black : :white)
    end
    ax, hm
end

_, h1 = heat!(1, Ms, "slab interface")
_, h2 = heat!(2, Mu, "3D uniform")
Colorbar(fig[1, 3], h2; label = "per-RHS speedup ×")

out = joinpath(figdir, "fmm_speedup.png")
save(out, fig)
println("saved ", out)
