# Performance heatmaps of the FMM3D kernel: throughput (million source·RHS evaluations / s,
# = N*nd / t_fmm) over OMP_NUM_THREADS x N, at nd = 36. Compares the slab BIE interface
# (figs//fmm_bench.csv) against a 3D-uniform point cloud at the same N
# (figs//fmm_bench_uniform.csv, if present).
#
# julia --project plot_fmm_heatmap.jl  ->  figs//fmm_heatmap.png

using CairoMakie, DelimitedFiles

read5(path) = begin
    d, _ = readdlm(path, ','; header = true)
    (omp = Int.(d[:, 1]), N = Int.(d[:, 3]), nd = Int.(d[:, 4]), t = Float64.(d[:, 5]))
end

figdir  = joinpath(@__DIR__, "figs")
datadir = joinpath(@__DIR__, "data")
slab = read5(joinpath(datadir, "fmm_bench.csv"))
uni_path = joinpath(datadir, "fmm_bench_uniform.csv")
have_uni = isfile(uni_path)
uni = have_uni ? read5(uni_path) : nothing

ND = 36
threads = sort(unique(slab.omp))
Ns = sort(unique(slab.N))

# throughput matrix (Mevals/s): rows = N, cols = threads, at nd=ND
function tput(b)
    M = fill(NaN, length(Ns), length(threads))
    for i in eachindex(Ns), j in eachindex(threads)
        k = findfirst(q -> b.N[q] == Ns[i] && b.omp[q] == threads[j] && b.nd[q] == ND, eachindex(b.N))
        k === nothing && continue
        M[i, j] = Ns[i] * ND / b.t[k] / 1e6
    end
    M
end

Mslab = tput(slab)
panels = have_uni ? 3 : 1
fig = Figure(size = (430 * panels, 430))

function heat!(col, M, title)
    ax = Axis(fig[1, col]; title = title,
        xlabel = "OMP_NUM_THREADS", ylabel = "N (interface unknowns)",
        xticks = (1:length(threads), string.(threads)),
        yticks = (1:length(Ns), string.(Ns)))
    hm = heatmap!(ax, 1:length(threads), 1:length(Ns), permutedims(M);
        colormap = :viridis)
    for i in eachindex(Ns), j in eachindex(threads)
        isnan(M[i, j]) && continue
        text!(ax, j, i; text = string(round(M[i, j], digits = 2)),
            align = (:center, :center), color = :white, fontsize = 11)
    end
    Colorbar(fig[1, col + 1], hm; label = "M evals/s")
    ax
end

if have_uni
    Muni = tput(uni)
    heat!(1, Mslab, "slab interface  (throughput, nd=36)")
    heat!(3, Muni, "3D uniform  (throughput, nd=36)")
    # ratio uniform/slab
    axr = Axis(fig[1, 5]; title = "uniform / slab  (speed ratio)",
        xlabel = "OMP_NUM_THREADS", ylabel = "N",
        xticks = (1:length(threads), string.(threads)),
        yticks = (1:length(Ns), string.(Ns)))
    R = Muni ./ Mslab
    hmr = heatmap!(axr, 1:length(threads), 1:length(Ns), permutedims(R); colormap = :magma)
    for i in eachindex(Ns), j in eachindex(threads)
        isnan(R[i, j]) && continue
        text!(axr, j, i; text = string(round(R[i, j], digits = 1)) * "x",
            align = (:center, :center), color = :white, fontsize = 11)
    end
    Colorbar(fig[1, 6], hmr; label = "ratio")
else
    heat!(1, Mslab, "slab interface  (throughput, nd=36)")
end

out = joinpath(figdir, "fmm_heatmap.png")
save(out, fig)
println("saved ", out)
