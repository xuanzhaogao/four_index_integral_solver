# FMM3D throughput (M source·RHS evaluations / s = N*nd / t_fmm) as a function of each of the
# three knobs — number of unknowns N, OMP threads, and number of RHS nd — comparing the slab BIE
# interface (figs//fmm_bench.csv) against a 3D-uniform cloud (figs//fmm_bench_uniform.csv).
# Each panel sweeps ONE knob with the other two held at a representative value (threads=64, nd=36,
# N=461,376). slab = solid, uniform = dashed.
#
# julia --project=scripts plot_throughput.jl  ->  figs//fmm_throughput.png

using CairoMakie, DelimitedFiles

figdir  = joinpath(@__DIR__, "figs")
datadir = joinpath(@__DIR__, "data")
load(path) = begin
    d, _ = readdlm(path, ','; header = true)
    Dict((Int(d[i,1]), Int(d[i,3]), Int(d[i,4])) => Float64(d[i,5]) for i in axes(d,1))
end
slab = load(joinpath(datadir, "fmm_bench.csv"))
uni  = load(joinpath(datadir, "fmm_bench_uniform.csv"))
tput(D, o, N, nd) = haskey(D,(o,N,nd)) ? N*nd/D[(o,N,nd)]/1e6 : NaN

Ns   = sort(unique(k[2] for k in keys(slab)))
omps = sort(unique(k[1] for k in keys(slab)))
nds  = sort(unique(k[3] for k in keys(slab)))
THR, ND, NBIG = 64, 36, maximum(Ns)

fig = Figure(size = (1320, 460))
Label(fig[0, 1:3], "FMM3D throughput  (M source·RHS evals / s)  —  slab (solid) vs 3D-uniform (dashed)";
      fontsize = 16, font = :bold)

function panel(col, xs, ys_slab, ys_uni, xlabel, title; xlog = true)
    ax = Axis(fig[1, col]; xlabel = xlabel, ylabel = "throughput (M evals/s)",
        xscale = xlog ? log10 : identity, title = title,
        xticks = (xs, string.(xs)))
    scatterlines!(ax, xs, ys_slab; color = :crimson, markersize = 11, linewidth = 2.4, label = "slab")
    scatterlines!(ax, xs, ys_uni;  color = :royalblue, markersize = 11, linewidth = 2.4,
        linestyle = :dash, marker = :rect, label = "uniform")
    axislegend(ax; position = :lt)
    ax
end

# Panel 1: vs N   (threads=THR, nd=ND)
panel(1, Ns,  [tput(slab,THR,N,ND) for N in Ns], [tput(uni,THR,N,ND) for N in Ns],
      "N (interface unknowns)", "vs N   (threads=$THR, nd=$ND)")
# Panel 2: vs threads   (N=NBIG, nd=ND)
panel(2, omps, [tput(slab,o,NBIG,ND) for o in omps], [tput(uni,o,NBIG,ND) for o in omps],
      "OMP_NUM_THREADS", "vs threads   (N=$NBIG, nd=$ND)")
# Panel 3: vs nd   (N=NBIG, threads=THR)
panel(3, nds, [tput(slab,THR,NBIG,nd) for nd in nds], [tput(uni,THR,NBIG,nd) for nd in nds],
      "number of RHS  (nd)", "vs nd   (N=$NBIG, threads=$THR)")

out = joinpath(figdir, "fmm_throughput.png")
save(out, fig)
println("saved ", out)
