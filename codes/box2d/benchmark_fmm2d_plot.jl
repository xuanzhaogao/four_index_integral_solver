using CSV, DataFrames
using CairoMakie

df = CSV.read("data/benchmark_fmm2d.csv", DataFrame)

n_points = df.n_points
time_map = df.time_map
time_singleiter = df.time_singleiter
time_alliter = df.time_alliter

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Number of points", ylabel = "Time (s)", title = "FMM2D Run Time", xscale = log10, yscale = log10)
    scatter!(ax, n_points, time_map, color = :blue, marker = :diamond, markersize = 10, label = "Linear operator")
    # scatter!(ax, n_points, time_singleiter, color = :red, marker = :diamond, markersize = 10, label = "Single GMRES iteration")
    # scatter!(ax, n_points, time_alliter, color = :green, marker = :triangle, markersize = 10, label = "All GMRES iterations")
    axislegend(ax, position = :lt)
    save("figs/benchmark_fmm2d_plot.svg", fig)
    fig
end

df_fortran = CSV.read("data/rfmm2d_fortran_benchmark.csv", DataFrame)
df_julia = CSV.read("data/rfmm2d_julia_benchmark.csv", DataFrame)

begin
    fig = Figure(size = (650, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Number of points", ylabel = "Time (s)", title = "FMM2D Run Time", xscale = log10, yscale = log10)

    for n_threads in unique(df_fortran.num_thread)
        ns = df_fortran[df_fortran.num_thread .== n_threads, :ns]
        time = df_fortran[df_fortran.num_thread .== n_threads, :elapsed]
        scatterlines!(ax, ns, time, label = "Fortran $(n_threads) threads", linewidth = 2)
    end

    ns = df_julia.ns
    time = df_julia.time
    scatterlines!(ax, ns, time, label = "Julia -t auto", color = :red, linewidth = 2)

    Legend(fig[1, 2], ax, position = :lt, labelsize = 12, nbanks = 1)
    save("figs/benchmark_fmm2d_comparison.svg", fig)
    fig
end