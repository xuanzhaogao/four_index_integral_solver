using CSV, DataFrames
using CairoMakie

df = CSV.read("data/benchmark_fmm2d.csv", DataFrame)

n_points = df.n_points
time_map = df.time_map
time_singleiter = df.time_singleiter
time_alliter = df.time_alliter

begin
    fig = Figure(size = (700, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Number of points", ylabel = "Time (s)", title = "FMM2D Run Time", xscale = log10, yscale = log10)
    scatter!(ax, n_points, time_map, color = :blue, marker = :diamond, markersize = 10, label = "Linear operator")
    scatter!(ax, n_points, time_singleiter, color = :red, marker = :circle, markersize = 10, label = "Single GMRES iteration")
    scatter!(ax, n_points, time_alliter, color = :green, marker = :utriangle, markersize = 10, label = "All GMRES iterations")
    Legend(fig[1, 2], ax, nbanks = 1, labelsize = 15)
    save("figs/benchmark_fmm2d_plot.svg", fig)
    fig
end

df_fortran_1 = CSV.read("data/rfmm2d_fortran_benchmark_1.csv", DataFrame)
df_fortran_2 = CSV.read("data/rfmm2d_fortran_benchmark_2.csv", DataFrame)
df_fortran_4 = CSV.read("data/rfmm2d_fortran_benchmark_4.csv", DataFrame)
df_fortran_8 = CSV.read("data/rfmm2d_fortran_benchmark_8.csv", DataFrame)
df_fortran_16 = CSV.read("data/rfmm2d_fortran_benchmark_16.csv", DataFrame)
df_julia_1 = CSV.read("data/rfmm2d_julia_benchmark_1.csv", DataFrame)
df_julia_2 = CSV.read("data/rfmm2d_julia_benchmark_2.csv", DataFrame)
df_julia_4 = CSV.read("data/rfmm2d_julia_benchmark_4.csv", DataFrame)
df_julia_8 = CSV.read("data/rfmm2d_julia_benchmark_8.csv", DataFrame)
df_julia_16 = CSV.read("data/rfmm2d_julia_benchmark_16.csv", DataFrame)

begin
    fig = Figure(size = (650, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Number of points", ylabel = "Time (s)", title = "FMM2D Run Time", xscale = log10, yscale = log10)

    scatterlines!(ax, df_fortran_1.ns, df_fortran_1.elapsed, label = "Fortran 1 thread", linewidth = 2, color = :blue)
    scatterlines!(ax, df_fortran_2.ns, df_fortran_2.elapsed, label = "Fortran 2 threads", linewidth = 2, color = :red)
    scatterlines!(ax, df_fortran_4.ns, df_fortran_4.elapsed, label = "Fortran 4 threads", linewidth = 2, color = :green)
    scatterlines!(ax, df_fortran_8.ns, df_fortran_8.elapsed, label = "Fortran 8 threads", linewidth = 2, color = :orange)
    scatterlines!(ax, df_fortran_16.ns, df_fortran_16.elapsed, label = "Fortran 16 threads", linewidth = 2, color = :purple)
    scatterlines!(ax, df_julia_1.ns, df_julia_1.time, label = "Julia 1 thread", linewidth = 2, marker = :diamond, color = :blue)
    scatterlines!(ax, df_julia_2.ns, df_julia_2.time, label = "Julia 2 threads", linewidth = 2, marker = :diamond, color = :red)
    scatterlines!(ax, df_julia_4.ns, df_julia_4.time, label = "Julia 4 threads", linewidth = 2, marker = :diamond, color = :green)
    scatterlines!(ax, df_julia_8.ns, df_julia_8.time, label = "Julia 8 threads", linewidth = 2, marker = :diamond, color = :orange)
    scatterlines!(ax, df_julia_16.ns, df_julia_16.time, label = "Julia 16 threads", linewidth = 2, marker = :diamond, color = :purple)

    Legend(fig[1, 2], ax, position = :lt, labelsize = 12, nbanks = 1)
    save("figs/benchmark_fmm2d_comparison.svg", fig)
    fig
end