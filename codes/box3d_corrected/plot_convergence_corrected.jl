using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/convergence_corrected.csv"), DataFrame)

ps = unique(df.p)

begin
    fig = Figure(size = (1000, 450), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "r", ylabel = "Error Flux", yscale = log10, title = "Uncorrected")
    ax_2 = Axis(fig[1, 2], xlabel = "r", ylabel = "Error Flux", yscale = log10, title = "Corrected")

    for p in ps
        df_p = df[df.p .== p, :]
        scatter!(ax_1, df_p.r, abs.(df_p.total_flux_u .- 1.0), label = "p = $p", markersize = 10)
        scatter!(ax_2, df_p.r, abs.(df_p.total_flux_c .- 1.0), label = "p = $p", markersize = 10)
    end
    
    axislegend(ax_1, position = :lb)

   ylims!(ax_1, 1e-6, 1e-0)
   ylims!(ax_2, 1e-6, 1e-0)

    fig
end