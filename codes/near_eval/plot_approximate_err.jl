using CSV, DataFrames
using CairoMakie
import BoundaryIntegral as BI

df = CSV.read("data/doublelayer_represent_error.csv", DataFrame)


colors = [:red, :blue, :green, :purple]
markers = [:circle, :rect, :utriangle, :diamond]

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "n", ylabel = "error", yscale = log10)

    pt_zs = unique(df.pt_z)
    for (i, pt_z) in enumerate(pt_zs)
        df_pt_z = df[df.pt_z .== pt_z, :]
        scatter!(ax, df_pt_z.n, df_pt_z.err, label = "z = $pt_z", markersize = 10, color = colors[i], marker = markers[i])

        rho = BI.bernstein_rho_2d(0.2, 0.3, pt_z)
        ratio = df_pt_z.err[1] / rho .^ (-2)
        lines!(ax, df_pt_z.n, ratio .* rho .^ (- df_pt_z.n), color = colors[i])
    end

    axislegend(ax, position = :lb)
    ylims!(ax, 1e-16, 1e4)

    fig
end

save("figs/doublelayer_represent_error.svg", fig)
