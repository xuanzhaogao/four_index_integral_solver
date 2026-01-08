using CSV, DataFrames
using CairoMakie
import BoundaryIntegral as BI

df = CSV.read("data/doublelayer_upsampling_error.csv", DataFrame)

pt_zs = unique(df.pt_z)
trg_zs = unique(df.trg_z)
ns = unique(df.n)

colors = [:red, :blue, :green, :purple, :orange]
markers = [:circle, :rect, :utriangle, :diamond, :star4]

begin
    fig = Figure(size = (1000, 800), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "n_up", ylabel = "error", yscale = log10, title = "trg_z = 0.1")
    ax_2 = Axis(fig[1, 2], xlabel = "n_up", ylabel = "error", yscale = log10, title = "trg_z = 0.2")
    ax_3 = Axis(fig[2, 1], xlabel = "n_up", ylabel = "error", yscale = log10, title = "trg_z = 0.5")
    ax_4 = Axis(fig[2, 2], xlabel = "n_up", ylabel = "error", yscale = log10, title = "trg_z = 1.0")

    axs = [ax_1, ax_2, ax_3, ax_4]

    for (i, trg_z) in enumerate(trg_zs)
        df_trg_z = df[df.trg_z .== trg_z, :]
        rho = BI.bernstein_rho_2d(0.0, 0.0, trg_z)
        for (j, n) in enumerate(ns)
            df_n = df_trg_z[df_trg_z.n .== n, :]
            scatter!(axs[i], df_n.n_up, df_n.err, label = "n = $n", color = colors[j], marker = markers[j])
        end
        lines!(axs[i], df_trg_z.n_up, rho .^ (-2 .* df_trg_z.n_up), color = :black, linestyle = :dot)
    end

    axislegend(axs[1], position = :lb)

    for ax in axs
        ylims!(ax, 1e-8, 1e3)
    end

    fig
end

save("figs/doublelayer_upsampling_error.svg", fig)