using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/total_flux.csv"), DataFrame)

Ls = sort(unique(df.L))
ps = sort(unique(df.p))
rs = sort(unique(df.r))
ids = sort(unique(df.id))

begin
    fig = Figure(size = (300 * length(ps), 100 + 250 * length(Ls)), fontsize = 16)

    axs = []

    markers_id1 = [:xcross, :circle, :rect]
    markers_id2 = [:cross, :diamond, :rtriangle]
    colors = [:green, :red, :blue]
    ms = 10

    for (i, L) in enumerate(Ls)
        for (j, p) in enumerate(ps)
            ax = Axis(
                fig[i, j],
                xlabel = "desired accuracy",
                ylabel = "total flux error",
                xscale = log10,
                xreversed = true,
                yscale = log10,
                title = "L = $L, p = $p",
            )

            for (k, r) in enumerate(rs)
                # id = 1
                df_lpr1 = df[(df.L .== L) .& (df.p .== p) .& (df.r .== r) .& (df.id .== 1), :]
                scatter!(ax, df_lpr1.tol, abs.(df_lpr1.total_flux), label = "center, r=$r", markersize = ms, marker = markers_id1[k], color = colors[k])
                
                # id = 2
                df_lpr2 = df[(df.L .== L) .& (df.p .== p) .& (df.r .== r) .& (df.id .== 2), :]
                scatter!(ax, df_lpr2.tol, abs.(df_lpr2.total_flux), label = "corner, r=$r", markersize = ms, marker = markers_id2[k], color = colors[k])
            end

            # x = y line
            lines!(ax, [1e-10, 1e0], [1e-10, 1e0], color = :gray, linestyle = :dash)

            push!(axs, ax)
        end
    end

    Legend(fig[0, :], axs[1], orientation = :horizontal, nbanks = 1)

    fig
end

save(joinpath(@__DIR__, "figs/total_flux.svg"), fig)
