using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/convergence_corrected_far.csv"), DataFrame)

Ls = sort(unique(df.L))
ps = sort(filter(p -> p in (2, 4, 6), unique(df.p)))

begin
    fig = Figure(size = (900, 900), fontsize = 16)

    for (i, L) in enumerate(Ls)
        for (j, p) in enumerate(ps)
            ax = Axis(
                fig[i, j],
                xlabel = "l_min",
                ylabel = "Error total flux",
                xscale = log10,
                xreversed = true,
                yscale = log10,
                title = "L = $L, p = $p",
            )

            df_lp = df[(df.L .== L) .& (df.p .== p), :]
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.(df_lp.total_flux_u), label = "uncorrected", markersize = 8)
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.(df_lp.total_flux_c), label = "corrected", markersize = 8)

            if i == 1 && j == 1
                axislegend(ax, position = :lb)
            end

            ylims!(ax, 1e-6, 1e0)
        end
    end

    fig
end

save(joinpath(@__DIR__, "figs/convergence_corrected_far.svg"), fig)