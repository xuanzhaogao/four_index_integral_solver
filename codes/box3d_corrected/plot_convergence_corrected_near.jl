using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/convergence_corrected_near.csv"), DataFrame)
df_varquad = CSV.read(joinpath(@__DIR__, "data/convergence_corrected_near_varquad.csv"), DataFrame)

Ls = sort(unique(df.L))
ps = sort(filter(p -> p in (4, 6), unique(df.p)))
rs = sort(unique(df.r))
dzs = sort(unique(df.dz))

df_ref = df[df.p .== 7, :]

markers = [:xcross, :circle, :rtriangle, :diamond, :utriangle, :pentagon, :hexagon, :heptagon]
colors = [:black, :red, :blue, :green, :purple, :orange, :cyan, :magenta, :brown]
ms = 10

begin
    fig = Figure(size = (300 * length(dzs) + 100, 400), fontsize = 16)

    axs = []

    dzss = [0.01, 0.1, 1.0]


    for (i, dz) in enumerate(dzs)

        ref_pot = df_ref[(df_ref.dz .== dz), :].potential_cff[1]
        # ref_pot = df[(df.L .== L) .& (df.p .== 6) .& (df.r .== 6), :].potential_ctt[1]

        ax = Axis(
            fig[1, i],
            xlabel = "l_min",
            ylabel = "relative error of potential",
            xscale = log10,
            xreversed = true,
            yscale = log10,
            title = "dz = $(dzss[i])",
        )

        for (j, p) in enumerate([4, 6])
            
            df_p = df[(df.dz .== dz) .& (df.p .== p), :]
            scatter!(ax, df_p.l_ec, abs.((df_p.potential_u .- ref_pot) ./ ref_pot), label = "uncorrected, p=$p", markersize = ms, marker = markers[3 * (j - 1) + 1], color = colors[2 * (j - 1) + 1])
            scatter!(ax, df_p.l_ec, abs.((df_p.potential_cff .- ref_pot) ./ ref_pot), label = "corrected, p=$p", markersize = ms, marker = markers[3 * (j - 1) + 2], color = colors[2 * (j - 1) + 2])
            # scatter!(ax, df_p.l_ec, abs.((df_p.potentail_cft .- ref_pot) ./ ref_pot), label = "corrected (ft), p=$p", markersize = ms, marker = markers[3 * (j - 1) + 3], color = colors[3])

            # xlims!(ax, 10^(0.1), 10^(-2.1))
            # ylims!(ax, 1e-5, 1e1)
            push!(axs, ax)
        end
    end

    Legend(fig[0, :], axs[1], orientation = :horizontal, nbanks = 2)

    fig
end

save(joinpath(@__DIR__, "figs/convergence_corrected_near_r.svg"), fig)

begin
    fig = Figure(size = (500, 400), fontsize = 16)
    ax = Axis(
        fig[1, 1],
        xlabel = "l_min",
        ylabel = "relative error of potential",
        xscale = log10,
        xreversed = true,
        yscale = log10,
        title = "Variable quadrature points",
    )

    ref_pot = df_ref[(df_ref.dz .≈ 0.1), :].potential_cff[1]

    for (j, p) in enumerate([4, 6, 8])
        
        df_p = df_varquad[(df_varquad.dz .≈ 0.01) .& (df_varquad.p .== p), :]
        scatterlines!(ax, df_p.l_ec, abs.((df_p.potential_cff .- ref_pot) ./ ref_pot), label = "p=$p", markersize = ms, marker = markers[j], color = colors[j])
    end

    axislegend(ax, position = :rb)

    # xlims!(ax, 10^(0.1), 10^(-2.1))
    # ylims!(ax, 1e-5, 1e1)

    fig
end

save(joinpath(@__DIR__, "figs/convergence_corrected_near_varquad_r.svg"), fig)