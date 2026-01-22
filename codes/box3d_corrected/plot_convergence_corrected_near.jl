using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/convergence_corrected_near.csv"), DataFrame)
df_ref = CSV.read(joinpath(@__DIR__, "data/convergence_corrected_near_ref.csv"), DataFrame)

Ls = sort(unique(df_ref.L))
ps = sort(filter(p -> p in (2, 3, 4, 5, 6), unique(df.p)))
rs = sort(unique(df.r))

begin
    fig = Figure(size = (900, 100 + 250 * length(Ls)), fontsize = 16)

    axs = []

    markers = [:xcross, :circle, :rect, :rtriangle, :diamond]
    colors = [:black, :red, :blue, :green, :purple]
    ms = 10

    for (i, L) in enumerate(Ls)

        ref_pot = df_ref[(df_ref.L .== L), :].potential_cff[1]
        # ref_pot = df[(df.L .== L) .& (df.p .== 6) .& (df.r .== 6), :].potential_ctt[1]

        for (j, p) in enumerate([2, 4, 6])
            ax = Axis(
                fig[i, j],
                xlabel = "l_min",
                ylabel = "relative error of potential",
                xscale = log2,
                xreversed = true,
                yscale = log10,
                title = "L = $L, p = $p",
            )

            df_lp = df[(df.L .== L) .& (df.p .== p), :]
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_u .- ref_pot) ./ ref_pot), label = "uncorrected", markersize = ms, marker = markers[1], color = colors[1])
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_cff .- ref_pot) ./ ref_pot), label = "corrected (ff)", markersize = ms, marker = markers[2], color = colors[2])
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_ctf .- ref_pot) ./ ref_pot), label = "corrected (tf)", markersize = ms, marker = markers[3], color = colors[3])
            scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_cft .- ref_pot) ./ ref_pot), label = "corrected (ft)", markersize = ms, marker = markers[4], color = colors[4])
            # scatter!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_ctt .- ref_pot) ./ ref_pot), label = "corrected (tt)", markersize = ms, marker = markers[5], color = colors[5])

            # if i == 1 && j == 1
                # axislegend(ax, position = :lb)
            # end

            xlims!(ax, 2^(0.5), 2^(-6.5))
            ylims!(ax, 1e-5, 1e1)

            push!(axs, ax)
        end
    end

    Legend(fig[0, :], axs[1], orientation = :horizontal)

    fig
end

save(joinpath(@__DIR__, "figs/convergence_corrected_near_r.svg"), fig)

# Second plot: fix r, use p as x-axis
begin
    fig2 = Figure(size = (900, 100 + 250 * length(Ls)), fontsize = 16)

    axs2 = []

    markers = [:xcross, :circle, :rect, :rtriangle, :diamond]
    colors = [:black, :red, :blue, :green, :purple]
    ms = 10

    for (i, L) in enumerate(Ls)

        ref_pot = df_ref[(df_ref.L .== L), :].potential_cff[1]

        for (j, r) in enumerate([2, 4, 6])
            ax = Axis(
                fig2[i, j],
                xlabel = "p",
                ylabel = "relative error of potential",
                yscale = log10,
                title = "L = $L, r = $r",
            )

            df_lr = df[(df.L .== L) .& (df.r .== r), :]
            scatter!(ax, df_lr.p, abs.((df_lr.potential_u .- ref_pot) ./ ref_pot), label = "uncorrected", markersize = ms, marker = markers[1], color = colors[1])
            scatter!(ax, df_lr.p, abs.((df_lr.potential_cff .- ref_pot) ./ ref_pot), label = "corrected (ff)", markersize = ms, marker = markers[2], color = colors[2])
            scatter!(ax, df_lr.p, abs.((df_lr.potential_ctf .- ref_pot) ./ ref_pot), label = "corrected (tf)", markersize = ms, marker = markers[3], color = colors[3])
            scatter!(ax, df_lr.p, abs.((df_lr.potential_cft .- ref_pot) ./ ref_pot), label = "corrected (ft)", markersize = ms, marker = markers[4], color = colors[4])
            # scatter!(ax, df_lr.p, abs.((df_lr.potential_ctt .- ref_pot) ./ ref_pot), label = "corrected (tt)", markersize = ms, marker = markers[5], color = colors[5])

            ylims!(ax, 1e-5, 1e1)

            push!(axs2, ax)
        end
    end

    Legend(fig2[0, :], axs2[1], orientation = :horizontal)

    fig2
end

save(joinpath(@__DIR__, "figs/convergence_corrected_near_p.svg"), fig2)
