using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/near_source_thin_L1.csv"), DataFrame)
df_corrected = CSV.read(joinpath(@__DIR__, "data/near_source_thin_corrected.csv"), DataFrame)

Ls = sort(unique(df.L))
ps = sort(filter(p -> p in (2, 3, 4, 5, 6, 7), unique(df.p)))
rs = sort(filter(r -> r in (0, 2, 4, 6, 8), unique(df.r)))

df_ref = df[(df.p .== 7) .& (df.r .== 9), :]
ref_pot = df_ref.potential_u[1]

begin
    fig = Figure(size = (1500, 500), fontsize = 20)
    
    markers = [:xcross, :circle, :rect, :rtriangle, :diamond]
    colors = [:black, :red, :blue, :green, :purple]
    ms = 10

    ax = Axis(fig[1, 2], xlabel = "l_min", ylabel = "relative error of potential", xscale = log10, xreversed = true, yscale = log10, title = "adaptive mesh")

    ax_corrected = Axis(fig[1, 3], xlabel = "l_min", ylabel = "relative error of potential", xscale = log10, xreversed = true, yscale = log10, title = "adaptive mesh, corrected")

    ax2 = Axis(fig[1, 1], xlabel = "l_min", ylabel = "relative error of potential", xscale = log10, xreversed = true, yscale = log10, title = "non-adaptive mesh")

    for (j, p) in enumerate(ps)
        df_lp = df[(df.L .== 10.0) .& (df.p .== p) .& (df.r .!= 9), :]
        df_lp_corrected = df_corrected[(df_corrected.L .== 10.0) .& (df_corrected.p .== p), :]
        # df_lp_corrected = df_corrected[(df_corrected.L .== 1.0) .& (df_corrected.p .== p), :]
        scatterlines!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_u .- ref_pot) ./ ref_pot), label = "p = $p", markersize = ms, marker = markers[j], color = colors[j])

        scatterlines!(ax_corrected, 1.0 ./ (2 .^ df_lp_corrected.r), abs.((df_lp_corrected.potential .- ref_pot) ./ ref_pot), label = "p = $p", markersize = ms, marker = markers[j], color = colors[j])

        scatterlines!(ax2, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_n .- ref_pot) ./ ref_pot), label = "p = $p", markersize = ms, marker = markers[j], color = colors[j])
    end

    # axislegend(ax, position = :lb, nbanks = 3)
    Legend(fig[0, :], ax, orientation = :horizontal)
    ylims!(ax, 1e-6, 10^(-1))
    # ylims!(ax2, 1e-6, 10^(-1))
    ylims!(ax_corrected, 1e-6, 10^(-1))

    fig
end

save(joinpath(@__DIR__, "figs/thin_near_source_r.svg"), fig)