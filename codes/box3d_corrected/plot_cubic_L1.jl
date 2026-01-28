using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/near_source_cubic_L1.csv"), DataFrame)

Ls = sort(unique(df.L))
ps = sort(filter(p -> p in (2, 3, 4, 5, 6), unique(df.p)))
rs = sort(unique(df.r))

df_ref = df[(df.p .== 7), :]
ref_pot = df_ref.potential_u[1]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    
    markers = [:xcross, :circle, :rect, :rtriangle, :diamond]
    colors = [:black, :red, :blue, :green, :purple]
    ms = 10

    ax = Axis(fig[1, 1], xlabel = "l_min", ylabel = "relative error of potential", xscale = log10, xreversed = true, yscale = log10, title = "adaptive mesh")

    ax2 = Axis(fig[1, 2], xlabel = "l_min", ylabel = "relative error of potential", xscale = log10, xreversed = true, yscale = log10, title = "non-adaptive mesh")

    for (j, p) in enumerate(ps)
        df_lp = df[(df.L .== 1.0) .& (df.p .== p), :]
        scatterlines!(ax, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_u .- ref_pot) ./ ref_pot), label = "p = $p", markersize = ms, marker = markers[j], color = colors[j])
        scatterlines!(ax2, 1.0 ./ (2 .^ df_lp.r), abs.((df_lp.potential_n .- ref_pot) ./ ref_pot), label = "p = $p", markersize = ms, marker = markers[j], color = colors[j])
    end

    axislegend(ax, position = :rt, nbanks = 2)
    ylims!(ax, 1e-4, 10^(-0.5))
    ylims!(ax2, 1e-4, 10^(1.0))

    fig
end

save(joinpath(@__DIR__, "figs/cubic_L1_r.svg"), fig)