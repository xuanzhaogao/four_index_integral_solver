using CSV, DataFrames
using CairoMakie

df = CSV.read("data/single_box2d_gi.csv", DataFrame)

n_quads = unique(df.n_quad)
n_adapts = unique(df.n_adapt)

begin
    fig = Figure(size = (1500, 500), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "N_adapt", ylabel = "L1 Error", yscale = log10)
    ax2 = Axis(fig[1, 2], xlabel = "N_adapt", ylabel = "L2 Error", yscale = log10)
    ax3 = Axis(fig[1, 3], xlabel = "N_adapt", ylabel = "Linf Error", yscale = log10)

    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df)
        scatterlines!(ax1, df_quad.n_adapt, df_quad.L1_err, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
    end

    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df)
        scatterlines!(ax2, df_quad.n_adapt, df_quad.L2_err, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
    end
    
    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df)
        scatterlines!(ax3, df_quad.n_adapt, df_quad.Linf_err, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
    end

    Legend(fig[0, :], ax1, nbanks = 1, orientation = :horizontal)

    save("figs/single_box2d_gi.svg", fig)

    fig
end