using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/single_box3d_nonadapt.csv"), DataFrame)

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"p", ylabel = L"\mathcal{E}_r", yscale = log10, title = "Relative Error in Potential", xticks = [0:8:72...])
    ax2 = Axis(fig[1, 2], xlabel = L"p", ylabel = L"\mathcal{E}_r", yscale = log10, title = "Error in Flux", xticks = [0:8:72...])

    scatterlines!(ax, df.n_quad, df.pot_relerr, marker = :circle, markersize = 12)
    xlims!(ax, 0, 73)
    ylims!(ax, 1e-3, 1e1)

    scatterlines!(ax2, df.n_quad, abs.(1 .- df.flux), marker = :circle, markersize = 12)
    xlims!(ax2, 0, 73)
    ylims!(ax2, 1e-3, 1e1)

    save(joinpath(@__DIR__, "figs/single_box3d_nonadapt.svg"), fig)
    fig
end