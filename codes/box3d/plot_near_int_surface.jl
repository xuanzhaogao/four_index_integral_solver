using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/near_int_surface.csv"), DataFrame)

ds = unique(df.δ)
ps = unique(df.p)

begin
    fig = Figure(size = (600, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"p", ylabel = L"\mathcal{E}_r", title = "Near Integral", yscale = log10)
    for d in ds
        df_d = filter(row -> row.δ == d, df)
        errors = abs.((df_d.t .- df_d.t[end]) ./ df_d.t[end])
        lines!(ax, df_d.p, errors, label = "$(d)")
    end
    # axislegend(ax, position = :lb)
    Legend(fig[1, 2], ax, L"\delta")
    xlims!(ax, 0, 256)
    ylims!(ax, 1e-16, 1e1)
    save(joinpath(@__DIR__, "figs/near_int_surface.svg"), fig)
    fig
end