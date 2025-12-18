using CSV, DataFrames, JLD2
using CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/single_thin_box_convergence.csv"), DataFrame)

# L == 5.0
for L in [5.0, 10.0, 20.0]
    df_L = filter(row -> row.L == L, df)
    begin
        fig = Figure(size = (1000, 1200), fontsize = 20)

        for (i,nxy) in enumerate([5, 10, 20] .* Int(L / 5.0))
            nz = nxy ÷ Int(L)
            ax1 = Axis(fig[i, 1], xlabel = L"r", ylabel = "error", yscale = log10, title = "nxy = $(nxy), nz = $(nz)")
            ax2 = Axis(fig[i, 2], xlabel = L"p", ylabel = "error", yscale = log10, title = "nxy = $(nxy), nz = $(nz)")
            df_nxy = filter(row -> row.nxy == nxy, df_L)
            for p in unique(df_nxy.p)
                df_p = filter(row -> row.p == p, df_nxy)
                scatterlines!(ax1, df_p.r, abs.(1 .- df_p.gi), markersize = 12, label = "$(p)")
            end
            axislegend(ax1, position = :lb, "p", nbanks = 2)

            for r in unique(df_nxy.r)
                df_r = filter(row -> row.r == r, df_nxy)
                scatterlines!(ax2, df_r.p, abs.(1 .- df_r.gi), markersize = 12, label = "$(r)")
            end
            axislegend(ax2, position = :lb, "r", nbanks = 2)

            ylims!(ax1, 1e-4, 1e1)
            ylims!(ax2, 1e-4, 1e1)
        end

        save(joinpath(@__DIR__, "figs/single_thin_box_convergence_L$(Int(L)).svg"), fig)
        fig
    end
end

begin
    fig = Figure(size = (1000, 900), fontsize = 20)

    for (i, L) in enumerate([5.0, 10.0, 20.0])
        df_L = filter(row -> row.L == L, df)
        for (j, nxy) in enumerate([5, 10, 20] .* Int(L / 5.0))
            nz = nxy ÷ Int(L)
            ax1 = Axis(fig[i, j], xlabel = L"r", ylabel = "error", yscale = log10, title = "L = $(L), nxy = $(nxy), nz = $(nxy ÷ Int(L))")
            df_nxy = filter(row -> row.nxy == nxy, df_L)
            for p in unique(df_nxy.p)
                df_p = filter(row -> row.p == p, df_nxy)
                scatterlines!(ax1, df_p.r, abs.(1 .- df_p.gi), markersize = 12, label = "$(p)")
            end
            # axislegend(ax1, position = :lb, "p", nbanks = 2)
            ylims!(ax1, 1e-5, 1e0)
        end
    end

    Legend(fig[0, :], ax1, L"p", orientation = :horizontal)

    save(joinpath(@__DIR__, "figs/single_thin_box_convergence.svg"), fig)
    fig
end