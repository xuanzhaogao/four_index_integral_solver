using CSV, DataFrames
using CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/hubbard_graphene.csv"), DataFrame)

# reference values from 10.1103/PhysRevLett.106.236805 
refs = [17.0, 8.5, 5.4, 4.7]
labels = ["U_00", "U_01", "U_02", "U_03"]
colors = [:blue, :orange, :green, :red]
markers = [:circle, :xcross, :diamond, :utriangle]

begin
    fig = Figure(size = (1000, 450), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "N_FFT", ylabel = "U(eV)")
    ax_2 = Axis(fig[1, 2], xlabel = "N_FFT", ylabel = "self-convergence rel_err", yscale = log10)

    for (i, label) in enumerate(labels)
        df_subset = df[df.pair .== label, :]
        scatter!(ax_1, df_subset.N_FFT, df_subset.U_ev, label = label, markersize = 12, color = colors[i], marker = markers[i])
        ref_value = refs[i]
        hlines!(ax_1, [ref_value], color = colors[i], linestyle = :dash)

        @show abs.(df_subset.U_ev .- df_subset.U_ev[end]) ./ df_subset.U_ev[end]

        scatter!(ax_2, df_subset.N_FFT, abs.(df_subset.U_ev .- df_subset.U_ev[end]) ./ df_subset.U_ev[end], label = label, markersize = 12, color = colors[i], marker = markers[i])

        xlims!(ax_2, 0, 250)
        ylims!(ax_2, 1e-6, 1)
    end
    
    Legend(fig[0, :], ax_1, orientation = :horizontal)

    fig
end

save(joinpath(@__DIR__, "figs/hubbard_graphene_convergence.svg"), fig)