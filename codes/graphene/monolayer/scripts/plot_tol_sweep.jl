using CSV
using DataFrames
using CairoMakie

const CSV_PATH = joinpath(@__DIR__, "..", "data", "bare_monolayer_tol_sweep.csv")
const FIG_DIR  = joinpath(@__DIR__, "..", "..", "figs")

const CHANNELS = ("onsite", "nn", "hund_sf", "hund_ph")
const CHANNEL_LABEL = Dict(
    "onsite"  => "onsite (intra-orbital U)",
    "nn"      => "nearest-neighbor V",
    "hund_sf" => "Hund spin-flip J",
    "hund_ph" => "Hund pair-hopping J'",
)
const CHANNEL_COLOR = Dict(
    "onsite"  => :crimson,
    "nn"      => :royalblue,
    "hund_sf" => :seagreen,
    "hund_ph" => :darkorange,
)

const REFERENCE_EV = Dict(
    :onsite  => 17.434191,
    :nn      => 8.839615,
    :hund_sf => 0.130804,
    :hund_ph => 0.130804,
)

const PRL_REF_EV = Dict(
    :onsite  => 17.0,
    :nn      => 8.5,
)

# function main()
#     df = CSV.read(CSV_PATH, DataFrame)
#     sort!(df, [:channel, :tol], rev = [false, true])  # tol descending so we sweep tighter

#     # --- Figure 1: u_ev vs tol, one panel per channel, with CoQui ref overlaid ---
#     fig1 = Figure(; size = (1100, 800))
#     for (i, ch) in enumerate(CHANNELS)
#         sub = df[df.channel .== ch, :]
#         sort!(sub, :tol, rev = true)
#         row, col = fldmod1(i, 2)
#         ax = Axis(fig1[row, col];
#             title = CHANNEL_LABEL[ch],
#             xlabel = "source_tol = volume_tol",
#             ylabel = "U (eV)",
#             xscale = log10,
#             xreversed = true,
#         )
#         scatterlines!(ax, sub.tol, sub.u_ev;
#             color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
#             label = "this run")
#         ref = sub.u_ref_ev[1]
#         hlines!(ax, [ref]; color = :black, linestyle = :dash, label = "CoQui cRPA")
#         axislegend(ax; position = :rb)
#     end
#     Label(fig1[0, :], "Bare-channel U vs integration tolerance (level 0, monolayer graphene)";
#         fontsize = 18, font = :bold)
#     save(joinpath(FIG_DIR, "monolayer_tol_sweep_u_vs_tol.png"), fig1; px_per_unit = 2)

#     # --- Figure 2: relative error % vs tol, all channels overlaid ---
#     fig2 = Figure(; size = (900, 550))
#     ax2 = Axis(fig2[1, 1];
#         title = "Relative error vs CoQui (level 0)",
#         xlabel = "source_tol = volume_tol",
#         ylabel = "rel. error (%)",
#         xscale = log10,
#         xreversed = true,
#     )
#     for ch in CHANNELS
#         sub = df[df.channel .== ch, :]
#         sort!(sub, :tol, rev = true)
#         scatterlines!(ax2, sub.tol, sub.rel_err_pct;
#             color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
#             label = CHANNEL_LABEL[ch])
#     end
#     hlines!(ax2, [0.0]; color = :black, linestyle = :dash)
#     axislegend(ax2; position = :rb)
#     save(joinpath(FIG_DIR, "monolayer_tol_sweep_relerr_vs_tol.png"), fig2; px_per_unit = 2)

#     # --- Figure 3: self-convergence — |U(tol) - U(tol_min)| (log-log) ---
#     fig3 = Figure(; size = (900, 550))
#     ax3 = Axis(fig3[1, 1];
#         title = "Self-convergence: |U(tol) − U(tol_min)|",
#         xlabel = "source_tol = volume_tol",
#         ylabel = "|ΔU| from tightest tol (eV)",
#         xscale = log10,
#         yscale = log10,
#         xreversed = true,
#     )
#     for ch in CHANNELS
#         sub = df[df.channel .== ch, :]
#         sort!(sub, :tol, rev = true)
#         u_min = sub.u_ev[end]                   # value at tightest tol
#         deltas = abs.(sub.u_ev .- u_min)
#         keep = deltas .> 0                      # log scale: drop the zero point
#         scatterlines!(ax3, sub.tol[keep], deltas[keep];
#             color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
#             label = CHANNEL_LABEL[ch])
#     end
#     axislegend(ax3; position = :lb)
#     save(joinpath(FIG_DIR, "monolayer_tol_sweep_self_convergence.png"), fig3; px_per_unit = 2)

#     # --- Figure 4: all 4 channel energies on one axis, refs as vlines ---
#     fig4 = Figure(; size = (950, 600))
#     ax4 = Axis(fig4[1, 1];
#         title = "Monolayer bare U: all channels vs CoQui reference",
#         xlabel = "source_tol = volume_tol",
#         ylabel = "U (eV)",
#         xscale = log10,
#         xreversed = true,
#     )
#     for ch in CHANNELS
#         sub = df[df.channel .== ch, :]
#         sort!(sub, :tol, rev = true)
#         scatterlines!(ax4, sub.tol, sub.u_ev;
#             color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
#             label = CHANNEL_LABEL[ch])
#         ref = sub.u_ref_ev[1]
#         hlines!(ax4, [ref]; color = CHANNEL_COLOR[ch], linestyle = :dash, linewidth = 1.5)
#     end
#     axislegend(ax4; position = :rb)
#     save(joinpath(FIG_DIR, "monolayer_tol_sweep_all_channels.png"), fig4; px_per_unit = 2)

#     println("Wrote 4 figures to ", FIG_DIR)
# end

# main()


begin
    df = CSV.read(CSV_PATH, DataFrame)
    sort!(df, [:channel, :tol], rev = [false, true])  # tol descending so we sweep tighter

    # --- Figure 4: all 4 channel energies on one axis, refs as vlines ---
    fig4 = Figure(; size = (950, 500))
    ax4 = Axis(fig4[1, 1];
        title = "Monolayer bare U: all channels vs CoQui reference",
        xlabel = "source_tol = volume_tol",
        ylabel = "U (eV)",
        xscale = log10,
        xreversed = true,
    )
    ax5 = Axis(fig4[1, 2];
        title = "Monolayer bare U: all channels vs CoQui reference",
        xlabel = "source_tol = volume_tol",
        ylabel = "relative error vs CoQui (%)",
        xscale = log10,
        xreversed = true,
    )
    for ch in CHANNELS[1:2]
        sub = df[df.channel .== ch, :]
        sort!(sub, :tol, rev = true)
        scatterlines!(ax4, sub.tol, sub.u_ev;
            color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
            label = CHANNEL_LABEL[ch])
        ref = sub.u_ref_ev[1]
        ref_prl = PRL_REF_EV[Symbol(ch)]
        hlines!(ax4, [ref]; color = CHANNEL_COLOR[ch], linestyle = :dash, linewidth = 1.5, label = "CoQui cRPA")
        hlines!(ax4, [ref_prl]; color = CHANNEL_COLOR[ch], linewidth = 1.5, label = "PRL reference")
    end
    axislegend(ax4; position = :rb)
    ylims!(ax4, 0, 20)

    for ch in CHANNELS[3:4]
        sub = df[df.channel .== ch, :]
        sort!(sub, :tol, rev = true)
        scatterlines!(ax5, sub.tol, sub.u_ev;
            color = CHANNEL_COLOR[ch], markersize = 11, linewidth = 2,
            label = CHANNEL_LABEL[ch])
        ref = sub.u_ref_ev[1]
        hlines!(ax5, [ref]; color = CHANNEL_COLOR[ch], linestyle = :dash, label = "CoQui cRPA")
    end
    axislegend(ax5; position = :rb)
    ylims!(ax5, 0, 0.2)

    save(joinpath(FIG_DIR, "monolayer_tol_sweep_all_channels.png"), fig4; px_per_unit = 2)

    println("Wrote 4 figures to ", FIG_DIR)
end