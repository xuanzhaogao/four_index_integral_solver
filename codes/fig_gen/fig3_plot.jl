#=
Figure 3 plot script. Loads pre-computed data from fig3_data.jls and
generates fig3_edge_singularity.png.

Layout:
  (a) 3D surface plot of σ on the +z face of the cube (deepest edge-local mesh).
      σ value is the height field, also colored by σ.
  (b) Self-convergence of the far-field single-layer potential. For each
      refinement level we evaluate u_ℓ at a set of targets on a sphere of
      radius R_far outside the cube; the deepest edge-local level serves as
      reference, and we plot ||u_ℓ − u_ref||_2 / ||u_ref||_2 vs N.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

const datapath = joinpath(@__DIR__, "fig3_data.jls")
const data = open(deserialize, datapath, "r")

begin
    fig = Figure(size = (1300, 540), fontsize = 18)

    # ----- Panel (a): σ on +z face as a 3D surface ------------------------
    ax_a = Axis3(fig[1, 1];
                 aspect = (1, 1, 0.6),
                 xlabel = "x", ylabel = "y", zlabel = L"|\sigma|",
                 azimuth = 0.30π, elevation = 0.18π,
                 protrusions = (40, 20, 25, 25))
    surface!(ax_a, data.surface_xs, data.surface_ys, abs.(data.surface_sigma);
             colormap = :viridis, shading = NoShading)

    # ----- Panel (b): far-field potential self-convergence ----------------
    ax_b = Axis(fig[1, 2];
                xscale = log10, yscale = log10,
                xlabel = "DOF",
                ylabel = L"\mathcal{E}_u")

    Ns_e = [r.N for r in data.edge_results]
    Ns_u = [r.N for r in data.uniform_results]
    e_e  = data.edge_pot_errors
    e_u  = data.uniform_pot_errors

    keep_e = isfinite.(e_e) .& (e_e .> 0)
    keep_u = isfinite.(e_u) .& (e_u .> 0)

    scatterlines!(ax_b, Ns_e[keep_e], e_e[keep_e]; color = :crimson,
                  marker = :circle, markersize = 12, linewidth = 2,
                  label = "edge-local")
    scatterlines!(ax_b, Ns_u[keep_u], e_u[keep_u]; color = :royalblue,
                  marker = :rect, markersize = 12, linewidth = 2,
                  label = "uniform")
    axislegend(ax_b; position = :rt)

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "fig3_edge_singularity.png")
    save(outpath, fig; px_per_unit = 4)
    @info "Saved figure" outpath
end
