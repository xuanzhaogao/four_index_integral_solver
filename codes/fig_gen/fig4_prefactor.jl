#=
Diagnostic on the *absolute* near-field quadrature error
    E_abs(p, d) = ‖I_std − I_ref‖
together with the prefactor
    C_abs(p, d) = E_abs(p, d) · ρ(d)^(2p),   ρ(d) = d + √(1+d²).

Left  panel : absolute error vs d for each p in fig4_data.jls.
Right panel : the same data divided by the canonical Bernstein exponential
              ρ^(-2p). If the canonical scaling held with a d-independent
              prefactor, the right-panel curves would be flat in d.

Reading C_abs against C from the relative-error diagnostic isolates the
contribution of the normalization ‖I_ref‖ ~ 1/d in the far-field regime.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

const datapath = joinpath(@__DIR__, "fig4_data.jls")
const data     = open(deserialize, datapath, "r")

const d_list   = data.d_list
const p_values = sort(collect(data.p_values))
const sweeps   = Dict(s.p => s for s in data.sweeps)

ρ_of(d) = d + sqrt(1 + d^2)

# ---------------------------------------------------------------------------
# Print table and empirical slope of C_abs vs d (log-log)
# ---------------------------------------------------------------------------
function loglog_slope(xs, ys; mask = trues(length(xs)))
    n = count(mask)
    n < 2 && return NaN
    lx = log10.(xs[mask]); ly = log10.(ys[mask])
    x̄ = sum(lx)/n; ȳ = sum(ly)/n
    return sum((lx .- x̄) .* (ly .- ȳ)) / sum((lx .- x̄).^2)
end

println("\n=== absolute error E_abs(p, d) ===")
println(@sprintf("%-9s %-12s %-12s %-12s", "d", "E_abs(p=4)", "E_abs(p=6)", "E_abs(p=8)"))
for (k, d) in enumerate(d_list)
    es = [sweeps[p].E_abs[k] for p in p_values]
    @printf("%-9.4f  %-12.3e %-12.3e %-12.3e\n", d, es[1], es[2], es[3])
end

println("\n=== prefactor C_abs(p, d) = E_abs(p, d) · ρ(d)^(2p) ===")
println(@sprintf("%-9s %-12s %-12s %-12s", "d", "C(p=4)", "C(p=6)", "C(p=8)"))
for (k, d) in enumerate(d_list)
    cs = [sweeps[p].E_abs[k] * ρ_of(d)^(2p) for p in p_values]
    @printf("%-9.4f  %-12.3e %-12.3e %-12.3e\n", d, cs[1], cs[2], cs[3])
end

println("\n=== empirical d-slope α(p) of log C_abs(p, d) vs log d ===")
println(@sprintf("%-6s  %-14s  %-18s  %s", "p",
                 "α(p), all d", "α(p), d > 1.5", "α(p), well-resolved"))
for p in p_values
    sw   = sweeps[p]
    Cd   = sw.E_abs .* ρ_of.(d_list).^(2p)
    α_all = loglog_slope(d_list, Cd)
    α_far = loglog_slope(d_list, Cd; mask = d_list .> 1.5)
    # Well-resolved: above ~10× the round-off / HCubature floor
    floor_p = max(sw.R_ref[end] * 1e-13, 1e-30)
    α_ok = loglog_slope(d_list, Cd;
                        mask = (sw.E_abs .> 10 * floor_p) .& isfinite.(Cd))
    @printf("%-6d  %+10.3f      %+10.3f          %+10.3f\n",
            p, α_all, α_far, α_ok)
end

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig = Figure(size = (1100, 460), fontsize = 18)

palette = cgrad(:viridis, length(p_values) + 1, categorical = true)
markers = [:circle, :rect, :utriangle]

# Panel (a): absolute error vs d --------------------------------------------
ax_a = Axis(fig[1, 1];
            xscale = log10, yscale = log10,
            xlabel = L"d",
            ylabel = L"\Vert I_p - I_{\mathrm{ref}} \Vert_2")
for (i, p) in enumerate(p_values)
    sw = sweeps[p]
    y  = clamp.(sw.E_abs, 1e-30, 1e30)
    scatterlines!(ax_a, d_list, y;
                  color = palette[i], marker = markers[i],
                  markersize = 11, linewidth = 2,
                  label = L"p = %$p")
end
axislegend(ax_a; position = :rt)
xlims!(ax_a, 10^(-1.1), 10^(1.1))

# Panel (b): absolute error · ρ^(2p) ----------------------------------------
ax_b = Axis(fig[1, 2];
            xscale = log10, yscale = log10,
            xlabel = L"d",
            ylabel = L"\Vert I_p - I_{\mathrm{ref}} \Vert_2 \cdot \rho^{2p}")
for (i, p) in enumerate(p_values)
    sw  = sweeps[p]
    Cd  = sw.E_abs .* ρ_of.(d_list).^(2p)
    Cd  = clamp.(Cd, 1e-30, 1e30)
    scatterlines!(ax_b, d_list, Cd;
                  color = palette[i], marker = markers[i],
                  markersize = 11, linewidth = 2,
                  label = L"p = %$p")
end
axislegend(ax_b; position = :lt)
xlims!(ax_b, 10^(-1.1), 10^(1.1))

colgap!(fig.layout, 1, 28)

outpath = joinpath(@__DIR__, "figs/fig4_prefactor.pdf")
save(outpath, fig; px_per_unit = 4)
png_out = replace(outpath, ".pdf" => ".png")
save(png_out, fig; px_per_unit = 4)
@info "Saved diagnostic" outpath png_out

fig
