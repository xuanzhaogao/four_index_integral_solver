# Extract the on-site Hubbard U (diagonal V[ρ_ii, ρ_ii]) from an assembled four-index tensor
# V_full_eV.jls, report its spatial variation across the Si|SiO₂ junction, and plot U(x).
# Run: julia --project=codes/lattice_scale scripts/extract_onsite_from_V.jl [campaign.toml]
using BoundaryIntegral, Serialization, LinearAlgebra, Printf, CairoMakie
const BI = BoundaryIntegral
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

toml = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "campaigns", "lattice_10x10_het3x.toml")
c = BI.load_campaign(toml)
d = open(deserialize, joinpath(c.root, "V_full_eV.jls"))
pair_ids, V = d.pair_ids, d.V
XJUNC = 5.547

recs = Tuple{Int,Float64,Float64,Float64}[]          # (orbital, x, y, U_eV)
for (a, p) in enumerate(pair_ids)
    p[1] == p[2] || continue
    pos = c.orbitals[p[1]].pos
    push!(recs, (p[1], pos[1], pos[2], V[a, a]))
end
sort!(recs, by = r -> r[2])
Us = [r[4] for r in recs]
si = [r for r in recs if r[2] < XJUNC]; ox = [r for r in recs if r[2] >= XJUNC]
@printf("onsite U: n=%d  min=%.4f  max=%.4f  mean=%.4f eV   (Si %d: mean %.4f | SiO2 %d: mean %.4f)\n",
        length(recs), minimum(Us), maximum(Us), sum(Us)/length(Us),
        length(si), sum(r->r[4], si)/length(si), length(ox), sum(r->r[4], ox)/length(ox))

open(joinpath(c.root, "onsite_U_diag.tsv"), "w") do io
    println(io, "orbital\tx\ty\tU_eV")
    for r in recs; @printf(io, "%d\t%.6g\t%.6g\t%.8g\n", r...); end
end

fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1]; xlabel = "orbital x (Å)   Si ← junction → SiO₂",
          ylabel = "onsite U  (eV)", title = "$(c.name): onsite Hubbard U across the junction")
vlines!(ax, [XJUNC]; color = (:gray, 0.6), linestyle = :dash, linewidth = LW_DATA)
scatter!(ax, [r[2] for r in si], [r[4] for r in si]; color = QUAL.blue, markersize = MS, label = "over Si (ε=11.9)")
scatter!(ax, [r[2] for r in ox], [r[4] for r in ox]; color = QUAL.orange, markersize = MS, label = "over SiO₂ (ε=3.9)")
axislegend(ax; position = :lt, framevisible = true)
out = joinpath(@__DIR__, "..", "figs", "fig_benchmark_het3x_onsiteU.pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
println("wrote $out")
