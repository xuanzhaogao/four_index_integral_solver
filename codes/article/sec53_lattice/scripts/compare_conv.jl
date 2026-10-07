# Compare the two l_ec self-convergence runs (lattice_conv_l2 vs lattice_conv_l3): on-site U
# (diagonal of V) matched by orbital, and each run's exchange-symmetry bound. Reports max/mean
# |ΔU| and relative change, plots U(x) for both l_ec, and prints the symmetry from report.txt.
# Run: julia --project=codes/article scripts/compare_conv.jl
using Serialization, Printf, CairoMakie
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
const CEPH = "/mnt/ceph/users/xgao1/four_index"

# orbital → (x, U_eV) from the diagonal of V_full_eV.jls
function onsite(name)
    d = open(deserialize, joinpath(CEPH, name, "V_full_eV.jls"))
    pos_of = Dict{Int,Float64}()      # orbital → x, read from campaign toml
    for ln in eachline(joinpath(@__DIR__, "..", "campaigns", "$(name).toml"))
        # not needed: positions taken from V diag via pair (i,i); x from the toml orbital order
    end
    U = Dict{Int,Float64}()
    for (a, p) in enumerate(d.pair_ids); p[1] == p[2] && (U[p[1]] = d.V[a, a]); end
    return U
end

# orbital x positions from a campaign toml (order of [[orbital]] blocks = orbital id)
function xpos(name)
    xs = Float64[]
    for ln in eachline(joinpath(@__DIR__, "..", "campaigns", "$(name).toml"))
        m = match(r"^x = (.+)$", ln); m !== nothing && push!(xs, parse(Float64, m.captures[1]))
    end
    return xs
end

U2, U3 = onsite("lattice_conv_l2"), onsite("lattice_conv_l3")
X = xpos("lattice_conv_l2")
ids = sort(collect(intersect(keys(U2), keys(U3))))
isempty(ids) && error("no matching orbitals — did both runs finish assemble?")

u2 = [U2[i] for i in ids]; u3 = [U3[i] for i in ids]; xs = [X[i] for i in ids]
del = abs.(u2 .- u3); rel = del ./ abs.(u2)
@printf("n=%d orbitals.  l_ec 2.27 vs 1.14:  max|ΔU|=%.2e eV  mean=%.2e  max rel=%.2e (3-digit target=1e-3)\n",
        length(ids), maximum(del), sum(del)/length(del), maximum(rel))
for nm in ("lattice_conv_l2", "lattice_conv_l3")
    rp = joinpath(CEPH, nm, "report.txt")
    isfile(rp) && println("  $nm: ", strip(join(filter(l -> occursin("asymmetry", l), readlines(rp)), "")))
end

o = sortperm(xs)
fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1]; xlabel = "orbital x (Å)", ylabel = "onsite U (eV)",
          title = "l_ec self-convergence (all tol = 1e-5)")
scatterlines!(ax, xs[o], u2[o]; color = QUAL.blue, marker = :circle, markersize = MS, linewidth = LW_DATA, label = "l_ec = 2.27")
scatterlines!(ax, xs[o], u3[o]; color = QUAL.orange, marker = :xcross, markersize = MS, linewidth = LW_DATA, label = "l_ec = 1.14")
axislegend(ax; position = :lt, framevisible = true)
out = joinpath(@__DIR__, "..", "figs", "fig_lec_convergence.pdf")
save(out, fig; px_per_unit = PX_PER_UNIT)
println("wrote $out")
