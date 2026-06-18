# Re-normalize the eV four-index tensor by ORBITAL norms instead of per-pair density norms.
# assemble_v divided V_raw by Na*Nb with Na = ∫ρ_a (the pair density norm). That is correct
# only for onsite pairs (∫ρ_ii = ‖φ_i‖²); for off-diagonal pairs ∫ρ_ij = the orbital overlap,
# which → 0 for distant pairs and blows up the eV value. The physical four-index element uses
# orbital-normalized densities ρ̂_ij = φ_iφ_j/(‖φ_i‖‖φ_j‖), i.e. normalize by ‖φ_i‖‖φ_j‖.
# Recovers V_raw from the stored (wrong) eV + per-pair norms, re-applies the orbital-norm
# factor (orbital norms = the onsite-pair norms), and overwrites V_full_eV.{jls,tsv}.
using Serialization, LinearAlgebra
root = "/mnt/ceph/users/xgao1/four_index/lattice_10x10_multicube"
d = deserialize(joinpath(root, "V_full_eV.jls"))
pid, VeV, Na, E2 = d.pair_ids, d.V, d.norms, d.E2_eV_Ang
n = length(pid)

# orbital norm^2 ‖φ_i‖² = the onsite pair (i,i) density norm
orbn2 = Dict{Int,Float64}()
for (a, (i, j)) in enumerate(pid); i == j && (orbn2[i] = Na[a]); end
@assert length(orbn2) == length(unique(vcat(first.(pid), last.(pid)))) "missing onsite pair for some orbital"
orbnrm = [sqrt(orbn2[pid[a][1]] * orbn2[pid[a][2]]) for a in 1:n]   # ‖φ_i‖‖φ_j‖ for pair a=(i,j)

# Vfix[a,b] = V_raw[a,b] * 4π E2 / (orbnrm[a] orbnrm[b]); recover V_raw = VeV*(Na[a]Na[b])/(4π E2)
Vfix = Matrix{Float64}(undef, n, n)
for a in 1:n, b in 1:n
    Vfix[a, b] = VeV[a, b] * (Na[a] * Na[b]) / (orbnrm[a] * orbnrm[b])
end

onn = sqrt(sum(values(orbn2)) / length(orbn2))          # rms orbital norm (≈‖φ‖)
asym = maximum(abs.(Vfix .- transpose(Vfix))) / maximum(abs, Vfix)
println("orbitals=", length(orbn2), "  ‖φ‖²≈", round(sum(values(orbn2))/length(orbn2); digits=3),
        "  min per-pair Na=", minimum(Na))
println("BEFORE: max|V_eV|=", maximum(abs, VeV), "   AFTER: max|V_eV|=", maximum(abs, Vfix),
        "  onsite[1,1]=", round(Vfix[1,1]; digits=4), " eV  max_rel_asym=", round(asym; digits=4))

serialize(joinpath(root, "V_full_eV.jls"),
    (; pair_ids = pid, V = Vfix, orbital_norms = orbn2, unit = "eV", E2_eV_Ang = E2, normalization = "orbital"))
open(joinpath(root, "V_full_eV.tsv"), "w") do io
    println(io, "i\tj\tk\tl\tV")
    for r in 1:n, cc in 1:n
        (i, j) = pid[r]; (k, l) = pid[cc]
        println(io, i, '\t', j, '\t', k, '\t', l, '\t', Vfix[r, cc])
    end
end
println("rewrote V_full_eV.{jls,tsv}  (orbital-normalized eV)")
