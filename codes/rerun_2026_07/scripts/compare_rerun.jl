#=
Diff a rerun campaign against the published one.

This is what the paper needs in order to state the effect of the two LHS fixes:
how far the tensor moved, and how much the iteration count fell. It reads only
finished campaign outputs, so it cannot be wrong about BI internals.

Reports:
  * on-site U statistics both ways (min, max, per-side means) and the shift
  * max and RMS relative change over the whole 2623x2623 tensor
  * GMRES iterations per batch, published vs rerun
  * exchange-symmetry error both ways

Usage:
    BASE=lattice_conv_l3 TAG=_v2 \
    julia --project=codes/lattice_scale codes/rerun_2026_07/scripts/compare_rerun.jl
=#

using BoundaryIntegral, Serialization, Printf, Statistics
const BI = BoundaryIntegral

const CEPH = "/mnt/ceph/users/xgao1/four_index"
const BASE = get(ENV, "BASE", "lattice_conv_l3")
const TAG  = get(ENV, "TAG", "_v2")
const XJ   = 5.547   # junction plane, Angstrom

const CAMPAIGNS = normpath(joinpath(@__DIR__, "..", "..", "lattice_scale", "campaigns"))

pub_dir = joinpath(CEPH, BASE)
new_dir = joinpath(CEPH, BASE * TAG)
for d in (pub_dir, new_dir)
    isdir(d) || error("missing campaign directory: $d")
end

c = BI.load_campaign(joinpath(CAMPAIGNS, BASE * ".toml"))
px = [o.pos[1] for o in c.orbitals]
norb = length(c.orbitals)

function onsite(dir)
    d = open(deserialize, joinpath(dir, "V_full_eV.jls"))
    p = sortperm(d.pair_ids)
    U = fill(NaN, norb)
    for (a, pr) in enumerate(d.pair_ids)
        pr[1] == pr[2] && (U[pr[1]] = d.V[a, a])
    end
    return U, d.V[p, p]
end

@info "loading tensors" pub_dir new_dir
Upub, Vpub = onsite(pub_dir)
Unew, Vnew = onsite(new_dir)

ok = .!isnan.(Upub) .& .!isnan.(Unew)
up, un, x = Upub[ok], Unew[ok], px[ok]

stats(u, m) = (minimum(u[m]), maximum(u[m]), mean(u[m]))
si = x .< XJ; ox = x .>= XJ

println("\n=== on-site U (eV), $(BASE) ===")
@printf("%-26s %8s %8s %10s %10s\n", "", "min", "max", "Si mean", "SiO2 mean")
for (lbl, u) in (("published (both off)", up), ("rerun$(TAG) (both on)", un))
    @printf("%-26s %8.4f %8.4f %10.4f %10.4f\n",
            lbl, minimum(u), maximum(u), mean(u[si]), mean(u[ox]))
end
dU = un .- up
@printf("%-26s %8.4f %8.4f %10.4f %10.4f\n", "shift (rerun - published)",
        minimum(dU), maximum(dU), mean(dU[si]), mean(dU[ox]))
@printf("\nmax |dU|/U = %.3e   RMS |dU|/U = %.3e\n",
        maximum(abs.(dU ./ up)), sqrt(mean((dU ./ up) .^ 2)))
@printf("Si/SiO2 contrast:  published %.4f eV   rerun %.4f eV\n",
        mean(up[ox]) - mean(up[si]), mean(un[ox]) - mean(un[si]))

# whole-tensor difference, on the shared sparsity pattern
scale = maximum(abs.(Vpub))
D = abs.(Vnew .- Vpub) ./ scale
@printf("\n=== full tensor (%d x %d) ===\nmax |dV|/max|V| = %.3e   RMS = %.3e\n",
        size(Vpub, 1), size(Vpub, 2), maximum(D), sqrt(mean(D .^ 2)))

asym(V) = maximum(abs.(V .- transpose(V))) / maximum(abs.(V))
@printf("exchange symmetry error:  published %.3e   rerun %.3e\n", asym(Vpub), asym(Vnew))

# iteration counts, sampled across batches
println("\n=== GMRES iterations per batch ===")
function niters(dir)
    bd = joinpath(dir, "batches")
    isdir(bd) || return Int[]
    fs = sort(readdir(bd))
    [deserialize(joinpath(bd, f)).stats["niter"] for f in fs[1:max(1, cld(length(fs), 12)):end]]
end
np, nn = niters(pub_dir), niters(new_dir)
if !isempty(np) && !isempty(nn)
    @printf("published: min %d  median %.1f  max %d   (n=%d sampled)\n",
            minimum(np), median(np), maximum(np), length(np))
    @printf("rerun%s:   min %d  median %.1f  max %d   (n=%d sampled)\n",
            TAG, minimum(nn), median(nn), maximum(nn), length(nn))
    @printf("median speedup in iterations: %.1fx\n", median(np) / median(nn))
else
    println("(batch files unavailable on one side; skipped)")
end
