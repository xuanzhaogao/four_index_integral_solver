# Diagnostic for the multi-RHS n_src question: why does the union support grow
# with K, given supp(rho_1j) = supp(phi_1) ∩ supp(phi_j) ⊆ supp(phi_1)?
#
# Confirms: (1) the union is ALWAYS ⊆ supp(phi_1) (no point outside phi_1's
# footprint); (2) the growth is the rss-truncation revealing phi_1's tail —
# onsite rho_11 = phi_1^2 (quadratic) is kept only down to phi_1 ~ 1% of peak,
# while a neighbor pair rho_1j = phi_1*phi_j (linear in phi_1 near site j, where
# phi_j ~ its peak) is kept down to phi_1 ~ 0.01% of peak, so it lifts phi_1's
# faint tail above threshold. Bounded by |supp(phi_1)|.
#
# Geometry only (no solve). Run: julia -t 4 --project=<BI> diag_support.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Printf, Statistics

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const TEMPL = [BI.read_xsf(joinpath(REF_DIR, "graphene_00001.xsf"))[2],
               BI.read_xsf(joinpath(REF_DIR, "graphene_00002.xsf"))[2]]
const CENTROID = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)
addv(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
scalev(s, a) = (s * a[1], s * a[2], s * a[3])
distv(a, b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)

const dg1 = TEMPL[1]
const NX, NY, NZ = dg1.nx, dg1.ny, dg1.nz
phi1_at(g) = (1 <= g[1] <= NX && 1 <= g[2] <= NY && 1 <= g[3] <= NZ) ? dg1.values[g[1], g[2], g[3]] : 0.0
const PHI1_PEAK = maximum(abs, dg1.values)

function geom(cutoff; rtol = 1e-4, nrange = 3)
    center = CENTROID[1]
    sites = Tuple{Int, NTuple{3,Float64}}[(1, center)]
    for t in 1:2, n1 in -nrange:nrange, n2 in -nrange:nrange
        (t == 1 && n1 == 0 && n2 == 0) && continue
        pos = addv(CENTROID[t], addv(scalev(Float64(n1), A1), scalev(Float64(n2), A2)))
        distv(pos, center) <= cutoff && push!(sites, (t, pos))
    end
    insts = Dict(id => BI.OrbitalInstance(id, t, BI.snap_orbital(TEMPL[t], CENTROID[t], pos))
                 for (id, (t, pos)) in enumerate(sites))
    pairs = [(1, id) for id in 1:length(sites)]
    return BI.assemble_lattice_batch(TEMPL, insts, pairs; support_rtol = rtol), pairs
end

ngrid = NX * NY * NZ
nz_phi1 = count(!iszero, dg1.values)
@printf("grid total = %d   |supp(phi_1)| (phi_1 != 0) = %d  (%.2f%% of grid)  -- the hard bound on n_src\n\n",
        ngrid, nz_phi1, 100 * nz_phi1 / ngrid)

@printf("%-8s %-4s %-10s %-22s\n", "cutoff", "K", "n_src", "union points outside supp(phi_1)")
for cutoff in [0.0, 1.5, 2.5, 2.9, 3.9]
    b, pairs = geom(cutoff)
    outside = count(g -> phi1_at(g) == 0.0, b.gidx)
    @printf("%-8.1f %-4d %-10d %-22d\n", cutoff, length(pairs), length(b.gidx), outside)
end

println("\n-- where do the points ADDED going K=1 -> K=19 live in phi_1? --")
b1, _ = geom(0.0)
core = Set(b1.gidx)
b19, _ = geom(3.9)
rowof = Dict(g => r for (r, g) in enumerate(b19.gidx))
added = [g for g in b19.gidx if !(g in core)]
rel_phi1_core = [abs(phi1_at(g)) / PHI1_PEAK for g in b1.gidx]
rel_phi1_added = [abs(phi1_at(g)) / PHI1_PEAK for g in added]
ratios = Float64[]
for g in added
    d = @view b19.densities[rowof[g], :]
    r11 = abs(d[1])
    r11 > 0 && push!(ratios, maximum(abs, d) / r11)
end
@printf("onsite core (K=1): %d pts   phi_1/peak  min=%.2e  median=%.2e\n",
        length(b1.gidx), minimum(rel_phi1_core), median(rel_phi1_core))
@printf("added at K=19    : %d pts   phi_1/peak  min=%.2e  median=%.2e  max=%.2e\n",
        length(added), minimum(rel_phi1_added), median(rel_phi1_added), maximum(rel_phi1_added))
@printf("at added pts: max_k|rho_1k| / |rho_11|   median=%.1f  max=%.1f  (neighbor pair >> onsite)\n",
        median(ratios), maximum(ratios))
println("\nDIAG SUPPORT DONE")
