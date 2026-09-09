# Build and cache the FIXED evaluation target set for the Sec. 5.2 evaluation benchmark:
# the union of the quadrature points of ALL orbitals of the Sec. 5.3 model (198 of them).
#
# Why a fixed set, and why cached:
#   * Fixed -- the earlier block-eval used the K-batch's OWN grid as targets, which grows with
#     K (46,878 at K = 1 to 634,750 at K = 46). Time per source then mixes two effects and is
#     uninterpretable. With N_p held constant, t_eval/K measures exactly one thing.
#   * Cached -- assembling 198 orbitals costs minutes and tens of GB. Rebuilding it inside
#     every array task would dominate the benchmark it is supposed to support.
#
# Each orbital is snapped as its own onsite pair (i,i): that support is the widest any pair
# density involving orbital i can have, so their union covers every possible target.
#
# Writes targets_model.jls (positions + weights + provenance) next to the data tree.
# Run:  RERUN_TAG=_sec53evalpilot julia --project=codes/numerical_results \
#         exp66_multicube/scripts/prep_model_targets.jl
using BoundaryIntegral, TOML, Serialization, Printf, Dates
import BoundaryIntegral as BI

const TAG      = get(ENV, "RERUN_TAG", "")
const EXP      = normpath(joinpath(@__DIR__, ".."))
const CAMPAIGN = get(ENV, "MULTICUBE_TARGET_CAMPAIGN",
    "/mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale/campaigns/lattice_conv_l3_eps2.4.toml")
const OUT = get(ENV, "MULTICUBE_TARGET_CACHE", joinpath(EXP, "data" * TAG, "targets_model.jls"))
mkpath(dirname(OUT))

if isfile(OUT)
    d = deserialize(OUT)
    @printf("cache already present: %s\n  N_p = %d  (from %s)\n", OUT, size(d.positions, 2), d.campaign)
    exit(0)
end

c = TOML.parsefile(CAMPAIGN)
TEMPL = [BI.read_xsf(t)[2] for t in c["templates"]]
CEN = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
orbs = c["orbital"]
@printf("campaign %s\n  %d orbitals\n", basename(CAMPAIGN), length(orbs))
flush(stdout)

insts = Dict{Int, BI.OrbitalInstance}()
for (i, o) in enumerate(orbs)
    t = o["type"]
    insts[i] = BI.OrbitalInstance(i, t, BI.snap_orbital(TEMPL[t], CEN[t], (o["x"], o["y"], o["z"])))
end

t0 = time()
b = BI.assemble_lattice_batch(TEMPL, insts, [(i, i) for i in 1:length(orbs)]; support_rtol = 1e-4)
t_asm = time() - t0
src = BI.batch_volume_sources(b)[1]        # shared-positions contract: one grid for all
positions = Matrix{Float64}(src.positions)
weights = Vector{Float64}(src.weights)
N_p = size(positions, 2)
@printf("assembled in %.1f s   N_p = %d   maxrss %.1f GB\n", t_asm, N_p, Sys.maxrss() / 2^30)
@printf("extent: x %.2f..%.2f  y %.2f..%.2f  z %.2f..%.2f\n",
        extrema(positions[1, :])..., extrema(positions[2, :])..., extrema(positions[3, :])...)

# Only positions and weights are cached. The 198 target-side densities would be another
# N_p x 198 matrix (GBs) and are not needed to TIME the field evaluation, which is what the
# benchmark measures; assembling the V entries themselves is a separate, cheap contraction.
serialize(OUT, (; positions, weights, N_p, campaign = CAMPAIGN,
                  n_orbitals = length(orbs), support_rtol = 1e-4,
                  built = string(Dates.now())))
@printf("wrote %s  (%.1f MB)\n", OUT, filesize(OUT) / 2^20)
