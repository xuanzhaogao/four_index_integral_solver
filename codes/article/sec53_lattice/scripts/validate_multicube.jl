# Quick validation of the box-based multi-region screening on a small multicube substrate
# (demo_2x2_multicube): runs the full four-index pipeline in-memory (no ceph files) and
# checks the V tensor is finite, sane, and ~symmetric. Confirms the BI change before the
# production 10x10 campaign commits nodes.
#   julia --project=codes/article codes/article/sec53_lattice/scripts/validate_multicube.jl
using BoundaryIntegral
using LinearAlgebra

c = load_campaign(joinpath(@__DIR__, "..", "campaigns", "demo_2x2_multicube.toml"))
println("campaign=", c.name, "  boxes=", length(c.boxes), "  epses=", c.epses, "  eps_out=", c.eps_out)
flush(stdout)

r = four_index_integrals(c)          # solve_batch_core + eval_batch_core (box-based screening)
V, pid = r.V, r.pair_ids
asym = maximum(abs.(V .- transpose(V))) / maximum(abs, V)
@info "validation" n_pairs=length(pid) V11_raw=V[1,1] finite=all(isfinite, V) maxabs=maximum(abs, V) max_rel_asym=asym
println("VALIDATION DONE  (n_pairs=$(length(pid)), V[1,1]_raw=$(V[1,1]), max_rel_asym=$(round(asym; digits=4)))")
