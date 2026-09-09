# Validation: does a pipeline change reproduce the tensor?
#
# Batching alters which pairs share an interface, the triangle change alters which entries are
# evaluated rather than mirrored, and the pottrg fix alters only threading. None of them touch
# the mathematics, so V_full_eV must agree with the reference to solver tolerance
# (gmres_rtol = 1e-5 here). A larger difference is a bug, not a tolerance effect.
#
# max_order is the exception: it changes the near-field quadrature, so a comparison ACROSS a
# max_order change is measuring a real accuracy shift, not checking an invariance. Read the
# verdict accordingly -- it only means "the two agree", which is the wrong question there.
#
# Run: julia --project=codes/lattice_scale scripts/compare_kb31.jl [ref_root] [new_root]
#   K-independence of the max_order = 64 tensor:
#     ... compare_kb31.jl lattice_conv_l3_eps2.4_kb31 lattice_conv_l3_eps2.4_k46
#   effect of the max_order fix (default):
#     ... compare_kb31.jl
#   full-matrix vs symmetry-restricted evaluation, same root:
#     ... compare_kb31.jl lattice_conv_l3_eps2.4_k46/V_full_eV_triangle.jls \
#         lattice_conv_l3_eps2.4_k46/V_full_eV.jls
using Serialization, Printf, Statistics, LinearAlgebra
const CEPH = "/mnt/ceph/users/xgao1/four_index"
const REF_ROOT = length(ARGS) >= 1 ? ARGS[1] : "lattice_conv_l3_eps2.4"
const NEW_ROOT = length(ARGS) >= 2 ? ARGS[2] : "lattice_conv_l3_eps2.4_kb31"
println("ref = $(REF_ROOT)\nnew = $(NEW_ROOT)\n")
# An argument may be a campaign root name (whose V_full_eV.jls is used) or a path to a .jls
# directly -- the latter is needed to compare two tensors sitting in the SAME root, e.g.
# V_full_eV.jls against a V_full_eV_triangle.jls kept from a previous evaluation.
_resolve(a) = endswith(a, ".jls") ? (isabspath(a) ? a : joinpath(CEPH, a)) :
                                    joinpath(CEPH, a, "V_full_eV.jls")
ref = deserialize(_resolve(REF_ROOT))
new = deserialize(_resolve(NEW_ROOT))

# The two campaigns hold the SAME pairs in a different order: consolidate writes pair_ids in
# manifest batch order, and k_target batching partitions the pairs differently from the old
# per-anchor grouping. Align by pair before comparing rather than assuming a common ordering.
Set(ref.pair_ids) == Set(new.pair_ids) ||
    error("pair SETS differ, not just their order: $(length(symdiff(Set(ref.pair_ids), Set(new.pair_ids)))) pairs unmatched")
perm = indexin(ref.pair_ids, new.pair_ids)
any(isnothing, perm) && error("unmatched pairs after alignment")
A, B = ref.V, new.V[perm, perm]
d = abs.(A .- B); scale = maximum(abs, A)
@printf("pairs               %d\n", length(ref.pair_ids))
@printf("max|V_ref|          %.6g eV\n", scale)
@printf("max abs difference  %.3e eV\n", maximum(d))
@printf("max rel difference  %.3e   (vs gmres_rtol 1e-5)\n", maximum(d) / scale)
@printf("rms rel difference  %.3e\n", sqrt(mean(d .^ 2)) / scale)

# the on-site U values are what Sec. 5.3 quotes, so check them explicitly
ons = [i for (i, p) in enumerate(ref.pair_ids) if p[1] == p[2]]
ua = [A[i, i] for i in ons]; ub = [B[i, i] for i in ons]
@printf("\non-site U (%d orbitals)\n", length(ons))
@printf("  ref  min %.5f  max %.5f  mean %.5f eV\n", minimum(ua), maximum(ua), mean(ua))
@printf("  new  min %.5f  max %.5f  mean %.5f eV\n", minimum(ub), maximum(ub), mean(ub))
@printf("  max |dU|  %.3e eV   (%.2e relative)\n",
        maximum(abs.(ua .- ub)), maximum(abs.(ua .- ub) ./ abs.(ua)))
verdict = maximum(d) / scale < 1e-4 ? "PASS" : "FAIL -- investigate before trusting the rerun"
println("\nverdict: ", verdict)
