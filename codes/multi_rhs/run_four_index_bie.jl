# Run the full Step 0–7 pipeline from a .bie file and print the four-index integrals V for a
# center group. The operator A and per-source u_inc are built once; V is contracted from the
# block-GMRES layer densities (see four_index_integrals).
#
# Usage:
#   julia --project scripts/run_four_index_bie.jl [bie_path] [center_id]
# Defaults: scripts/graphene_8neighbors.bie, center 1.

using BoundaryIntegral
using Printf

bie_path  = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "graphene_8neighbors.bie")
center_id = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

println("Reading ", bie_path, "  (center ", center_id, ")")
si = read_system_input(bie_path)

println("group(", center_id, ") = ", sort(si.groups[center_id]))
for id in sort(collect(keys(si.orbitals)))
    o = si.orbitals[id]
    @printf("  orb %d  shift=%s  center=(%.3f, %.3f, %.3f)\n", id, o.grid_shift, o.center...)
end

println("Solving (Step 0–7) ...")
res = four_index_integrals(si, center_id)
K = length(res.labels)
println("  labels   = ", res.labels)
println("  interface: ", length(res.interface.panels), " panels, ",
        BoundaryIntegral.num_points(res.interface), " pts;  Sigma ", size(res.sigma))

println("\nV_ab (raw), rows/cols = ", join(res.labels, ", "), ":")
for a in 1:K
    @printf("  %-12s", res.labels[a])
    for b in 1:K
        @printf(" % .5g", res.V[a, b])
    end
    println()
end

# convenience: the first row (center-on-site against the whole group)
println("\nV[center, ·] (on-site + neighbors):")
for b in 1:K
    @printf("  %-12s % .6g\n", res.labels[b], res.V[1, b])
end
println("done.")
