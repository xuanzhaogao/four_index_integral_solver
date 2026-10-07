# Analysis of the 6.1 Fig.-1 system sweep: errors vs the reference run.
# Writes data/fig1_errors.csv and prints a summary table.

include(joinpath(@__DIR__, "..", "..", "..", "common", "Harness.jl"))
using .Harness
using LinearAlgebra, Printf

# RERUN_TAG: read the rerun output tree instead of the published one.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "data" * TAG)

ref = load_ref(joinpath(DATA, "raw", "fig1_ref.jls"))
@printf("reference: p=%d r=%d eps=%.0e N=%d niter=%d V=%.12e\n\n",
        ref.p, ref.r, ref.eps, ref.N, ref.niter, ref.V)

rel2(a, b) = norm(a .- b) / norm(b)

rows = NamedTuple[]
@printf("%-3s %-3s %9s %6s | %10s | %10s %10s %10s | %8s\n",
        "p", "r", "N", "niter", "V_relerr", "phi_near", "phi_supp", "phi_far", "wall[s]")
for p in (2, 4, 6), r in 1:5
    f = joinpath(DATA, "raw", "fig1_p$(p)_r$(r).jls")
    isfile(f) || continue
    d = load_ref(f)
    ev = abs(d.V - ref.V) / abs(ref.V)
    en = rel2(d.phi_near, ref.phi_near)
    es = rel2(d.phi_supp, ref.phi_supp)
    ef = rel2(d.phi_far, ref.phi_far)
    @printf("%-3d %-3d %9d %6d | %10.3e | %10.3e %10.3e %10.3e | %8.1f\n",
            p, r, d.N, d.niter, ev, en, es, ef, d.wall)
    push!(rows, (; p, r, N = d.N, niter = d.niter, V_relerr = ev,
                 phi_near = en, phi_supp = es, phi_far = ef, V = d.V, wall = d.wall))
end
@printf("\nreference timings: %s\n",
        join(["$k=$(round(v; digits=1))s" for (k, v) in sort(collect(ref.times))], "  "))

out = joinpath(DATA, "fig1_errors.csv")
isfile(out) && rm(out)
for row in rows
    append_csv_row(out, row)
end
println("wrote $out")
