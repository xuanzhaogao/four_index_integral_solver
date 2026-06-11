# 6.3 analysis: per-eps2 convergence tables (p = 6, r = 1..5) against each
# eps2's own reference (p = 8, eps = 1e-6, r = 6). Prints accuracy / N_iter /
# DOF per run, writes data/contrast_errors.csv.

include(joinpath(@__DIR__, "..", "..", "common", "Lite.jl"))
using .Lite
using Printf

const DATA = joinpath(@__DIR__, "..", "data")
const RAW = joinpath(DATA, "raw")

e2s = sort([parse(Float64, match(r"^ratio_ref_e2_(.+)\.jls$", f).captures[1])
            for f in readdir(RAW) if occursin(r"^ratio_ref_e2_.*\.jls$", f)])
rs_of(e2) = sort([parse(Int, match(r"_r(\d+)\.jls$", f).captures[1])
                  for f in readdir(RAW) if startswith(f, "ratio_test_e2_$(e2)_r")])
Vvac = load_ref(joinpath(RAW, "v_vacuum.jls")).Vvac
@printf("V_vacuum = %.12e\n", Vvac)

rows = NamedTuple[]
for e2 in e2s
    rf = load_ref(joinpath(RAW, "ratio_ref_e2_$(e2).jls"))
    g = (e2 - 4.0) / (e2 + 4.0)
    @printf("\n=== eps2 = %g   (gamma12 = %+.4f)   V_ref = %.12e   (N = %d, %d iters, %.0fs)\n",
            e2, g, rf.V, rf.N, rf.niter, rf.wall)
    println("   r          N  niter      V_relerr        wall")
    for r in rs_of(e2)
        t = load_ref(joinpath(RAW, "ratio_test_e2_$(e2)_r$(r).jls"))
        relerr = abs(t.V - rf.V) / abs(rf.V)
        @printf("  %2d  %9d  %5d    %10.3e    %6.0fs\n", r, t.N, t.niter, relerr, t.wall)
        push!(rows, (; eps2 = e2, gamma12 = g, r, N = t.N,
                     niter = t.niter, V = t.V, V_ref = rf.V, V_relerr = relerr,
                     niter_ref = rf.niter, N_ref = rf.N,
                     V_over_Vvac = t.V / Vvac, wall = t.wall))
    end
end

csv = joinpath(DATA, "contrast_errors.csv")
isfile(csv) && rm(csv)
foreach(r -> append_csv_row(csv, r), rows)
println("\nwrote $(csv)")
