# 6.3 analysis: per-eps2 error vs reference, iteration counts, screening curve.
# Reads data/raw/ratio_{test,ref}_e2_*.jls + v_vacuum.jls, prints the table,
# writes data/contrast_errors.csv.
# Cross-check: eps2 = 12 reference must reproduce the 6.1 fig1 reference
# V = 0.055678122198 (same protocol, same mesh construction).

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const DATA = joinpath(@__DIR__, "..", "data")
const RAW = joinpath(DATA, "raw")
const V61_REF = 0.055678122198

e2s = sort([parse(Float64, match(r"^ratio_test_e2_(.+)\.jls$", f).captures[1])
            for f in readdir(RAW) if occursin(r"^ratio_test_e2_.*\.jls$", f)])
Vvac = load_ref(joinpath(RAW, "v_vacuum.jls")).Vvac
@printf("V_vacuum = %.12e   (%d cases)\n\n", Vvac, length(e2s))

println("   eps2    gamma12       N_test  it_t  it_r        V_test         V_ref    V_relerr     V/V_vac   wall_t   wall_r")
rows = NamedTuple[]
for e2 in e2s
    t = load_ref(joinpath(RAW, "ratio_test_e2_$(e2).jls"))
    rf = isfile(joinpath(RAW, "ratio_ref_e2_$(e2).jls")) ?
         load_ref(joinpath(RAW, "ratio_ref_e2_$(e2).jls")) : nothing
    g = (e2 - 4.0) / (e2 + 4.0)
    relerr = rf === nothing ? NaN : abs(t.V - rf.V) / abs(rf.V)
    @printf("%7.4g  %9.6f  %11d  %4d  %4s  %12.6e  %12s  %10.3e  %10.4e  %7.0f  %7s\n",
            e2, g, t.N, t.niter,
            rf === nothing ? "-" : string(rf.niter),
            t.V, rf === nothing ? "-" : @sprintf("%.6e", rf.V),
            relerr, t.V / Vvac, t.wall,
            rf === nothing ? "-" : @sprintf("%.0f", rf.wall))
    push!(rows, (; eps2 = e2, gamma12 = g, N_test = t.N,
                 niter_test = t.niter, niter_ref = rf === nothing ? -1 : rf.niter,
                 V_test = t.V, V_ref = rf === nothing ? NaN : rf.V,
                 V_relerr = relerr, V_over_Vvac = t.V / Vvac,
                 wall_test = t.wall, wall_ref = rf === nothing ? NaN : rf.wall))
end

# cross-check against the 6.1 fig1 reference (identical protocol at eps2 = 12)
if any(r -> r.eps2 == 12.0 && !isnan(r.V_ref), rows)
    V12 = first(r for r in rows if r.eps2 == 12.0).V_ref
    @printf("\ncross-check eps2=12 ref vs 6.1 fig1 ref: %.12e vs %.12e  (rel dev %.2e)\n",
            V12, V61_REF, abs(V12 - V61_REF) / V61_REF)
end

csv = joinpath(DATA, "contrast_errors.csv")
isfile(csv) && rm(csv)
foreach(r -> append_csv_row(csv, r), rows)
println("\nwrote $(csv)")
