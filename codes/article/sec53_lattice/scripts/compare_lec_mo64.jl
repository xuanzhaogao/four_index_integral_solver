# The l_ec self-convergence sweep at max_order = 8 vs 64.
#
# The two questions are different and both matter:
#   1. Does U still converge in l_ec?  That is the figure's claim, and the max_order = 8 sweep
#      answered it -- flat to 4e-4 eV over seven levels. Resolving the near field must not break
#      that, and if it does the figure's claim was an artifact.
#   2. Where did the plateau move to?  max_order = 8 clamped 92.5% of the non-touching near
#      corrections below their own tolerance, at every level equally, so the sweep converged to
#      a biased value it could not see. The shift is the size of that bias.
#
# Run: julia --project=codes/article scripts/compare_lec_mo64.jl
using Printf
const FIGS = normpath(joinpath(@__DIR__, "..", "figs"))

read_tsv(f) = begin
    rows = Dict{Int,NTuple{4,Float64}}()
    for ln in eachline(f)
        (isempty(ln) || startswith(ln, "level")) && continue
        c = split(ln, '\t')
        rows[parse(Int, c[1])] = (parse(Float64, c[2]), parse(Float64, c[3]),
                                  parse(Float64, c[4]), parse(Float64, c[5]))
    end
    rows
end

const A = read_tsv(joinpath(FIGS, "lec_single_conv_eps2.4.tsv"))        # max_order = 8
const BF = joinpath(FIGS, "lec_single_conv_eps2.4_mo64.tsv")
isfile(BF) || error("$(BF) not written yet -- the max_order = 64 sweep has not finished")
const B = read_tsv(BF)

levels = sort(collect(intersect(keys(A), keys(B))))
isempty(levels) && error("no levels in common")

println("level  l_ec(Å)      dof    U(mo8)      U(mo64)     ΔU(eV)     Δ/U      t8(s)  t64(s)")
for lv in levels
    a, b = A[lv], B[lv]
    @printf("%5d  %7.4f  %9d  %.7f  %.7f  %+.2e  %+.2e  %6.0f %6.0f\n",
            lv, a[1], round(Int, a[3]), a[2], b[2], b[2] - a[2], (b[2] - a[2]) / a[2], a[4], b[4])
end

flat(d) = begin
    us = [d[lv][2] for lv in levels]
    (maximum(us) - minimum(us), maximum(abs, diff(us)))
end
s8, d8 = flat(A); s64, d64 = flat(B)
@printf("\nl_ec convergence over levels %d-%d\n", first(levels), last(levels))
@printf("  max_order = 8   spread %.2e eV   largest level-to-level step %.2e eV\n", s8, d8)
@printf("  max_order = 64  spread %.2e eV   largest level-to-level step %.2e eV\n", s64, d64)
println(s64 <= 10 * s8 ? "  -> convergence in l_ec survives the near-field fix" :
                         "  -> WARNING: the resolved sweep is NOT as flat; the old plateau may have been an artifact")

fin = last(levels)
@printf("\nbias at the finest level (%d): %.7f -> %.7f eV  (%+.2e eV, %+.3f%%)\n",
        fin, A[fin][2], B[fin][2], B[fin][2] - A[fin][2],
        100 * (B[fin][2] - A[fin][2]) / A[fin][2])
println("\nSec. 5.3 quotes the level-4 value; Sec. 5.2's cross-validation compares against it.")
