# Verify: is (413,375,325) the OUTPUT mode grid or the internal upsampled fine grid?
# FINUFFT debug=1 prints the fine grid nf1/nf2/nf3 that the FFT actually runs on.
# Also contrasts upsampfac 2.0 vs 1.25 (cost + the nf it produces).
#   julia --project=. exp65_orbital_bench/scripts/diag_upsampfac.jl
import FINUFFT
using Printf, Random
Random.seed!(1)

const N = (413, 375, 325)            # the grid reported in issue #869
const NSRC = 350_000
x = π .* (2 .* rand(NSRC) .- 1)
y = π .* (2 .* rand(NSRC) .- 1)
z = π .* (2 .* rand(NSRC) .- 1)
c = complex.(rand(NSRC))

@printf("requested OUTPUT modes N = %s  (= %d)\n", N, prod(N))
println("physical cores on this host: ", Sys.CPU_THREADS)

for ups in (2.0, 1.25)
    println("\n===================  upsampfac = $ups  ===================")
    # debug=1 -> FINUFFT prints "nf1=.. nf2=.. nf3=.." (the upsampled fine grid it FFTs)
    plan = FINUFFT.finufft_makeplan(1, collect(N), -1, 1, 1e-3; upsampfac = ups, nthreads = 8, debug = 1)
    FINUFFT.finufft_setpts!(plan, x, y, z)
    out = FINUFFT.finufft_exec(plan, c)
    @printf("  -> OUTPUT array size = %s   (FFT runs on the nf* grid printed above)\n", size(out))
    FINUFFT.finufft_destroy!(plan)

    for nthr in (1, 8, 16, 32)
        nthr > Sys.CPU_THREADS && continue
        p = FINUFFT.finufft_makeplan(1, collect(N), -1, 1, 1e-3; upsampfac = ups, nthreads = nthr)
        FINUFFT.finufft_setpts!(p, x, y, z)
        FINUFFT.finufft_exec(p, c)                 # warm
        t = @elapsed FINUFFT.finufft_exec(p, c)
        @printf("  nthreads=%2d  exec=%.3f s\n", nthr, t)
        FINUFFT.finufft_destroy!(p)
    end
end
println("DONE")
