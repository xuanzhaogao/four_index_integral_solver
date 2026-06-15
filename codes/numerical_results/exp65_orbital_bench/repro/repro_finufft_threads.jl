# Standalone reproducer — FINUFFT 3D transform ~50x SLOWER at 96 threads than
# at 16 for certain mode-grid dimensions, while a nearby grid is unaffected.
# Depends ONLY on FINUFFT.jl 3.5.2 (finufft_jll 2.5.1, FFTW). No BIE/TKM3D/
# orbital data — synthetic uniform-random nonuniform points.
#
# Discovered via the graphene BIE pipeline (exp65): the screened-density
# volume-potential call ltkm3dc landed on mode grid (413,375,325) and took
# ~54 s at 96 threads; the same transform is ~1 s at 16 threads.
#
#   JULIA_NUM_THREADS=1 julia --project=exp65_orbital_bench/repro \
#       exp65_orbital_bench/repro/repro_finufft_threads.jl
# (Julia's own thread count is irrelevant; FINUFFT nthreads is set per call.)

import FINUFFT as FN
using Random, Printf

# min-of-ntrials wall time for a 3D type-1 nufft at the given mode grid +
# FINUFFT thread count (simple interface: plan + spread + FFT bundled).
function time_type1(x, y, z, c, ms, tol, nthreads; ntrials = 2)
    FN.nufft3d1(x, y, z, c, -1, tol, ms...; nthreads = nthreads)            # warm
    minimum(@elapsed(FN.nufft3d1(x, y, z, c, -1, tol, ms...; nthreads = nthreads)) for _ in 1:ntrials)
end

# type-2: time makeplan+setpts and exec separately (the FFT is in exec).
function time_type2(x, y, z, fk, ms, tol, nthreads; ntrials = 2)
    function once()
        tp = @elapsed begin
            plan = FN.finufft_makeplan(2, FN.BIGINT[ms...], 1, 1, tol; dtype = Float64, nthreads = nthreads)
            FN.finufft_setpts!(plan, x, y, z)
        end
        te = @elapsed FN.finufft_exec(plan, fk)
        FN.finufft_destroy!(plan)
        return tp, te
    end
    once()                                                                   # warm
    res = [once() for _ in 1:ntrials]
    return minimum(first.(res)), minimum(last.(res))
end

function main()
    N = 356_000
    tol = 1e-3
    Random.seed!(1)
    x = 2π .* rand(N) .- π; y = 2π .* rand(N) .- π; z = 2π .* rand(N) .- π
    c = complex.(randn(N))
    grids = [("BAD  ltkm3dc-box", (413, 375, 325)),
             ("GOOD field-box  ", (441, 401, 353))]
    threads = [1, 8, 16, 32, 64, 96]

    @printf("FINUFFT %s thread reproducer | N=%d  tol=%.0e  julia %s\n",
            "3.5.2 / finufft_jll 2.5.1", N, tol, string(VERSION))
    @printf("%-26s %5s  %10s  %10s  %10s\n", "grid (modes)", "nthr", "type1(s)", "t2_plan(s)", "t2_exec(s)")
    for (name, ms) in grids
        fk = FN.nufft3d1(x, y, z, c, -1, tol, ms...; nthreads = 16)          # coeff for the type-2 test
        for nth in threads
            t1 = time_type1(x, y, z, c, ms, tol, nth)
            tp, te = time_type2(x, y, z, fk, ms, tol, nth)
            @printf("%-26s %5d  %10.2f  %10.2f  %10.2f\n",
                    name * " " * string(ms), nth, t1, tp, te); flush(stdout)
            GC.gc()
        end
    end
    println("REPRO DONE")
end

main()
