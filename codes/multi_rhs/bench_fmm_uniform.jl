# FMM3D throughput baseline on UNIFORM random points (3D), at the SAME particle counts as the
# slab-interface benchmark (bench_fmm.jl). Isolates whether the slab's thin/anisotropic surface
# point distribution hurts FMM performance vs an FMM-friendly 3D-uniform cloud.
#
# Same settings: fmm_tol=1e-4, pg=2, nd in {1,4,16,36}; uses env OMP_NUM_THREADS.
# Sweep:  for t in 1 2 4 8 16 32 64 96; do OMP_NUM_THREADS=$t julia --project scripts/bench_fmm_uniform.jl; done

using FMM3D, Printf, Random

const FMM_TOL = 1e-4
const NS  = [33696, 93312, 215136, 461376]   # match the slab interface sizes (levels 0..3)
const NDS = [1, 4, 16, 36]
const NREP = 3

Random.seed!(1)
best(f) = (f(); minimum(begin t = time_ns(); f(); (time_ns() - t) / 1e9 end for _ in 1:NREP))

omp = get(ENV, "OMP_NUM_THREADS", "unset")
@printf("# UNIFORM  OMP_NUM_THREADS=%s   fmm_tol=%g\n", omp, FMM_TOL)
@printf("# %4s %10s %5s %14s\n", "lvl", "N", "nd", "t_fmm(s)")
for (lvl, N) in enumerate(NS)
    sources = rand(3, N)                       # uniform in [0,1]^3
    for nd in NDS
        charges = rand(nd, N)
        tf = best(() -> lfmm3d(FMM_TOL, sources; charges = charges, pg = 2, nd = nd))
        @printf("  %4d %10d %5d %14.5f\n", lvl - 1, N, nd, tf)
        flush(stdout)
    end
end
