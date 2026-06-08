# Micro-benchmark of the dielectric matvec's two pieces — the FMM3D kernel and the near
# correction — as a function of OMP_NUM_THREADS, number of RHS (nd), and number of unknowns N.
#
# Monolayer slab (L=90, Lz=3.35, eps_in=3.5) with a FAR point source (smooth RHS -> the
# interface refinement is geometry-driven). Edge-refinement level 0..3 -> different N. FMM
# tolerance fixed at 1e-4. NO solve: one matvec's worth of work, timing the FMM
# (lfmm3d nd=K, pg=2) and the near correction (corrections * X) separately (best of NREP).
#
# The operators are built once and serialized so that the per-thread runs only load + time
# (the build is not part of the measurement). Usage, sweeping threads:
#   for t in 1 2 4 8 16 32 64 96; do OMP_NUM_THREADS=$t julia --project scripts/bench_fmm.jl; done

using BoundaryIntegral
using FMM3D, SparseArrays, LinearAlgebra, Printf, Serialization
const BI = BoundaryIntegral

const FMM_TOL  = 1e-4
const LEVELS   = [0, 1, 2, 3]          # edge_refine_level -> N
const NDS      = [1, 4, 16, 36]        # number of right-hand sides (FMM nd)
const NREP     = 3
const OPS_FILE = joinpath(@__DIR__, "bench_fmm_ops.jls")

function build_ops()
    boxes  = [(center = (0.0, 0.0, 7.5), Lx = 90.0, Ly = 90.0, Lz = 3.35)]
    epses  = [3.5]; eps_out = 1.0
    ps     = PointSource((0.0, 0.0, 50.0), 1.0)   # far source -> geometry-driven refinement only
    n_quad = 6; max_order = 8
    ops = NamedTuple[]
    for level in LEVELS
        l_ec = 3.35 / 2.0^level
        interface = multi_dielectric_box3d_rhs_adaptive(
            n_quad, l_ec, boxes, epses, ps, 1.0, 1e-3, eps_out; max_depth = 20)
        op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, FMM_TOL, FMM_TOL, max_order)
        push!(ops, (level = level, n = op.n, sources = op.sources, weights = op.weights,
                    norms = op.norms, thresh = op.thresh,
                    corrections = SparseMatrixCSC{Float64,Int}(op.corrections)))
        @info "built level=$level  N=$(op.n)  nnz(corr)=$(nnz(op.corrections))"
    end
    return ops
end

if "--build" in ARGS
    serialize(OPS_FILE, build_ops())
    @info "built + serialized operators -> $OPS_FILE"
    exit()
end
# timing mode: the serialized operators must already exist (build once with --build first,
# so concurrent job-array tasks don't race to build/serialize).
isfile(OPS_FILE) || error("$OPS_FILE not found — run `julia --project scripts/bench_fmm.jl --build` first")
ops = deserialize(OPS_FILE)

best(f) = (f(); minimum(begin t = time_ns(); f(); (time_ns() - t) / 1e9 end for _ in 1:NREP))

omp = get(ENV, "OMP_NUM_THREADS", "unset")
@printf("# OMP_NUM_THREADS=%s   fmm_tol=%g   nthreads(julia)=%d\n", omp, FMM_TOL, Threads.nthreads())
@printf("# %4s %10s %5s %14s %14s\n", "lvl", "N", "nd", "t_fmm(s)", "t_near(s)")
for o in ops
    n = o.n
    for nd in NDS
        X = randn(n, nd)
        charges = Matrix{Float64}(undef, nd, n)
        @inbounds for k in 1:nd, i in 1:n
            charges[k, i] = o.weights[i] * X[i, k]
        end
        Y = similar(X)
        t_fmm  = best(() -> lfmm3d(o.thresh, o.sources; charges = charges, pg = 2, nd = nd))
        t_near = best(() -> mul!(Y, o.corrections, X))
        @printf("  %4d %10d %5d %14.5f %14.5f\n", o.level, n, nd, t_fmm, t_near)
        flush(stdout)
    end
end
