# One-off: extend the FMM kernel sweep to nd=64 at N=461,376 (slab, edge-refine level 3), so the
# nd-amortization in plot_speedup_grid.jl can be checked for convergence. Same setup as bench_fmm.jl
# (level 3). Build the operator once, then time nd=64 at each pinned OMP_NUM_THREADS:
#   julia --project bench_nd64.jl --build
#   for t in 1 2 4 8 16 32 64 96; do OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t julia --project bench_nd64.jl; done
# Emits "ROW,omp,level,N,nd,t_fmm,t_near" lines to append to data/fmm_bench.csv.

using BoundaryIntegral
using FMM3D, SparseArrays, LinearAlgebra, Printf, Serialization
const BI = BoundaryIntegral
const FMM_TOL = 1e-4
const ND  = 64
const NREP = 1                       # nd=64 at N=461k is slow; one timed call (+warmup)
const OPS = joinpath(@__DIR__, "bench_nd64_ops.jls")

if "--build" in ARGS
    boxes = [(center = (0.0, 0.0, 7.5), Lx = 90.0, Ly = 90.0, Lz = 3.35)]
    epses = [3.5]; eps_out = 1.0
    ps = PointSource((0.0, 0.0, 50.0), 1.0)
    n_quad = 6; max_order = 8; l_ec = 3.35 / 2.0^3      # level 3
    interface = multi_dielectric_box3d_rhs_adaptive(
        n_quad, l_ec, boxes, epses, ps, 1.0, 1e-3, eps_out; max_depth = 20)
    op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, FMM_TOL, FMM_TOL, max_order)
    serialize(OPS, (n = op.n, sources = op.sources, weights = op.weights, norms = op.norms,
                    thresh = op.thresh, corrections = SparseMatrixCSC{Float64,Int}(op.corrections)))
    @info "built nd64 ops  N=$(op.n)"
    exit()
end

o = deserialize(OPS); n = o.n
best(f) = (f(); minimum(begin t = time_ns(); f(); (time_ns() - t) / 1e9 end for _ in 1:NREP))
X = randn(n, ND)
charges = Matrix{Float64}(undef, ND, n)
@inbounds for k in 1:ND, i in 1:n
    charges[k, i] = o.weights[i] * X[i, k]
end
Y = similar(X)
t_fmm  = best(() -> lfmm3d(o.thresh, o.sources; charges = charges, pg = 2, nd = ND))
t_near = best(() -> mul!(Y, o.corrections, X))
omp = get(ENV, "OMP_NUM_THREADS", "unset")
@printf("ROW,%s,3,%d,%d,%.5f,%.6f\n", omp, n, ND, t_fmm, t_near)
