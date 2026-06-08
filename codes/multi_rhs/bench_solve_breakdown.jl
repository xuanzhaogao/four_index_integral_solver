# Break down the block-GMRES solve cost into: lfmm3d (FMM), near-correction matmul
# (corrections*X), dense glue (charges build + gradient contraction), and GMRES overhead
# (Krylov orthogonalization/restart = total - sum of matvec pieces). Same 2-shell graphene
# K=36 system as the per-RHS benchmark, with pinned OMP_NUM_THREADS.
#
#   julia --project scripts/bench_solve_breakdown.jl --build      # build + serialize once
#   OMP_NUM_THREADS=t julia --project scripts/bench_solve_breakdown.jl   # load + time

using BoundaryIntegral
using Krylov, LinearAlgebra, SparseArrays, FMM3D, Printf, Serialization
const BI = BoundaryIntegral
const OPS = joinpath(@__DIR__, "bench_solve_ops.jls")
const BIE = joinpath(@__DIR__, "graphene_2shell.bie")
const FMM_TOL = 1e-4   # FMM + near-correction upsampling tolerance (match bench_fmm.jl)

if "--build" in ARGS
    si = read_system_input(BIE); sp = si.solve
    group = BI.assemble_rhs_group(si, 1; support_rtol = sp.support_rtol)
    sources = group_volume_sources(group)
    interface = build_group_interface(si, group; n_quad = sp.n_quad, rhs_atol = sp.rhs_tol,
                                       l_ec = resolved_l_ec(si), eps_out = si.eps_out, max_depth = sp.max_depth)
    op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, FMM_TOL, FMM_TOL, sp.max_order)
    F  = rhs_dielectric_box3d_fmm3d(interface, sources, FMM_TOL)
    serialize(OPS, (op = op, F = F, rtol = FMM_TOL, npan = length(interface.panels)))
    @info "built+serialized  N=$(op.n)  panels=$(length(interface.panels))  K=$(size(F,2))"
    exit()
end

# instrumented operator: same matvec as BatchedDielectricOperator, timing FMM and matmul.
mutable struct TimedOp
    op::BatchedDielectricOperator
    t_fmm::Float64
    t_mm::Float64
    t_mul::Float64
    nmv::Int
end
Base.size(o::TimedOp) = size(o.op)
Base.size(o::TimedOp, d::Integer) = size(o.op, d)
Base.eltype(::TimedOp) = Float64

function LinearAlgebra.mul!(Y::AbstractMatrix, o::TimedOp, X::AbstractMatrix)
    tm0 = time_ns()
    op = o.op; n = op.n; K = size(X, 2)
    charges = Matrix{Float64}(undef, K, n)
    @inbounds for k in 1:K, i in 1:n
        charges[k, i] = op.weights[i] * X[i, k]
    end
    t = time_ns(); vals = lfmm3d(op.thresh, op.sources; charges = charges, pg = 2, nd = K); o.t_fmm += (time_ns() - t) / 1e9
    grad = reshape(vals.grad, K, 3, n)
    t = time_ns(); C = op.corrections * X; o.t_mm += (time_ns() - t) / 1e9
    @inbounds for k in 1:K, i in 1:n
        gn = op.norms[1, i] * grad[k, 1, i] + op.norms[2, i] * grad[k, 2, i] + op.norms[3, i] * grad[k, 3, i]
        Y[i, k] = -gn / (4π) + C[i, k] + op.diag[i] * X[i, k]
    end
    o.nmv += 1
    o.t_mul += (time_ns() - tm0) / 1e9
    return Y
end
LinearAlgebra.mul!(y::AbstractVector, o::TimedOp, x::AbstractVector) =
    (mul!(reshape(y, o.op.n, 1), o, reshape(x, length(x), 1)); y)
Base.:*(o::TimedOp, X::AbstractMatrix) = mul!(Matrix{Float64}(undef, o.op.n, size(X, 2)), o, X)
Base.:*(o::TimedOp, x::AbstractVector) = mul!(Vector{Float64}(undef, o.op.n), o, x)

d = deserialize(OPS); op = d.op; F = d.F; rtol = d.rtol; Kall = size(F, 2)
op * F[:, 1]   # warm up

omp = get(ENV, "OMP_NUM_THREADS", "default")
@printf("# OMP_NUM_THREADS=%s   N=%d   panels=%d   Kall=%d   rtol=%g\n", omp, op.n, d.npan, Kall, rtol)
@printf("# %3s %5s %5s %9s | %9s %10s %8s %9s | shares (fmm / matmul / glue / gmres)\n",
        "K", "iter", "nmv", "total(s)", "fmm(s)", "matmul(s)", "glue(s)", "gmres(s)")
for K in filter(k -> k <= Kall, [1, 36])
    GC.gc()
    to = TimedOp(op, 0.0, 0.0, 0.0, 0)
    t0 = time_ns()
    Σ, st = Krylov.block_gmres(to, F[:, 1:K]; rtol = rtol, itmax = 500)
    tt = (time_ns() - t0) / 1e9
    tglue = to.t_mul - to.t_fmm - to.t_mm
    tgm   = tt - to.t_mul
    @printf("  %3d %5d %5d %9.2f | %9.2f %10.3f %8.2f %9.2f | %4.0f%% / %4.1f%% / %3.0f%% / %4.0f%%\n",
            K, st.niter, to.nmv, tt, to.t_fmm, to.t_mm, tglue, tgm,
            100 * to.t_fmm / tt, 100 * to.t_mm / tt, 100 * tglue / tt, 100 * tgm / tt)
    flush(stdout)
end
