# 6.5b — detailed benchmark of the RHS-adaptive mesh generation (the 201 s /
# 68% phase found by run_single_rhs.jl), and a prototype of the
# precompute-once evaluation strategy.
#
# Current implementation (single_dielectric_box3d_rhs_adaptive): at EVERY
# refinement depth, for the unsolved panels' ~136 test points each,
#   1. a KDTree over all source points is rebuilt to classify near/far,
#   2. far targets -> lfmm3d over ALL sources (tree + multipole setup per call),
#   3. near targets -> ltkm3dc over ALL sources (type-1 NUFFT + kernel scaling
#      + type-2 per call; the Fourier box depends on the targets, which is why
#      nothing is reusable).
# In the monolayer geometry the 5h-criterion classifies every target far, so
# the per-depth cost is dominated by the repeated FMM setup over the 356k-point
# density (plus the KDTree rebuild).
#
# Proposed strategy (prototype in part C): fix the near-evaluation region
# B = density bbox + margin*h ONCE (target-independent Fourier box), do the
# type-1 NUFFT + truncated-kernel diagonal scaling + spectral gradient
# coefficients ONCE, store them; per refinement round only a type-2 NUFFT at
# the new targets inside B, and a direct threaded sum for targets outside B
# (exact, no per-call setup). Validated against part A values; the refinement
# path is REPLAYED from part A's resolved bitmaps so both evaluators see
# identical target sets (decision agreement is reported separately).
#
# Parts (env PARTS, default "ABC"):
#   A  instrumented production build: per depth t(targets/classify/far FMM/
#      near TKM/interp check), counts; edge refinement + discretization time.
#   B  microbenchmarks: lfmm3d call floor vs n_targets; KDTree build+query;
#      direct-sum rate; TKM decomposition (type-1 / scaling / grad coeffs /
#      type-2) at the production kmax and at the spectral cutoff
#      estimate_kcut3dc(tol); kcut value and mode-grid sizes.
#   C  strategy prototype, end to end against part A.
#
# Run:  JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#         julia --project=. exp65_orbital_bench/scripts/bench_meshgen.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Printf
using Serialization

const TKM = BI.TKM3D
const FN = TKM.FINUFFT

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const PARTS = uppercase(get(ENV, "PARTS", "ABC"))
const MARGIN_H = parse(Float64, get(ENV, "MARGIN_H", "5.0"))
const DATA = joinpath(@__DIR__, "..", "data")
mkpath(joinpath(DATA, "raw"))

# production parameters (identical to run_single_rhs.jl)
const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
const N_QUAD, EDGE_LEVEL, RHS_TOL = 6, 4, 1e-3
const FMM_TOL = RHS_TOL * 0.1
const L_EC = LZ / 2.0^EDGE_LEVEL * 1.01

# ---------------------------------------------------------------------------
# setup: orbital -> screened source (untimed here; benchmarked in run_single_rhs)
# ---------------------------------------------------------------------------
println(">>> setup: load orbital + screened source"); flush(stdout)
const SVS = let
    dg = load_squared_xsf(XSF_1)
    sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    vs1 = BI.VolumeSource(dg, tol = 1e-3)
    BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
end
const H = BI._estimate_source_spacing(SVS)
const KMAX_PROD = BI._estimate_tkm3dc_kmax(H)
const SRC, Q = BI._volume_source_fmm_sources(SVS)
const NSRC = length(Q)
const SRC_LO = vec(minimum(SRC; dims = 2))
const SRC_HI = vec(maximum(SRC; dims = 2))
@printf("  src points %d   h = %.4f   kmax(prod, Nyquist) = %.3f\n", NSRC, H, KMAX_PROD)
@printf("  density bbox  x [%.2f, %.2f]  y [%.2f, %.2f]  z [%.2f, %.2f]\n",
        SRC_LO[1], SRC_HI[1], SRC_LO[2], SRC_HI[2], SRC_LO[3], SRC_HI[3])
flush(stdout)

const NS, WS = BI.gausslegendre(N_QUAD)
const LAM = BI.gl_barycentric_weights(NS, WS)
const XS = range(-1.0, 1.0; length = 10)

# resolution check, replicated verbatim from _rhs_panel3d_resolved_volume_fmm
function resolve_flags(rhs_vals, n_panels)
    resolved = fill(false, n_panels)
    idx = 0
    for p in 1:n_panels
        quad_vals = Matrix{Float64}(undef, N_QUAD, N_QUAD)
        for i in 1:N_QUAD, j in 1:N_QUAD
            idx += 1
            quad_vals[i, j] = rhs_vals[idx]
        end
        err = 0.0
        for u in XS
            rx = BI.barycentric_row(NS, LAM, u)
            for v in XS
                ry = BI.barycentric_row(NS, LAM, v)
                approx = 0.0
                for i in 1:N_QUAD, j in 1:N_QUAD
                    approx += quad_vals[i, j] * rx[i] * ry[j]
                end
                idx += 1
                err = max(err, abs(rhs_vals[idx] - approx))
            end
        end
        resolved[p] = err <= RHS_TOL
    end
    return resolved
end

function subdivide(unsolved, resolved)
    solved_add = BI.TempPanel3D{Float64}[]
    nxt = BI.TempPanel3D{Float64}[]
    for i in eachindex(unsolved)
        resolved[i] ? push!(solved_add, unsolved[i]) :
                      append!(nxt, BI.divide_temp_panel3d(unsolved[i], 2, 2))
    end
    return solved_add, nxt
end

function direct_grad_rhs!(out, srcs, q, targets, normals, idxs)
    sx = vec(srcs[1, :]); sy = vec(srcs[2, :]); sz = vec(srcs[3, :])
    Threads.@threads for ii in eachindex(idxs)
        i = idxs[ii]
        x, y, z = targets[1, i], targets[2, i], targets[3, i]
        gx = 0.0; gy = 0.0; gz = 0.0
        @inbounds @simd for j in eachindex(q)
            dx = x - sx[j]; dy = y - sy[j]; dz = z - sz[j]
            r2 = dx * dx + dy * dy + dz * dz
            s = q[j] / (r2 * sqrt(r2))
            gx -= s * dx; gy -= s * dy; gz -= s * dz
        end
        out[i] = -(normals[1, i] * gx + normals[2, i] * gy + normals[3, i] * gz) / (4π)
    end
end

# fixed-box spectral data for the proposed strategy / part B decomposition
function spectral_box(km)
    lo = SRC_LO .- MARGIN_H * H
    hi = SRC_HI .+ MARGIN_H * H
    corners = hcat(lo, hi)
    lengths, center = TKM.combined_box_geometry_3xn(SRC, corners)
    Lbig = sqrt(sum(abs2, lengths))
    dks = [prevfloat(2π / (lengths[d] + Lbig)) for d in 1:3]
    kx = TKM.centered_mode_axis(dks[1], km)
    ky = TKM.centered_mode_axis(dks[2], km)
    kz = TKM.centered_mode_axis(dks[3], km)
    return (; lo, hi, center, Lbig, dks, kx, ky, kz,
            nm = (length(kx), length(ky), length(kz)))
end

function scaled_type1!(box)
    sxn = box.dks[1] .* (vec(SRC[1, :]) .- box.center[1])
    syn = box.dks[2] .* (vec(SRC[2, :]) .- box.center[2])
    szn = box.dks[3] .* (vec(SRC[3, :]) .- box.center[3])
    t1 = @elapsed coeff0 = FN.nufft3d1(sxn, syn, szn, complex.(Q), -1, FMM_TOL, box.nm...)
    coeff = ndims(coeff0) == 4 ? dropdims(coeff0; dims = 4) : coeff0
    km = maximum(abs, box.kx)  # cutoff actually used to build the axes
    t2 = @elapsed @inbounds for iz in eachindex(box.kz), iy in eachindex(box.ky), ix in eachindex(box.kx)
        k = sqrt(box.kx[ix]^2 + box.ky[iy]^2 + box.kz[iz]^2)
        coeff[ix, iy, iz] = k <= km ? coeff[ix, iy, iz] * TKM.truncated_laplace3d_hat(k, box.Lbig) :
                                      zero(eltype(coeff))
    end
    t3 = @elapsed gradc = TKM._spectral_gradient_coeffs_3d(coeff, box.kx, box.ky, box.kz)
    return coeff, gradc, t1, t2, t3
end

# ---------------------------------------------------------------------------
# Part A — instrumented production build
# ---------------------------------------------------------------------------
function part_A()
    println(">>> PART A: instrumented adaptive build (production path)"); flush(stdout)
    BI.lfmm3d(FMM_TOL, SRC[:, 1:5000]; charges = Q[1:5000],
              targets = rand(3, 64) .* 10.0, pgt = 2)   # warm

    depth_rows = NamedTuple[]
    rhs_store = Vector{Float64}[]
    resolved_store = Vector{Bool}[]

    t0 = time()
    unsolved = BI._box3d_rhs_adaptive_initial_panels(L, L, LZ, sqrt(2.0))
    solved = BI.TempPanel3D{Float64}[]
    depth = 0
    while !isempty(unsolved) && depth < 12
        t_t = @elapsed ((targets, normals, _) =
            BI._rhs_panel3d_refinement_targets(unsolved, NS, WS; n_pts = 10))
        t_c = @elapsed is_near = BI._classify_near_far_targets(targets, SVS, H)
        n_near = count(is_near); n_far = length(is_near) - n_near

        rhs_vals = Vector{Float64}(undef, size(targets, 2))
        t_far = 0.0; t_near = 0.0
        if n_far > 0
            fidx = findall(!, is_near)
            t_far = @elapsed begin
                vals = BI.lfmm3d(FMM_TOL, SRC; charges = Q,
                                 targets = targets[:, fidx], pgt = 2)
                for (k, i) in enumerate(fidx)
                    rhs_vals[i] = -dot(view(normals, :, i), view(vals.gradtarg, :, k)) / (4π)
                end
            end
        end
        if n_near > 0
            nidx = findall(is_near)
            t_near = @elapsed begin
                vals = TKM.ltkm3dc(FMM_TOL, SRC; charges = Q,
                                   targets = targets[:, nidx], pgt = 2, kmax = KMAX_PROD)
                vals.ier == 0 || error("ltkm3dc failed")
                for (k, i) in enumerate(nidx)
                    rhs_vals[i] = -dot(view(normals, :, i), view(vals.gradtarg, :, k))
                end
            end
        end
        t_chk = @elapsed resolved = resolve_flags(rhs_vals, length(unsolved))

        push!(depth_rows, (; depth, n_panels = length(unsolved),
              n_targets = size(targets, 2), n_near, n_far,
              t_targets = t_t, t_classify = t_c, t_far, t_near, t_check = t_chk))
        push!(rhs_store, rhs_vals); push!(resolved_store, resolved)
        @printf("  depth %d: panels %5d  targets %6d (near %d / far %d)  classify %5.2fs  FMM %6.2fs  TKM %6.2fs  check %5.2fs\n",
                depth, length(unsolved), size(targets, 2), n_near, n_far, t_c, t_far, t_near, t_chk)
        flush(stdout)

        add, unsolved = subdivide(unsolved, resolved)
        append!(solved, add)
        depth += 1
    end
    append!(solved, unsolved)

    t_e0 = time()
    rough = copy(solved); refined = BI.TempPanel3D{Float64}[]
    while !isempty(rough)
        tpl = popfirst!(rough)
        has_ec = tpl.is_a_corner || tpl.is_b_corner || tpl.is_c_corner || tpl.is_d_corner ||
                 tpl.is_ab_edge || tpl.is_bc_edge || tpl.is_cd_edge || tpl.is_da_edge
        if has_ec && max(norm(tpl.b .- tpl.a), norm(tpl.a .- tpl.d)) > L_EC
            append!(rough, BI.divide_temp_panel3d(tpl, 2, 2))
        else
            push!(refined, tpl)
        end
    end
    total = time() - t0
    @printf("  edge refinement: %.2fs -> %d panels (%d points; production run had 960768)\n",
            time() - t_e0, length(refined), length(refined) * N_QUAD^2)
    @printf("  PART A total: %.1f s\n", total); flush(stdout)
    return depth_rows, rhs_store, resolved_store, total
end

# ---------------------------------------------------------------------------
# Part B — component microbenchmarks
# ---------------------------------------------------------------------------
function part_B()
    println(">>> PART B: component microbenchmarks"); flush(stdout)
    trg_probe(n) = SRC_LO .+ rand(3, n) .* (SRC_HI .- SRC_LO) .+ [0.0, 0.0, 3.0]

    BI.lfmm3d(FMM_TOL, SRC; charges = Q, targets = trg_probe(64), pgt = 2)  # warm
    for n in (136, 1_000, 10_000, 100_000)
        t = @elapsed BI.lfmm3d(FMM_TOL, SRC; charges = Q, targets = trg_probe(n), pgt = 2)
        @printf("  [B1] lfmm3d  %6d targets: %6.2f s\n", n, t); flush(stdout)
    end

    t_tree = @elapsed BI.NearestNeighbors.KDTree(SVS.positions)
    t = @elapsed BI._classify_near_far_targets(trg_probe(10_000), SVS, H)
    @printf("  [B2] KDTree build over %d pts: %.2f s   classify 10k targets (incl. rebuild): %.2f s\n",
            NSRC, t_tree, t); flush(stdout)

    let n = 10_000, trg = trg_probe(n), nrm = vcat(zeros(2, n), ones(1, n)), out = zeros(n)
        direct_grad_rhs!(out, SRC, Q, trg, nrm, collect(1:100))  # warm
        t = @elapsed direct_grad_rhs!(out, SRC, Q, trg, nrm, collect(1:n))
        @printf("  [B3] direct sum  %d targets x %d sources: %.2f s  (%.2e pair/s)\n",
                n, NSRC, t, n * NSRC / t); flush(stdout)
    end

    t_kcut = @elapsed kc = TKM.estimate_kcut3dc(SRC; charges = Q, tol = FMM_TOL)
    @printf("  [B4] estimate_kcut3dc(tol=%.0e): kcut = %.2f  (Nyquist %.2f, source-box modes %s)  [%.1f s]\n",
            FMM_TOL, kc.kcut, kc.kmax_nyquist, kc.nmodes, t_kcut); flush(stdout)

    for (lbl, km) in (("kmax=prod", KMAX_PROD), ("kmax=kcut", Float64(kc.kcut)))
        box = spectral_box(km)
        @printf("  [B4] %s (%.2f): modes %s = %.2e  coeff %.2f GB (grad x3)\n",
                lbl, km, box.nm, prod(box.nm), prod(box.nm) * 16 / 2^30); flush(stdout)
        if prod(box.nm) > 4e8
            println("       too many modes, skipping decomposition"); continue
        end
        coeff, gradc, t1, t2, t3 = scaled_type1!(box)
        for n in (1_000, 10_000)
            trg = trg_probe(n)
            txn = box.dks[1] .* (vec(trg[1, :]) .- box.center[1])
            tyn = box.dks[2] .* (vec(trg[2, :]) .- box.center[2])
            tzn = box.dks[3] .* (vec(trg[3, :]) .- box.center[3])
            t4 = @elapsed begin
                plan = TKM._finufft_make_type2_plan_3d(txn, tyn, tzn, 1, FMM_TOL, box.nm, 3, Float64)
                FN.finufft_exec(plan, gradc)
                FN.finufft_destroy!(plan)
            end
            @printf("       type-1 %.2fs  scale %.2fs  gradcoeff %.2fs  type-2(%d trg, ntrans=3) %.2fs\n",
                    t1, t2, t3, n, t4); flush(stdout)
        end
        coeff = nothing; gradc = nothing; GC.gc()
    end
end

# ---------------------------------------------------------------------------
# Part C — strategy prototype: precompute once, type-2 + direct per depth
# ---------------------------------------------------------------------------
function part_C(rhs_store, resolved_store, A_total)
    println(">>> PART C: precompute-once strategy (replaying part A refinement path)"); flush(stdout)
    km = haskey(ENV, "TKM_KMAX") ? parse(Float64, ENV["TKM_KMAX"]) :
         min(Float64(TKM.estimate_kcut3dc(SRC; charges = Q, tol = FMM_TOL).kcut), KMAX_PROD)

    local box, gradc, pref
    t_pre = @elapsed begin
        box = spectral_box(km)
        prod(box.nm) > 4e8 && error("mode grid too large: $(box.nm)")
        _, gradc, _, _, _ = scaled_type1!(box)
        pref = prod(box.dks) / (2π)^3
    end
    @printf("  precompute (type-1 + scale + grad coeffs): %.2f s   kmax=%.2f  modes %s\n",
            t_pre, km, box.nm); flush(stdout)

    in_B(targets, i) = (box.lo[1] <= targets[1, i] <= box.hi[1]) &&
                       (box.lo[2] <= targets[2, i] <= box.hi[2]) &&
                       (box.lo[3] <= targets[3, i] <= box.hi[3])

    total_eval = 0.0
    maxdiff = 0.0; maxref = 0.0; n_flips = 0
    unsolved = BI._box3d_rhs_adaptive_initial_panels(L, L, LZ, sqrt(2.0))
    for (d, resolvedA) in enumerate(resolved_store)
        targets, normals, _ = BI._rhs_panel3d_refinement_targets(unsolved, NS, WS; n_pts = 10)
        nt = size(targets, 2)
        rhs_vals = Vector{Float64}(undef, nt)

        t_split = @elapsed begin
            inb = [in_B(targets, i) for i in 1:nt]
            bidx = findall(inb); oidx = findall(!, inb)
        end
        t_t2 = 0.0
        if !isempty(bidx)
            t_t2 = @elapsed begin
                txn = box.dks[1] .* (vec(targets[1, bidx]) .- box.center[1])
                tyn = box.dks[2] .* (vec(targets[2, bidx]) .- box.center[2])
                tzn = box.dks[3] .* (vec(targets[3, bidx]) .- box.center[3])
                plan = TKM._finufft_make_type2_plan_3d(txn, tyn, tzn, 1, FMM_TOL, box.nm, 3, Float64)
                gc_ = FN.finufft_exec(plan, gradc)
                FN.finufft_destroy!(plan)
                g = pref .* real.(gc_)               # n_in x 3
                for (k, i) in enumerate(bidx)
                    rhs_vals[i] = -(normals[1, i] * g[k, 1] + normals[2, i] * g[k, 2] +
                                    normals[3, i] * g[k, 3])
                end
            end
        end
        t_dir = @elapsed direct_grad_rhs!(rhs_vals, SRC, Q, targets, normals, oidx)
        total_eval += t_split + t_t2 + t_dir

        dref = rhs_store[d]
        for i in 1:nt
            maxdiff = max(maxdiff, abs(rhs_vals[i] - dref[i]))
            maxref = max(maxref, abs(dref[i]))
        end
        resolvedC = resolve_flags(rhs_vals, length(unsolved))
        n_flips += count(resolvedC .!= resolvedA)
        @printf("  depth %d: targets %6d (inB %d / out %d)  type-2 %5.2fs  direct %6.2fs\n",
                d - 1, nt, length(bidx), length(oidx), t_t2, t_dir)
        flush(stdout)

        _, unsolved = subdivide(unsolved, resolvedA)   # replay part A path
    end
    @printf("  PART C: precompute %.1f s + per-depth eval %.1f s = %.1f s   (part A eval loop: %.1f s)\n",
            t_pre, total_eval, t_pre + total_eval, A_total)
    @printf("  validation: max |rhs_C - rhs_A| = %.3e  (max |rhs_A| = %.3e, rel %.2e)   decision flips: %d\n",
            maxdiff, maxref, maxdiff / maxref, n_flips)
    return (; t_pre, total_eval, maxdiff, maxref, n_flips, kmax_used = km)
end

# ---------------------------------------------------------------------------
function main()
    depth_rows = NamedTuple[]; rhs_store = Vector{Float64}[]
    resolved_store = Vector{Bool}[]; A_total = NaN
    if occursin("A", PARTS)
        depth_rows, rhs_store, resolved_store, A_total = part_A()
    end
    occursin("B", PARTS) && part_B()
    C = (occursin("C", PARTS) && !isempty(resolved_store)) ?
        part_C(rhs_store, resolved_store, A_total) : nothing
    serialize(joinpath(DATA, "raw", "bench_meshgen.jls"),
              (; depth_rows, A_total, C, margin_h = MARGIN_H,
                 nthreads = Threads.nthreads()))
    println("MESHGEN BENCH DONE")
end

main()
