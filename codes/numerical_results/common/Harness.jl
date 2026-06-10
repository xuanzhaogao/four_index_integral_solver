# Harness.jl — shared instrumented pipeline for the Section-6 experiments.
#
# Mirrors the production pipeline (ScreenedOrbitalSolve.solve_screened_mode /
# timing_density_solve.jl) with the body inlined so that every phase is timed
# and the near-correction structures (pair counts, p_up, sparse memory) are
# introspectable. See ../PLAN.md for the parameter mapping.

module Harness

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using LinearMaps
using SparseArrays
using Random
using Statistics
using Printf
using Serialization
using Dates

export SystemSpec, system1, system2, slab_system,
       gaussian_source, source_grid_n, screened_source,
       slab_internal, solve_system, eval_phi, eval_V, vacuum_V,
       zone_targets, eval_scatter_with_h0,
       uniform_refine, refine_to_dof,
       bernstein_rho_min, run_provenance, append_csv_row, run_cols, times_cols,
       save_ref, load_ref, FMM_FLOOR

const FMM_FLOOR = 1e-14          # FMM3D saturates at ~3e-14 rel. error (verified in smoke test)
const SEED_TARGETS = 20260609

fmm_tol_of(eps::Float64) = max(eps, FMM_FLOOR)

# ---------------------------------------------------------------------------
# Systems
# ---------------------------------------------------------------------------

struct SystemSpec
    name::String
    boxes::Vector{NamedTuple{(:center, :Lx, :Ly, :Lz), Tuple{NTuple{3, Float64}, Float64, Float64, Float64}}}
    epses::Vector{Float64}
    eps_out::Float64
    src_center::NTuple{3, Float64}
    src_sigma::Float64
    tgt_center::NTuple{3, Float64}
    tgt_sigma::Float64
end

box(center, Lx, Ly, Lz) = (center = center, Lx = Float64(Lx), Ly = Float64(Ly), Lz = Float64(Lz))

"System I: unit cube eps1 in vacuum; source Gaussian s=0.05 at standoff d above top face."
function system1(; eps1::Float64 = 10.0, d::Float64 = 0.1, s::Float64 = 0.05)
    src = (0.0, 0.0, 0.5 + d)
    tgt = (0.2, 0.0, 0.5 + d)
    SystemSpec("system1_eps$(eps1)_d$(d)",
        [box((0.0, 0.0, 0.0), 1.0, 1.0, 1.0)], [eps1], 1.0, src, s, tgt, s)
end

"System II: two coplanar substrate boxes + material box over the junction (paper Fig. 1)."
function system2(; eps1::Float64 = 4.0, eps2::Float64 = 12.0, epsm::Float64 = 2.0,
                  g::Float64 = 0.1, s::Float64 = 0.05)
    zm = 0.5 + g + 0.1          # material box center height (box Lz=0.2)
    SystemSpec("system2",
        [box((-0.5, 0.0, 0.25), 1.0, 1.0, 0.5),
         box((0.5, 0.0, 0.25), 1.0, 1.0, 0.5),
         box((0.0, 0.0, zm), 0.6, 0.6, 0.2)],
        [eps1, eps2, epsm], 1.0,
        (0.0, 0.0, zm), s, (0.2, 0.0, zm), s)
end

"6.1 slab: L x L x Lz slab (eps1) with the Gaussian source fully supported inside.
Source at the slab center; V-target displaced laterally by 0.2 at the same height."
function slab_internal(; L::Float64 = 10.0, Lz::Float64 = 1.0, eps1::Float64 = 10.0,
                        s::Float64 = 0.05)
    SystemSpec("slab_internal_L$(L)_Lz$(Lz)_eps$(eps1)",
        [box((0.0, 0.0, 0.0), L, L, Lz)], [eps1], 1.0,
        (0.0, 0.0, 0.0), s, (0.2, 0.0, 0.0), s)
end

"6.2/6.4 slab: A x A x 0.5 slab (eps 10) + material box (eps 2) at gap g above it."
function slab_system(; A::Float64 = 10.0, eps_slab::Float64 = 10.0, epsm::Float64 = 2.0,
                      g::Float64 = 0.05, s::Float64 = 0.05)
    zm = 0.25 + g + 0.1
    SystemSpec("slab_A$(A)",
        [box((0.0, 0.0, 0.0), A, A, 0.5),
         box((0.0, 0.0, zm), 0.6, 0.6, 0.2)],
        [eps_slab, epsm], 1.0,
        (0.0, 0.0, zm), s, (0.2, 0.0, zm), s)
end

min_box_dim(sys::SystemSpec) = minimum(min(b.Lx, b.Ly, b.Lz) for b in sys.boxes)
l_ec_of(sys::SystemSpec, r::Int) = min_box_dim(sys) / 2.0^r * 1.01

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

"""
Grid points per dimension so that (a) truncation tol = eps/10 and (b) the
midpoint-rule aliasing error exp(-(pi*sigma/h)^2/2) <= eps. margin: safety factor.
"""
function source_grid_n(eps::Float64; margin::Float64 = 1.25)
    L1 = log(10.0 / eps)   # truncation
    L2 = log(1.0 / eps)    # aliasing
    n = ceil(Int, margin * (4.0 / pi) * sqrt(L1 * L2))
    return max(n, 8)
end

"Normalized isotropic Gaussian volume source on an n^3 grid truncated at eps/10."
function gaussian_source(center::NTuple{3, Float64}, sigma::Float64, eps::Float64;
                         n::Union{Nothing, Int} = nothing, margin::Float64 = 1.25)
    n_res = isnothing(n) ? source_grid_n(eps; margin = margin) : n
    return BI.GaussianVolumeSource(center, sigma, n_res, eps / 10.0)
end

"Per-point dielectric screening rho -> rho/eps(x) (multi-box aware)."
function screened_source(sys::SystemSpec, vs::BI.VolumeSource{Float64, 3})
    return BI.screened_volume_source(sys.boxes, sys.epses, sys.eps_out, vs, BI.SharpScreening())
end

# ---------------------------------------------------------------------------
# Instrumented solve  (inlined production pipeline)
# ---------------------------------------------------------------------------

function _diag_coeffs(interface)
    n = BI.num_points(interface)
    t = Vector{Float64}(undef, n)
    offset = 0
    for i in 1:length(interface.panels)
        npts = BI.num_points(interface.panels[i])
        ti = 0.5 * (interface.eps_out[i] + interface.eps_in[i]) /
                   (interface.eps_out[i] - interface.eps_in[i])
        t[(offset + 1):(offset + npts)] .= ti
        offset += npts
    end
    return t
end

"""
    solve_system(sys; eps, p, r, ...) -> NamedTuple

Instrumented end-to-end solve. All component tolerances tied to `eps`
(FMM clamped at FMM_FLOOR). Returns interface, sigma, gmres stats, neighbor
data, timings, and the screened source for downstream evaluation.

`interface_override` (plus `screened_vs_override`/`tkm_kmax_override`) lets
variant-C / scaling experiments supply a custom interface while reusing the
identical assembly/solve path. `near_correction = false` gives variant B.
"""
function solve_system(sys::SystemSpec;
    eps::Float64, p::Int, r::Int,
    max_order::Int = 64, max_depth::Int = 128, itmax::Int = 500,
    src_n::Union{Nothing, Int} = nothing, src_margin::Float64 = 1.25,
    near_correction::Bool = true,
    interface_override = nothing,
    screened_vs_override = nothing,
    tkm_kmax_override = nothing,
    history::Bool = true,
    gmres_atol::Float64 = 1e-14,
    gmres_verbose::Int = 0,
)
    fmm_tol = fmm_tol_of(eps)
    tkm_tol = max(eps, 1e-13)
    times = Dict{String, Float64}()

    # 1. source prep (screened density + TKM kmax)
    t0 = time()
    if screened_vs_override === nothing
        vs = gaussian_source(sys.src_center, sys.src_sigma, eps; n = src_n, margin = src_margin)
        screened_vs = screened_source(sys, vs)
    else
        screened_vs = screened_vs_override
    end
    tkm_kmax = tkm_kmax_override === nothing ?
        BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(screened_vs)) : Float64(tkm_kmax_override)
    times["source_prep"] = time() - t0

    # 2. interface build (RHS-adaptive + edge refinement)
    t0 = time()
    if interface_override === nothing
        interface = BI.multi_dielectric_box3d_rhs_adaptive(
            p, l_ec_of(sys, r), sys.boxes, sys.epses, screened_vs, 1.0, eps, sys.eps_out;
            max_depth = max_depth, tkm_kmax = tkm_kmax)
    else
        interface = interface_override
    end
    times["interface_build"] = time() - t0
    N = BI.num_points(interface)

    # 3. RHS assembly (hybrid: TKM near + FMM far)
    t0 = time()
    rhs = BI.rhs_dielectric_box3d_hybrid(interface, screened_vs, 1.0, fmm_tol; tkm_kmax = tkm_kmax)
    times["rhs_assembly"] = time() - t0

    # 4. LHS: FMM base + sparse near correction (inlined for introspection)
    t0 = time()
    D_base = BI.laplace3d_DT_fmm3d(interface, fmm_tol)
    times["fmm_setup"] = time() - t0

    t0 = time()
    local nb, corr
    if near_correction
        nb = BI.build_neighbor_list(interface, max_order, eps)
        corr = BI.laplace3d_DT_corrections(interface, nb.upsample, nb.adaptive)
    else
        nb = (; upsample = Dict{Tuple{Int, Int}, Int}(), adaptive = Dict{Tuple{Int, Int}, BI.AdaptiveConfig}())
        corr = spzeros(Float64, N, N)
    end
    times["near_assembly"] = time() - t0

    tvec = _diag_coeffs(interface)
    apply = near_correction ?
        (x -> (D_base * x) .+ (corr * x) .+ (tvec .* x)) :
        (x -> (D_base * x) .+ (tvec .* x))
    A = LinearMap{Float64}(apply, N, N)

    # 5. GMRES
    t0 = time()
    sigma, stats = Krylov.gmres(A, rhs; atol = gmres_atol, rtol = eps,
                                itmax = itmax, history = history, verbose = gmres_verbose)
    times["gmres"] = time() - t0
    residual = norm(A * sigma - rhs) / max(norm(rhs), eps_float())
    @printf("    GMRES: niter=%d  solved=%s  relres=%.3e  t=%.1fs  (%s)\n",
            stats.niter, stats.solved, residual, times["gmres"], stats.status)
    flush(stdout)

    p_ups = collect(values(nb.upsample))
    return (
        sys = sys, eps = eps, p = p, r = r, max_order = max_order,
        fmm_tol = fmm_tol, tkm_tol = tkm_tol, tkm_kmax = tkm_kmax,
        interface = interface, screened_vs = screened_vs,
        rhs = rhs, sigma = sigma, stats = stats, residual = residual,
        N = N, n_panels = length(interface.panels),
        niter = stats.niter,
        gmres_history = history ? copy(stats.residuals) : Float64[],
        n_near_pairs = length(nb.upsample), n_adaptive_pairs = length(nb.adaptive),
        p_up_max = isempty(p_ups) ? 0 : maximum(p_ups),
        p_up_mean = isempty(p_ups) ? 0.0 : mean(p_ups),
        corr_nnz = nnz(corr), corr_bytes = Base.summarysize(corr),
        nb = nb, corr = corr, D_base = D_base, tvec = tvec,
        times = times,
    )
end

eps_float() = Base.eps(Float64)

# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

"""
Incident potential of the screened source at targets. TKM only for targets
near/inside the source support; direct smooth-kernel summation outside it
(TKM's Fourier domain covers sources AND targets, so far targets would blow
its grid up; outside the support the integrand is smooth and the direct sum
of the grid charges carries the same quadrature accuracy).
"""
function eval_incident(res, targets::Matrix{Float64})
    vs = res.screened_vs
    pos = vs.positions
    q = vs.weights .* vs.density
    h = BI._estimate_source_spacing(vs)
    lo = (minimum(@view pos[1, :]), minimum(@view pos[2, :]), minimum(@view pos[3, :])) .- 3h
    hi = (maximum(@view pos[1, :]), maximum(@view pos[2, :]), maximum(@view pos[3, :])) .+ 3h
    n = size(targets, 2)
    far = [!(lo[1] <= targets[1, j] <= hi[1] &&
             lo[2] <= targets[2, j] <= hi[2] &&
             lo[3] <= targets[3, j] <= hi[3]) for j in 1:n]
    out = Vector{Float64}(undef, n)
    near_ids = findall(!, far)
    if !isempty(near_ids)
        vals = BoundaryIntegral.TKM3D.ltkm3dc(
            res.tkm_tol, pos;
            charges = q, targets = targets[:, near_ids], pgt = 1, kmax = res.tkm_kmax)
        vals.ier == 0 || error("TKM3D.ltkm3dc failed with ier=$(vals.ier)")
        out[near_ids] .= real.(vals.pottarg)
    end
    far_ids = findall(far)
    if !isempty(far_ids)
        ns = length(q)
        Base.Threads.@threads for jj in eachindex(far_ids)
            j = far_ids[jj]
            acc = 0.0
            @inbounds for s in 1:ns
                dx = targets[1, j] - pos[1, s]
                dy = targets[2, j] - pos[2, s]
                dz = targets[3, j] - pos[3, s]
                acc += q[s] / sqrt(dx * dx + dy * dy + dz * dz)
            end
            out[j] = acc / (4pi)
        end
    end
    return out
end

"Scattered potential S[sigma] at targets (post-refined FMM + hcubature near)."
function eval_scatter(res, targets::Matrix{Float64})
    op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
        res.interface, targets, res.fmm_tol, res.eps, 5.0)
    return op * res.sigma
end

"Total potential phi = u_inc + S[sigma] at 3xM targets."
function eval_phi(res, targets::Matrix{Float64}; t_out::Union{Nothing, Dict{String, Float64}} = nothing)
    t0 = time()
    u_inc = eval_incident(res, targets)
    t1 = time()
    u_sc = eval_scatter(res, targets)
    t2 = time()
    if t_out !== nothing
        t_out["eval_incident"] = get(t_out, "eval_incident", 0.0) + (t1 - t0)
        t_out["eval_scatter"] = get(t_out, "eval_scatter", 0.0) + (t2 - t1)
    end
    return u_inc .+ u_sc
end

"V = int rho_tgt phi  (target density NOT screened)."
function eval_V(res; tgt_vs::Union{Nothing, BI.VolumeSource{Float64, 3}} = nothing,
                t_out = nothing, margin::Float64 = 1.25)
    vs_t = isnothing(tgt_vs) ?
        gaussian_source(res.sys.tgt_center, res.sys.tgt_sigma, res.eps; margin = margin) : tgt_vs
    phi = eval_phi(res, Matrix{Float64}(vs_t.positions); t_out = t_out)
    return dot(vs_t.weights .* vs_t.density, phi)
end

"Vacuum interaction V (no interfaces): direct TKM between unscreened densities."
function vacuum_V(sys::SystemSpec, eps::Float64)
    vs_s = gaussian_source(sys.src_center, sys.src_sigma, eps)
    vs_t = gaussian_source(sys.tgt_center, sys.tgt_sigma, eps)
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(vs_s))
    vals = BoundaryIntegral.TKM3D.ltkm3dc(
        max(eps, 1e-13), vs_s.positions;
        charges = vs_s.weights .* vs_s.density,
        targets = Matrix{Float64}(vs_t.positions), pgt = 1, kmax = kmax)
    vals.ier == 0 || error("TKM ier=$(vals.ier)")
    return dot(vs_t.weights .* vs_t.density, real.(vals.pottarg))
end

"""
Scattered-potential evaluation with caller-controlled post-refinement threshold
h0 (Inf disables refinement). Inlines laplace3d_pottrg_fmm3d_corrected_hcubature.
Returns (values, stats) where stats has n_refined, n_hcub, t_refine, t_fmm, t_near.
"""
function eval_scatter_with_h0(res, targets::Matrix{Float64}, h0::Float64;
                              hcub_atol::Float64 = res.eps, range_factor::Float64 = 5.0)
    interface = res.interface
    n_points = BI.num_points(interface)
    t0 = time()
    if isfinite(h0)
        refined, parent_ids, from_split = BI._refine_interface_for_targets(
            interface, targets, h0; range_factor = range_factor)
        prol = BI._refined_interface_prolongation(interface, refined, parent_ids, from_split)
    else
        refined = interface
        prol = sparse(1:n_points, 1:n_points, ones(n_points), n_points, n_points)
    end
    t_refine = time() - t0

    t0 = time()
    pot_base = BI.laplace3d_pottrg_fmm3d(refined, targets, res.fmm_tol)
    sig_ref = prol * res.sigma
    u = pot_base * sig_ref
    t_fmm = time() - t0

    t0 = time()
    tnl = BI.build_target_neighbor_list(refined, targets, false; range_factor = range_factor)
    n_hcub = isempty(tnl) ? 0 : sum(length(v) for v in values(tnl))
    corr = BI.laplace3d_pottrg_corrections_hcubature(refined, targets, tnl, hcub_atol)
    u .+= corr * sig_ref
    t_near = time() - t0

    return u, (; n_refined = BI.num_points(refined), n_hcub,
                 t_refine, t_fmm, t_near)
end

# ---------------------------------------------------------------------------
# Target zones (fixed RNG)
# ---------------------------------------------------------------------------

"""
Zone target sets for a system (3x200 matrices): :near (0.01 above the top
substrate surface), :support (ball of radius 2*sigma around the source),
:far (sphere of radius 5 around the origin).
"""
function zone_targets(sys::SystemSpec; n::Int = 200, seed::Int = SEED_TARGETS)
    rng = MersenneTwister(seed)
    # near-interface: 0.01 above the top face of the substrate box(es)
    sub = sys.boxes[1]
    xmin = minimum(b.center[1] - b.Lx / 2 for b in sys.boxes if b.center[3] ≈ sub.center[3])
    xmax = maximum(b.center[1] + b.Lx / 2 for b in sys.boxes if b.center[3] ≈ sub.center[3])
    ztop = sub.center[3] + sub.Lz / 2
    near = Matrix{Float64}(undef, 3, n)
    for i in 1:n
        near[1, i] = xmin + (xmax - xmin) * rand(rng)
        near[2, i] = sub.center[2] - sub.Ly / 2 + sub.Ly * rand(rng)
        near[3, i] = ztop + 0.01
    end
    # source support: uniform in ball of radius 2 sigma
    supp = Matrix{Float64}(undef, 3, n)
    for i in 1:n
        v = randn(rng, 3); v ./= norm(v)
        rad = 2 * sys.src_sigma * rand(rng)^(1 / 3)
        supp[:, i] .= collect(sys.src_center) .+ rad .* v
    end
    # far field: sphere at distance 5 beyond the system circumradius
    rfar = 5.0 + maximum(norm(collect(b.center)) + 0.5 * norm([b.Lx, b.Ly, b.Lz])
                         for b in sys.boxes)
    far = Matrix{Float64}(undef, 3, n)
    for i in 1:n
        v = randn(rng, 3); v ./= norm(v)
        far[:, i] .= rfar .* v
    end
    return (; near, supp, far)
end

# ---------------------------------------------------------------------------
# Refinement helpers (variant C, scaling)
# ---------------------------------------------------------------------------

"Split every panel into 4 children, k times (uniform refinement)."
function uniform_refine(interface, k::Int)
    k <= 0 && return interface
    panels = interface.panels
    eps_in = interface.eps_in
    eps_out = interface.eps_out
    for _ in 1:k
        new_panels = eltype(panels)[]
        new_in = Float64[]
        new_out = Float64[]
        for (i, panel) in enumerate(panels)
            for child in BI._split_panel4(panel)
                push!(new_panels, child)
                push!(new_in, eps_in[i])
                push!(new_out, eps_out[i])
            end
        end
        panels, eps_in, eps_out = new_panels, new_in, new_out
    end
    return BI.DielectricInterface(panels, eps_in, eps_out)
end

"Greedy largest-panel-first dyadic refinement until DOF >= target_dof (quasi-uniform)."
function refine_to_dof(interface, target_dof::Int)
    panels = collect(interface.panels)
    eps_in = collect(interface.eps_in)
    eps_out = collect(interface.eps_out)
    lens = [BI._panel_max_length(p) for p in panels]
    dof = sum(BI.num_points(p) for p in panels)
    while dof < target_dof
        i = argmax(lens)
        children = BI._split_panel4(panels[i])
        dof += 3 * BI.num_points(panels[i])
        ei, eo = eps_in[i], eps_out[i]
        panels[i] = children[1]
        lens[i] = BI._panel_max_length(children[1])
        for c in children[2:4]
            push!(panels, c); push!(eps_in, ei); push!(eps_out, eo)
            push!(lens, BI._panel_max_length(c))
        end
    end
    return BI.DielectricInterface(panels, eps_in, eps_out)
end

# ---------------------------------------------------------------------------
# Bernstein stagnation estimate (6.2 annotation)
# ---------------------------------------------------------------------------

"""
Minimum Bernstein radius over near pairs: for each upsample pair (i,j), map the
closest quadrature point of panel j into panel i's local [-1,1]^2 coordinates
and take rho of the 1D pole along the more singular axis. Stagnation ~ rho^(-2p).
"""
function bernstein_rho_min(interface, nb)
    rho_min = Inf
    for (i, j) in keys(nb.upsample)
        pi_ = interface.panels[i]
        pj = interface.panels[j]
        a, b, c, d = pi_.corners
        cc = (a .+ b .+ c .+ d) ./ 4
        bma = b .- a; dma = d .- a
        Lu = norm(bma); Lv = norm(dma)
        eu = bma ./ Lu; ev = dma ./ Lv
        nrm = pi_.normal
        for q in pj.points
            x = q .- cc
            u = 2 * dot(x, eu) / Lu
            v = 2 * dot(x, ev) / Lv
            w = dot(x, nrm)
            for (s, L) in ((u, Lu), (v, Lv))
                z = complex(s, 2 * abs(w) / L)
                rho = abs(z + sqrt(z^2 - 1))
                rho = max(rho, abs(z - sqrt(z^2 - 1)))
                rho_min = min(rho_min, rho)
            end
        end
    end
    return rho_min
end

# ---------------------------------------------------------------------------
# Logging / provenance / reference snapshots
# ---------------------------------------------------------------------------

function run_provenance()
    return (
        hostname = gethostname(),
        nthreads = Base.Threads.nthreads(),
        omp = get(ENV, "OMP_NUM_THREADS", ""),
        julia = string(VERSION),
        date = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
    )
end

"Append a NamedTuple as a CSV row (writes header if the file is new). Crash-safe."
function append_csv_row(path::AbstractString, row::NamedTuple)
    mkpath(dirname(path))
    newfile = !isfile(path)
    open(path, "a") do io
        if newfile
            println(io, join(string.(keys(row)), ","))
        end
        println(io, join([sprint(print, v) for v in values(row)], ","))
    end
end

times_cols(times::Dict{String, Float64}) = (
    t_source = get(times, "source_prep", 0.0),
    t_interface = get(times, "interface_build", 0.0),
    t_rhs = get(times, "rhs_assembly", 0.0),
    t_fmm_setup = get(times, "fmm_setup", 0.0),
    t_near = get(times, "near_assembly", 0.0),
    t_gmres = get(times, "gmres", 0.0),
    t_eval_inc = get(times, "eval_incident", 0.0),
    t_eval_sc = get(times, "eval_scatter", 0.0),
)

"Columns common to every run log."
function run_cols(res; extra...)
    pv = run_provenance()
    return merge((
        system = res.sys.name, eps = res.eps, p = res.p, r = res.r,
        N = res.N, n_panels = res.n_panels, niter = res.niter,
        residual = res.residual,
        n_near_pairs = res.n_near_pairs, n_adaptive_pairs = res.n_adaptive_pairs,
        p_up_max = res.p_up_max, p_up_mean = round(res.p_up_mean; digits = 3),
        corr_nnz = res.corr_nnz, corr_mb = round(res.corr_bytes / 1e6; digits = 2),
        max_order = res.max_order, fmm_tol = res.fmm_tol, tkm_kmax = res.tkm_kmax,
        hostname = pv.hostname, nthreads = pv.nthreads, date = pv.date,
    ), times_cols(res.times), NamedTuple(extra))
end

function save_ref(path::AbstractString, payload)
    mkpath(dirname(path))
    serialize(path, payload)
end

load_ref(path::AbstractString) = deserialize(path)

end # module
