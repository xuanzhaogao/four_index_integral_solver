# Point-charge screening probe + single-ε control (spec:
# docs/superpowers/specs/2026-07-09-pointcharge-screening-control-design.md).
#
# Places a unit point charge at each onsite_grid site (z=7.5) over a UNIFORM-ε substrate
# (geometry-identical to the het union: one 180×90×90 box) + the ε=10 slab, solves the BEM
# (ONE shared interface, multi-RHS block), and reports the induced reaction-field self-energy
# U_proxy = Φ_ind(r_a)·4π·E2 (diagonal of the induced-potential matrix; direct term excluded).
# If a uniform substrate ramps in x → the het U-ramp is a finite-size artifact; if it is flat
# with only a symmetric y-style edge bowl → the het ramp is real ε-contrast screening.
#
# Run (on worker7011):
#   julia --project=codes/lattice_scale -t <n> codes/lattice_scale/scripts/pointcharge_control.jl pilot
#   julia --project=codes/lattice_scale -t <n> codes/lattice_scale/scripts/pointcharge_control.jl full

using BoundaryIntegral, LinearAlgebra, Printf
const BI = BoundaryIntegral

# Julia 1.12 fully buffers stdout/stderr when redirected to a file — flush every 2 s so the
# log stays live during long runs (mirrors driver.jl).
Timer(_ -> (flush(stdout); flush(stderr)), 2; interval = 2)
const E2 = 14.3996
const GRID_TSV = "/mnt/ceph/users/xgao1/four_index/onsite_grid/onsite_U.tsv"

# solve params — verbatim from campaigns/onsite_grid.toml
const NQUAD = 6; const EDGE_LEVEL = 2; const RHS_TOL = 1e-3; const LHS_TOL = 1e-5
const GMRES_RTOL = 1e-5; const MAX_ORDER = 8; const MAX_DEPTH = 128
const RANGE_FACTOR = 5.0; const CHUNK = 128
const EPSES_FULL = (3.9, 11.9)

function read_grid(path)
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(ln, '\t')
        push!(xs, parse(Float64, f[2])); push!(ys, parse(Float64, f[3])); push!(zs, parse(Float64, f[4]))
    end
    return permutedims(hcat(xs, ys, zs))          # 3 × N
end

function build_boxes(eps_sub; scale = 1.0, cube = false)
    # cube=true: a genuine cube substrate scaled in ALL 3 dims (Lx=Ly=Lz=90·scale), top face
    #   fixed at z=3 (meets the slab bottom), so it grows outward AND downward — the "1 slab on
    #   1 cube" building block. cube=false: substrate spans the het union footprint (180 wide,
    #   Lz=90 fixed) — the geometry used by the size/cuts/full modes.
    if cube
        L = 90.0 * scale
        sub = (center = (5.547, 10.318, 3.0 - L / 2), Lx = L, Ly = L, Lz = L)  # top at z=3
    else
        sub = (center = (5.547, 10.318, -42.0), Lx = 180.0 * scale, Ly = 90.0 * scale, Lz = 90.0)
    end
    boxes = BI.BoxGeom[
        sub,
        (center = (5.547, 10.318, 7.5), Lx = 90.0 * scale, Ly = 90.0 * scale, Lz = 9.0),  # ε=10 slab, z∈[3,12]
    ]
    epses = [eps_sub, 10.0]
    return boxes, epses
end

l_ec_rule(boxes) = minimum(b.Lz for b in boxes) / 2.0^EDGE_LEVEL * 1.01

# Build the shared interface, LHS operator, and corrected pottrg map for target set `Ptgt`.
# The interface + operator depend only on geometry, so they are built ONCE and reused for
# every charge (the whole cost saving vs the orbital pipeline).
function build_shared(eps_sub, Penv, Ptgt; scale = 1.0, rhs_tol = RHS_TOL, edge_level = EDGE_LEVEL, cube = false)
    boxes, epses = build_boxes(eps_sub; scale = scale, cube = cube)
    eps_out = 1.0
    l_ec = minimum(b.Lz for b in boxes) / 2.0^edge_level * 1.01
    N = size(Penv, 2)
    env = BI.VolumeSource(copy(Penv), ones(N), ones(N))
    @info "build interface" eps_sub l_ec rhs_tol edge_level N_env=N
    t0 = time()
    interface = BI.multi_dielectric_box3d_rhs_adaptive(NQUAD, l_ec, boxes, epses, env, rhs_tol;
        eps_out = eps_out, max_depth = MAX_DEPTH)
    t_iface = time() - t0
    t1 = time()
    op = BI.batched_lhs_dielectric_box3d_fmm3d_corrected(interface, LHS_TOL, LHS_TOL, MAX_ORDER)
    t_op = time() - t1
    t2 = time()
    pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, Ptgt, LHS_TOL, LHS_TOL, RANGE_FACTOR)
    @info "interface + op + pottrg built" dof=size(pottrg, 2) nt=size(pottrg, 1) t_iface t_op t_pottrg=time()-t2
    return interface, op, pottrg, boxes, epses, eps_out
end

# U_proxy for the charges at columns `idx` of P. `pottrg`/target rows are aligned with idx order.
# The prebuilt LHS `op` is reused for every chunk; only the RHS + block-GMRES vary.
function solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
    U = fill(NaN, length(idx))
    for c0 in 1:CHUNK:length(idx)
        c1 = min(c0 + CHUNK - 1, length(idx))
        sub = idx[c0:c1]
        sources = [BI.VolumeSource(reshape(P[:, a], 3, 1), [1.0], [1.0]) for a in sub]
        tc = time()
        F = BI.rhs_dielectric_box3d_fmm3d(interface, sources, LHS_TOL;
            screen_boxes = boxes, screen_epses = epses, screen_eps_out = eps_out)
        Σ, stats = BI._block_gmres_solve(op, F; rtol = GMRES_RTOL, atol = 0.0, itmax = 500)
        for li in 1:length(sub)
            row = c0 + li - 1                       # pottrg row aligned with idx order
            phi = pottrg * Σ[:, li]                 # induced potential at ALL targets (lazy LinearMap apply)
            U[row] = phi[row] * 4π * E2             # self target = this charge's own location
        end
        @info "chunk" range="$(c0):$(c1)" niter=stats.niter t=round(time() - tc; digits = 1)
    end
    return U
end

function probe_indices(P)
    x = P[1, :]; y = P[2, :]
    cx, cy = sum(x) / length(x), sum(y) / length(y)
    center = argmin(@. hypot(x - cx, y - cy))
    return unique([center, argmin(x), argmax(x), argmin(y), argmax(y)])
end

# Two orthogonal cuts through the domain centroid: a horizontal band (probes U vs x at ~const y)
# and a vertical band (probes U vs y at ~const x). Their union answers the artifact question
# cheaply — the same axes used to analyze the heterojunction map. Returns (idx, tag) with tag ∈
# {:h, :v} per selected site (a site in both bands is tagged :h).
function cut_indices(P; halfwidth = 1.5)
    x = P[1, :]; y = P[2, :]
    cx = sum(x) / length(x); cy = sum(y) / length(y)
    hmask = abs.(y .- cy) .< halfwidth      # horizontal cut: ~const y
    vmask = abs.(x .- cx) .< halfwidth       # vertical cut:   ~const x
    idx = findall(hmask .| vmask)
    tag = [hmask[i] ? :h : :v for i in idx]
    return idx, tag, cx, cy
end

function main()
    mode = isempty(ARGS) ? "pilot" : ARGS[1]
    P = read_grid(GRID_TSV)
    @info "grid" N=size(P, 2) xspan=extrema(P[1, :]) yspan=extrema(P[2, :]) z=P[3, 1]

    if mode == "pilot"
        idx = probe_indices(P)
        Ptgt = P[:, idx]
        interface, op, pottrg, boxes, epses, eps_out = build_shared(11.9, P, Ptgt)
        U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
        println("\n=== PILOT (ε=11.9) probe self-energies ===")
        println("idx\tx\ty\tU_proxy")
        for (r, a) in enumerate(idx)
            @printf("%d\t%.2f\t%.2f\t%.5f\n", a, P[1, a], P[2, a], U[r])
        end
    elseif mode == "cuts"
        idx, tag, cx, cy = cut_indices(P)
        Ptgt = P[:, idx]
        @info "cuts" n=length(idx) n_h=count(==(:h), tag) n_v=count(==(:v), tag) cx cy
        outdir = joinpath(@__DIR__, "..", "figs")
        for eps_sub in EPSES_FULL
            interface, op, pottrg, boxes, epses, eps_out = build_shared(eps_sub, P, Ptgt)
            U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
            out = joinpath(outdir, "pointcharge_cuts_eps$(eps_sub).tsv")
            open(out, "w") do io
                println(io, "orbital\tcut\tx\ty\tz\tU_proxy")
                for (r, a) in enumerate(idx)
                    @printf(io, "%d\t%s\t%.10g\t%.10g\t%.10g\t%.10g\n",
                            a, tag[r], P[1, a], P[2, a], P[3, a], U[r])
                end
            end
            @info "wrote" out Umin=minimum(U) Umax=maximum(U)
        end
    elseif mode == "size"
        # Size-scaling test: a horizontal line of probes across x (mid-y), recomputed for a few
        # system sizes. Does the interior converge (edge upturn just moves out → boundary effect,
        # bigger box helps) or does the center itself drift (whole box too small)?
        eps_sub = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3.9
        scales = length(ARGS) >= 3 ? Tuple(parse(Float64, a) for a in ARGS[3:end]) : (1.0, 2.0, 3.0)
        x = P[1, :]; y = P[2, :]
        cy = sum(y) / length(y)
        line = findall(abs.(y .- cy) .< 1.5)
        line = line[sortperm(x[line])]
        idx = line[round.(Int, range(1, length(line); length = 15))]   # ~15 probes across x
        Ptgt = P[:, idx]
        @info "size test" eps_sub scales n_probes=length(idx) xspan=extrema(x[idx]) cy
        outdir = joinpath(@__DIR__, "..", "figs")
        for s in scales
            interface, op, pottrg, boxes, epses, eps_out = build_shared(eps_sub, P, Ptgt; scale = s)
            U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
            out = joinpath(outdir, "pointcharge_size_eps$(eps_sub)_s$(s).tsv")
            open(out, "w") do io
                println(io, "orbital\tx\ty\tU_proxy")
                for (r, a) in enumerate(idx)
                    @printf(io, "%d\t%.10g\t%.10g\t%.10g\n", a, P[1, a], P[2, a], U[r])
                end
            end
            @info "size done" scale=s Ucenter=U[fld(length(U) + 1, 2)] Uedge_lo=U[1] Uedge_hi=U[end]
        end
    elseif mode == "resconv"
        # Resolution-convergence at FIXED box size: refine the mesh (rhs_tol + edge_level) and
        # watch U on the probe line. If U drifts toward the same value the size test converged to,
        # the earlier "size effect" was really under-resolution (best_grid tied panel size to box
        # size), NOT finite size. `scale` fixed (default 1.0 = the onsite_grid baseline box).
        eps_sub = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3.9
        scale = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0
        # (rhs_tol, edge_level) rungs, coarse → fine
        rungs = [(1e-3, 2), (1e-4, 3), (1e-5, 4)]
        x = P[1, :]; y = P[2, :]
        cy = sum(y) / length(y)
        line = findall(abs.(y .- cy) .< 1.5)
        line = line[sortperm(x[line])]
        idx = line[round.(Int, range(1, length(line); length = 9))]   # ~9 probes across x
        Ptgt = P[:, idx]
        @info "resconv" eps_sub scale rungs n_probes=length(idx)
        outdir = joinpath(@__DIR__, "..", "figs")
        for (rt, el) in rungs
            interface, op, pottrg, boxes, epses, eps_out =
                build_shared(eps_sub, P, Ptgt; scale = scale, rhs_tol = rt, edge_level = el)
            U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
            out = joinpath(outdir, "pointcharge_resconv_eps$(eps_sub)_s$(scale)_rt$(rt)_el$(el).tsv")
            open(out, "w") do io
                println(io, "orbital\tx\ty\tU_proxy")
                for (r, a) in enumerate(idx)
                    @printf(io, "%d\t%.10g\t%.10g\t%.10g\n", a, P[1, a], P[2, a], U[r])
                end
            end
            @info "resconv rung done" rhs_tol=rt edge_level=el dof=size(pottrg, 2) Ucenter=U[fld(length(U) + 1, 2)] Uedge=U[1]
        end
    elseif mode == "singlecube"
        # "1 slab on 1 cube": genuine single cube substrate + ε=10 slab, unit point charge at the
        # box CENTER, swept over ε, at a large-enough (default ×2) size — the converged bulk
        # single-cube onsite-U reference (point-charge proxy). Parallels standalone_singlecube.jl.
        scale = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.0
        epses_sc = (3.9, 11.9)          # cube: SiO₂ (3.9) and Si (11.9); slab stays ε=10 (graphene)
        c0 = (5.547, 10.318, 7.5)
        Pc = reshape(collect(c0), 3, 1)                 # single charge / target at box center
        x = P[1, :]; y = P[2, :]
        envmask = (@. hypot(x - c0[1], y - c0[2])) .< 15.0
        Penv = P[:, envmask]                            # local env cloud → interface refinement near center
        @info "singlecube" scale n_env=size(Penv, 2) center=c0
        results = Tuple{Float64,Float64}[]
        for eps_sub in epses_sc
            interface, op, pottrg, boxes, epses, eps_out =
                build_shared(eps_sub, Penv, Pc; scale = scale, cube = true)
            U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, Pc, [1])
            push!(results, (eps_sub, U[1]))
            @info "singlecube done" eps_sub U=U[1]
        end
        println("\n=== single-cube (1 slab on 1 cube) point-charge U_proxy, scale=$(scale) ===")
        println("eps\tU_proxy")
        for (e, u) in results; @printf("%.1f\t%.6f\n", e, u); end
        out = joinpath(@__DIR__, "..", "figs", "pointcharge_singlecube_s$(scale).tsv")
        open(out, "w") do io
            println(io, "eps\tU_proxy")
            for (e, u) in results; @printf(io, "%.4g\t%.10g\n", e, u); end
        end
        @info "wrote" out
    elseif mode == "hetsweep"
        # The original 2-cube heterojunction (Si | SiO₂ meeting at x=0, graphene ε=10 slab
        # straddling the junction), scaled up so each cube is 90·scale cubed (default ×3 = 270³,
        # slab 270×270×9). Sweep a point-charge line across x at y=0, z=7.5 (slab center) and
        # report the induced self-energy from the Si bulk, through the junction, into the SiO₂
        # bulk — now at converged size, so plateaus/transition aren't finite-size contaminated.
        scale = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 3.0
        npts  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 31
        L = 90.0 * scale
        boxes = BI.BoxGeom[
            (center = (-L / 2, 0.0, 3.0 - L / 2), Lx = L, Ly = L, Lz = L),    # Si  (ε=11.9), x∈[-L,0]
            (center = ( L / 2, 0.0, 3.0 - L / 2), Lx = L, Ly = L, Lz = L),    # SiO₂(ε=3.9),  x∈[0,L]
            (center = (0.0,    0.0, 7.5),         Lx = L, Ly = L, Lz = 9.0),  # slab (ε=10),  x∈[-L/2,L/2]
        ]
        epses = [11.9, 3.9, 10.0]
        eps_out = 1.0
        l_ec = minimum(b.Lz for b in boxes) / 2.0^EDGE_LEVEL * 1.01
        half = L / 2 - 5.0                                    # stay just inside the slab footprint
        xs = collect(range(-half, half; length = npts))
        P = permutedims(hcat(xs, zeros(npts), fill(7.5, npts)))          # 3×npts sweep charges/targets
        # 2D env strip along the line so the interface refines over each charge's ~2D footprint
        exs = collect(range(-half, half; step = 3.0)); eys = collect(range(-9.0, 9.0; step = 3.0))
        ev = [(ex, ey, 7.5) for ex in exs for ey in eys]
        Penv = permutedims(hcat([p[1] for p in ev], [p[2] for p in ev], [p[3] for p in ev]))
        Nenv = size(Penv, 2)
        @info "hetsweep" scale L junction_x=0.0 n_sweep=npts n_env=Nenv sweep_x=extrema(xs)
        env = BI.VolumeSource(copy(Penv), ones(Nenv), ones(Nenv))
        t0 = time()
        interface = BI.multi_dielectric_box3d_rhs_adaptive(NQUAD, l_ec, boxes, epses, env, RHS_TOL;
            eps_out = eps_out, max_depth = MAX_DEPTH)
        op = BI.batched_lhs_dielectric_box3d_fmm3d_corrected(interface, LHS_TOL, LHS_TOL, MAX_ORDER)
        pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, P, LHS_TOL, LHS_TOL, RANGE_FACTOR)
        @info "interface built" dof=size(pottrg, 2) t=round(time() - t0; digits = 1)
        U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, collect(1:npts))
        out = joinpath(@__DIR__, "..", "figs", "pointcharge_hetsweep_s$(scale).tsv")
        open(out, "w") do io
            println(io, "x\ty\tz\tU_proxy")
            for a in 1:npts
                @printf(io, "%.10g\t%.10g\t%.10g\t%.10g\n", P[1, a], P[2, a], P[3, a], U[a])
            end
        end
        @info "wrote" out U_Si_bulk=U[1] U_junction=U[fld(npts + 1, 2)] U_SiO2_bulk=U[end]
    elseif mode == "full"
        outdir = joinpath(@__DIR__, "..", "figs")
        for eps_sub in EPSES_FULL
            idx = collect(1:size(P, 2))
            interface, op, pottrg, boxes, epses, eps_out = build_shared(eps_sub, P, P)
            U = solve_selfenergy(interface, op, pottrg, boxes, epses, eps_out, P, idx)
            out = joinpath(outdir, "pointcharge_ctrl_eps$(eps_sub).tsv")
            open(out, "w") do io
                println(io, "orbital\tx\ty\tz\tU_proxy")
                for a in idx
                    @printf(io, "%d\t%.10g\t%.10g\t%.10g\t%.10g\n", a, P[1, a], P[2, a], P[3, a], U[a])
                end
            end
            @info "wrote" out Umin=minimum(U) Umax=maximum(U) Umean=sum(U) / length(U)
        end
    else
        error("unknown mode $mode (use: pilot | full)")
    end
end

main()
