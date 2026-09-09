# Single-orbital l_ec self-convergence: one representative orbital (nearest the junction centre),
# EXACT benchmark geometry (×3 het cubes + tight ε=10 slab), all solver tolerances = 1e-5,
# neighbor_cutoff=0 (on-site only). Sweeps edge_refine_level ∈ {1,2,3,4} → l_ec ∈ {4.55, 2.27,
# 1.14, 0.57} and reports the on-site U (eV), interface DOF, and solve time at each. Cheap
# stand-in for the full-campaign l_ec convergence.
#
# Run (worker7011, FINUFFT-capped): TKM3D_FINUFFT_NTHREADS=16 julia --project -t8 scripts/lec_conv_single.jl
using BoundaryIntegral, LinearAlgebra, Printf, TOML
const BI = BoundaryIntegral
const E2 = 14.3996
Timer(_ -> (flush(stdout); flush(stderr)), 2; interval = 2)

const LS = normpath(joinpath(@__DIR__, ".."))
# RERUN_TAG suffixes the campaign name (hence its ceph root) and the output TSV.
const TAG = get(ENV, "RERUN_TAG", "")
# RERUN_BASE picks which campaign supplies templates/geometry/dielectrics. The eps=2.4 sweep must
# inherit that campaign's boxes, not the eps=10 ones.
const BASE = get(ENV, "RERUN_BASE", "lattice_10x10_het3x")
# DRY_RUN=1 writes the campaign TOMLs and stops before any solve (used to validate inputs).
const DRY = get(ENV, "DRY_RUN", "") == "1"
const LEVELS = isempty(ARGS) ? [1, 2, 3, 4] : parse.(Int, ARGS)
const XJ, YC = 5.547, 10.318

bench = TOML.parsefile(joinpath(LS, "campaigns", BASE * ".toml"))
templates = bench["templates"]; boxes = bench["dielectrics"]["boxes"]; orbs = bench["orbital"]
_, oi = findmin([hypot(o["x"] - XJ, o["y"] - YC) for o in orbs])
orb = orbs[oi]
@info "chosen orbital" idx = oi type = orb["type"] x = orb["x"] y = orb["y"]

function write_toml(level)
    nm = "lec_single_l$(level)$(TAG)"
    io = IOBuffer()
    println(io, "name = \"$nm\"")
    println(io, "root = \"/mnt/ceph/users/xgao1/four_index/$nm\"")
    println(io, "templates = ["); for t in templates; println(io, "  \"$t\","); end; println(io, "]\n")
    println(io, "[[orbital]]"); println(io, "type = $(orb["type"])")
    println(io, "x = $(orb["x"])"); println(io, "y = $(orb["y"])"); println(io, "z = $(orb["z"])\n")
    println(io, "[pairing]"); println(io, "neighbor_cutoff = 0.0\n")
    println(io, "[dielectrics]"); println(io, "eps_out = 1.0"); println(io, "boxes = [")
    for b in boxes; println(io, "  [", join(b, ", "), "],"); end; println(io, "]\n")
    println(io, "[solve]"); println(io, "n_quad = 6"); println(io, "edge_refine_level = $level")
    println(io, "rhs_tol = 1e-5"); println(io, "lhs_tol = 1e-5"); println(io, "gmres_rtol = 1e-5")
    println(io, "support_rtol = 1e-4"); println(io, "volume_tol = 1e-5")
    # max_order 64, not 8 -- the same fix the full-campaign generators took. p_up is clamped to
    # max_order with no adaptive fallback for non-touching pairs and no warning, and at 8 it left
    # 92.5% of a real conv_l3 interface's 12,208 upsample pairs pinned at the cap. The l_ec sweep
    # CANNOT detect that on its own: the bias is present at every refinement level, so U converges
    # cleanly to a slightly wrong value. Its flatness is only evidence of l_ec convergence.
    println(io, "max_order = 64"); println(io, "max_depth = 128\n")
    println(io, "[batching]"); println(io, "n_centers_per_batch = 1\n")
    println(io, "[eval]"); println(io, "c_pad = 5.0")
    path = joinpath(LS, "campaigns", "$nm.toml"); write(path, String(take!(io))); return path
end

function onsite_U(c)
    br = BI.load_batch_result(BI.batch_path(c, 1))
    dg = BI.load_templates!(c)[1][2]
    pos = BI.grid_positions(dg, br.gidx)
    Φ = BI.evaluate_batch_potential(br.interface, br.sigma,
            [BI.VolumeSource(copy(pos), copy(br.weights), br.densities[:, 1])], pos;
            lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], c_pad = c.c_pad,
            screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
    wρ = br.weights .* br.densities[:, 1]
    return dot(wρ, Φ[:, 1]) * 4π * E2 / sum(wρ)^2, size(br.sigma, 1)
end

# Merge results into the TSV, keyed by level, and write it out. Called after EVERY level, not
# only at the end: a later level that is OOM-killed or hits the walltime would otherwise discard
# the completed ones with it (level 6 alone is ~7 min of solve, level 7 considerably more).
const TSV = joinpath(LS, "figs", "lec_single_conv$(TAG).tsv")
function flush_tsv(results)
    rows = Dict{Int,Tuple{Float64,Float64,Int,Float64}}()   # level → (l_ec, U, dof, t)
    if isfile(TSV)
        for ln in eachline(TSV)
            (isempty(ln) || startswith(ln, "level")) && continue
            f = split(ln, '\t')
            rows[parse(Int, f[1])] = (parse(Float64,f[2]), parse(Float64,f[3]), parse(Int,f[4]), parse(Float64,f[5]))
        end
    end
    for r in results; rows[r[1]] = (r[2], r[3], r[4], r[5]); end
    open(TSV, "w") do io
        println(io, "level\tl_ec\tU_eV\tdof\tt_solve_s")
        for lv in sort(collect(keys(rows))); r = rows[lv]; @printf(io, "%d\t%.6g\t%.8g\t%d\t%.1f\n", lv, r...); end
    end
    println("wrote $TSV")
end

results = Tuple{Int,Float64,Float64,Int,Float64}[]
for lv in LEVELS
    lec = 9.0 / 2^lv * 1.01
    @info "=== edge_refine_level $lv  (l_ec = $(round(lec;digits=3))) ==="
    c = load_campaign(write_toml(lv))
    if DRY
        @info "dry run: campaign written, skipping solve" lv toml = c.toml_path
        continue
    end
    prepare(c)
    t0 = time(); solve_batch(c, 1); ts = time() - t0
    U, dof = onsite_U(c)
    push!(results, (lv, lec, U, dof, ts))
    @info "level done" lv lec U dof t_solve = round(ts; digits = 1)
    DRY || flush_tsv(results)      # checkpoint now; do not gamble it on later levels
end

println("\n=== single-orbital l_ec convergence (orbital idx=$oi, type $(orb["type"])) ===")
println("level\tl_ec\tU_eV\t\tdof\tt_solve(s)\tΔU")
for (k, r) in enumerate(results)
    du = k > 1 ? abs(r[3] - results[k-1][3]) : NaN
    @printf("%d\t%.3f\t%.6f\t%d\t%.0f\t\t%.2e\n", r[1], r[2], r[3], r[4], r[5], du)
end
# merge with any existing rows (so separate level runs accumulate), keyed by level
if DRY
    @info "dry run: campaigns written; leaving the TSV untouched"
    exit(0)
end
flush_tsv(results)
