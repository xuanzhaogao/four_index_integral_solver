# Orbital "1 slab on 1 cube" onsite U (in eV) at a given system size + cube ε.
# The REAL graphene orbital version of scripts/pointcharge_control.jl `singlecube` mode:
# same geometry (true cube substrate scaled in all 3 dims, top face fixed at z=3; ε=10 graphene
# slab z∈[3,12]; orbital at slab center (0,0,7.5)), but the source is the actual orbital density
# and the observable is the full (direct+induced) screened onsite U in eV — mirrors the original
# scripts/standalone_singlecube.jl onsite_U.
#
# Usage:  julia --project -t <n> scripts/orbital_singlecube.jl <scale> <eps_cube>
#   e.g.  julia --project -t 96 scripts/orbital_singlecube.jl 2.0 11.9
# One case per invocation → driven as a Slurm array (jobscripts/orbital_singlecube.sbatch).

using BoundaryIntegral, LinearAlgebra, Printf
const BI = BoundaryIntegral
const E2 = 14.3996
const TEMPLATE = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/graphene_00001.xsf"

# Julia 1.12 fully buffers stdout/stderr to a file — flush every 2 s to keep the log live.
Timer(_ -> (flush(stdout); flush(stderr)), 2; interval = 2)

function write_toml(scale, eps)
    L = 90.0 * scale                 # cube edge (all 3 dims)
    czc = 3.0 - L / 2                 # cube center-z so its top face sits at z=3 (meets slab)
    Lxy = 90.0 * scale               # slab lateral size (thickness stays 9, z∈[3,12])
    nm = "orbital_singlecube_eps$(eps)_s$(scale)"
    s = """
    name = "$(nm)"
    root = "/mnt/ceph/users/xgao1/four_index/$(nm)"
    templates = ["$(TEMPLATE)"]

    [[orbital]]
    type = 1
    x = 0.0
    y = 0.0
    z = 7.5

    [pairing]
    neighbor_cutoff = 0.0

    [dielectrics]
    eps_out = 1.0
    boxes = [
      [0.0, 0.0, $(czc), $(L), $(L), $(L), $(eps)],
      [0.0, 0.0, 7.5, $(Lxy), $(Lxy), 9.0, 10.0],
    ]

    [solve]
    n_quad = 6
    edge_refine_level = 2
    rhs_tol = 1e-3
    lhs_tol = 1e-5
    gmres_rtol = 1e-5
    support_rtol = 1e-4
    volume_tol = 1e-5
    max_order = 8
    max_depth = 128

    [batching]
    n_centers_per_batch = 1

    [eval]
    far_pad_steps = 2.0
    """
    path = joinpath(@__DIR__, "..", "campaigns", "$(nm).toml")
    write(path, s)
    return path
end

# Full (direct + induced) screened onsite U in eV — identical to standalone_singlecube.jl.
function onsite_U(c)
    br = BI.load_batch_result(BI.batch_path(c, 1))
    dg = BI.load_templates!(c)[1][2]
    pos = BI.grid_positions(dg, br.gidx)
    At, Bt, Ct = BI.true_cell_vectors(dg)
    far_pad = c.far_pad_steps * maximum((norm(collect(At)) / dg.nx,
                                         norm(collect(Bt)) / dg.ny,
                                         norm(collect(Ct)) / dg.nz))
    Φ = BI.evaluate_batch_potential(br.interface, br.sigma,
            [BI.VolumeSource(copy(pos), copy(br.weights), br.densities[:, 1])], pos;
            lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], far_pad = far_pad,
            screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
    wρ = br.weights .* br.densities[:, 1]
    return dot(wρ, Φ[:, 1]) * 4π * E2 / sum(wρ)^2
end

function main()
    length(ARGS) >= 2 || error("usage: orbital_singlecube.jl <scale> <eps_cube>")
    scale = parse(Float64, ARGS[1])
    eps   = parse(Float64, ARGS[2])
    @info "orbital singlecube" scale eps threads=Threads.nthreads() finufft_cap=get(ENV, "TKM3D_FINUFFT_NTHREADS", "(default 16)")

    c = load_campaign(write_toml(scale, eps))
    prepare(c)
    t0 = time()
    solve_batch(c, 1)
    @info "solve done" t_solve=round(time() - t0; digits = 1)
    t1 = time()
    U = onsite_U(c)
    @info "eval done" t_eval=round(time() - t1; digits = 1)
    @info "RESULT" scale eps U_eV=U

    out = joinpath(@__DIR__, "..", "figs", "orbital_singlecube_eps$(eps)_s$(scale).tsv")
    open(out, "w") do io
        println(io, "eps\tscale\tU_eV")
        @printf(io, "%.4g\t%.4g\t%.10g\n", eps, scale, U)
    end
    @info "wrote" out
end

main()
