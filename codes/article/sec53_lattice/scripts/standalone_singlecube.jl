# Standalone BULK-U test: an eps=10 slab on a SINGLE 90x90x90 cube, one graphene orbital at
# the CENTER of the system (0,0,7.5) — 45 Å from every lateral edge, so maximally bulk,
# no junction. Sweeps the cube permittivity eps in {3.9, 8.4, 11.9} and reports the onsite
# U (eV). In-process per eps: prepare -> solve_batch(1) -> local onsite eval.
#
# Run: julia --project=codes/article -t8 codes/article/sec53_lattice/scripts/standalone_singlecube.jl

using BoundaryIntegral, Serialization, LinearAlgebra, CairoMakie
const BI = BoundaryIntegral
const E2 = 14.3996
const TEMPLATE = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/graphene_00001.xsf"
const EPSES = (3.9, 8.4, 11.9)

function write_toml(eps)
    nm = "singlecube_eps$(eps)"
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
      [0.0, 0.0, -42.0, 90.0, 90.0, 90.0, $(eps)],
      [0.0, 0.0, 7.5, 90.0, 90.0, 9.0, 10.0],
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

function onsite_U(c)
    br = BI.load_batch_result(BI.batch_path(c, 1))
    dg = BI.load_templates!(c)[1][2]
    pos = BI.grid_positions(dg, br.gidx)
    At, Bt, Ct = BI.true_cell_vectors(dg)
    far_pad = c.far_pad_steps * maximum((norm(collect(At))/dg.nx,
                                         norm(collect(Bt))/dg.ny,
                                         norm(collect(Ct))/dg.nz))
    Φ = BI.evaluate_batch_potential(br.interface, br.sigma,
            [BI.VolumeSource(copy(pos), copy(br.weights), br.densities[:, 1])], pos;
            lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], far_pad = far_pad,
            screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
    wρ = br.weights .* br.densities[:, 1]
    return dot(wρ, Φ[:, 1]) * 4π * E2 / sum(wρ)^2
end

results = Tuple{Float64,Float64}[]
for eps in EPSES
    @info "=== single-cube test: cube eps = $eps ==="
    c = load_campaign(write_toml(eps))
    prepare(c)
    solve_batch(c, 1)
    U = onsite_U(c)
    @info "single-cube done" eps U
    push!(results, (eps, U))
end

println("\n=== BULK onsite U (orbital at center; eps=10 slab on single cube) ===")
for (e, u) in results; println("  cube eps = $(e):  U = $(round(u; digits=4)) eV"); end

open(joinpath(@__DIR__, "..", "figs", "singlecube_U.tsv"), "w") do io
    println(io, "eps\tU_eV"); for (e, u) in results; println(io, "$(e)\t$(u)"); end
end
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
fig = Figure(size = (FIG_W, FIG_H))
ax = Axis(fig[1, 1]; xlabel = "substrate ε", ylabel = "U (eV)",
          title = "bulk onsite U vs substrate ε (orbital at center, ε=10 slab)")
scatterlines!(ax, [r[1] for r in results], [r[2] for r in results];
              color = QUAL.blue, marker = :circle, linewidth = LW_DATA, markersize = MS)
save(joinpath(@__DIR__, "..", "figs", "fig_singlecube_U.pdf"), fig; px_per_unit = PX_PER_UNIT)
println("wrote figs/fig_singlecube_U.pdf + figs/singlecube_U.tsv")
