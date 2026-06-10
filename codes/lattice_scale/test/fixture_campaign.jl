# test/fixture_campaign.jl — a tiny self-contained campaign for pipeline tests:
# 6x6x6 template (3 cells/axis, 2 steps/cell), Gaussian blob, 2x1 lattice, 1 sublattice.
function write_fixture_xsf(path::AbstractString)
    open(path, "w") do io
        println(io, "CRYSTAL\nPRIMVEC\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0")
        println(io, "PRIMCOORD\n1 1\nX 0.75 0.75 0.75")
        println(io, "BEGIN_BLOCK_DATAGRID_3D\n g\nBEGIN_DATAGRID_3D_g")
        println(io, "6 6 6\n0.0 0.0 0.0")
        println(io, "2.5 0.0 0.0\n0.0 2.5 0.0\n0.0 0.0 2.5")
        vals = Float64[]
        for k in 1:6, j in 1:6, i in 1:6
            x = ((i-1)/6)*3.0; y = ((j-1)/6)*3.0; z = ((k-1)/6)*3.0
            push!(vals, exp(-((x-0.75)^2 + (y-0.75)^2 + (z-0.75)^2) / (2*0.4^2)))
        end
        for chunk in Iterators.partition(vals, 6)
            println(io, join(string.(chunk), " "))
        end
        println(io, "END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D")
    end
end

function write_fixture_campaign(dir::AbstractString)
    xsf = joinpath(dir, "orb.xsf")
    write_fixture_xsf(xsf)
    toml = joinpath(dir, "campaign.toml")
    write(toml, """
    [campaign]
    name = "mini"
    root = "$(joinpath(dir, "out"))"

    [orbitals]
    xsf = ["$xsf"]

    [lattice]
    nx = 2
    ny = 1
    neighbor_cutoff = 1.2

    [dielectrics]
    eps_out = 1.0
    boxes = [[1.5, 1.5, 1.5, 8.0, 8.0, 8.0, 2.0]]

    [solve]
    n_quad = 4
    edge_refine_level = 1
    rhs_tol = 1e-2
    lhs_tol = 1e-6
    gmres_rtol = 1e-8
    support_rtol = 1e-6
    volume_tol = 1e-8
    max_order = 8
    max_depth = 16

    [batching]
    n_centers_per_batch = 1

    [eval]
    far_pad_steps = 2.0
    """)
    return toml
end
