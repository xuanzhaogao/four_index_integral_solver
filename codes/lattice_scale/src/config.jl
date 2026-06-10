"""
    Campaign

Parsed campaign.toml + derived paths. Heavy template grids are NOT loaded here —
see `load_templates!`.
"""
struct Campaign
    name::String
    root::String                       # ceph output directory
    xsf::Vector{String}                # one template per sublattice
    nx::Int
    ny::Int
    neighbor_cutoff::Float64
    eps_out::Float64
    boxes::Vector{BoundaryIntegral.BoxGeom}
    epses::Vector{Float64}
    solve::Dict{String,Float64}        # n_quad, edge_refine_level, rhs_tol, lhs_tol,
                                       # gmres_rtol, support_rtol, volume_tol,
                                       # max_order, max_depth (numeric; Int-coerced on use)
    n_centers_per_batch::Int
    far_pad_steps::Float64
end

function load_campaign(toml_path::AbstractString)
    d = TOML.parsefile(toml_path)
    c, l, di, s = d["campaign"], d["lattice"], d["dielectrics"], d["solve"]
    haskey(d, "batching") || error("campaign.toml is missing the [batching] section")
    boxes = BoundaryIntegral.BoxGeom[]
    epses = Float64[]
    for row in di["boxes"]
        length(row) == 7 || error("dielectrics.boxes rows are [cx cy cz Lx Ly Lz eps]")
        push!(boxes, (center = (Float64(row[1]), Float64(row[2]), Float64(row[3])),
                      Lx = Float64(row[4]), Ly = Float64(row[5]), Lz = Float64(row[6])))
        push!(epses, row[7])
    end
    solve = Dict{String,Float64}(k => Float64(v) for (k, v) in s)
    return Campaign(c["name"], c["root"], String.(d["orbitals"]["xsf"]),
        Int(l["nx"]), Int(l["ny"]), Float64(l["neighbor_cutoff"]),
        Float64(get(di, "eps_out", 1.0)), boxes, epses, solve,
        Int(d["batching"]["n_centers_per_batch"]),
        Float64(get(get(d, "eval", Dict()), "far_pad_steps", 2.0)))
end

manifest_path(c::Campaign)  = joinpath(c.root, "manifest.tsv")
centers_path(c::Campaign)   = joinpath(c.root, "centers.tsv")
targets_path(c::Campaign)   = joinpath(c.root, "targets.jls")
rho_store_path(c::Campaign) = joinpath(c.root, "rho_store.jls")
logs_dir(c::Campaign)       = joinpath(c.root, "logs")
batch_path(c::Campaign, id::Int) = joinpath(c.root, "batches", @sprintf("batch_%04d.jls", id))
v_path(c::Campaign, id::Int)     = joinpath(c.root, "V", @sprintf("V_%04d.jls", id))

# NOTE: named campaign_l_ec, NOT resolved_l_ec — BoundaryIntegral exports
# resolved_l_ec(::SystemInput) and a same-named definition here would clash.
campaign_l_ec(c::Campaign) =
    minimum(b.Lz for b in c.boxes) / 2.0^Int(c.solve["edge_refine_level"]) * 1.01

# worker-local template cache: path => (structure, datagrid)
const TEMPLATE_CACHE = Dict{String,Any}()
function load_templates!(c::Campaign)
    return [get!(TEMPLATE_CACHE, p) do
                BoundaryIntegral.read_xsf(p)
            end for p in c.xsf]
end
