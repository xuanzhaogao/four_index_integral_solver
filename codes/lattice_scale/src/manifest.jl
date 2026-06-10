"""
    CenterInfo

One orbital site of the flake: lattice cell (Rx, Ry), sublattice template, the
global-frame integer offset `steps`, and the real-space center (centroid + R).
"""
struct CenterInfo
    id::Int
    template_id::Int
    Rx::Int
    Ry::Int
    steps::NTuple{3,Int}
    center::NTuple{3,Float64}
end

struct BatchSpec
    batch_id::Int
    anchors::Vector{Int}                 # center ids whose pair lists were merged
    pairs::Vector{Tuple{Int,Int}}
end

Base.:(==)(a::CenterInfo, b::CenterInfo) =
    a.id == b.id && a.template_id == b.template_id && a.Rx == b.Rx && a.Ry == b.Ry &&
    a.steps == b.steps && a.center == b.center
Base.:(==)(a::BatchSpec, b::BatchSpec) =
    a.batch_id == b.batch_id && a.anchors == b.anchors && a.pairs == b.pairs

"""
    enumerate_centers(nx, ny, primvec, centroids, steps_per_cell) -> Vector{CenterInfo}

Pure geometry (no file IO): center id = (Ry*nx + Rx)*n_sub + template_id, position =
template centroid + Rx·a1 + Ry·a2, steps = Rx·steps(a1) + Ry·steps(a2).
`steps_per_cell` = the per-template-grid steps of (a1, a2) from `lattice_grid_steps`.
`primvec` rows are the lattice vectors (a1 = primvec[1,:], a2 = primvec[2,:]).
"""
function enumerate_centers(nx::Int, ny::Int, primvec::AbstractMatrix,
        centroids::Vector{NTuple{3,Float64}}, steps_per_cell::NTuple{2,NTuple{3,Int}})
    nsub = length(centroids)
    a1 = primvec[1, :]; a2 = primvec[2, :]
    s1, s2 = steps_per_cell
    out = CenterInfo[]
    for Ry in 0:(ny-1), Rx in 0:(nx-1), t in 1:nsub
        id = (Ry * nx + Rx) * nsub + t
        steps = ntuple(d -> Rx * s1[d] + Ry * s2[d], 3)
        ctr = ntuple(d -> centroids[t][d] + Rx * a1[d] + Ry * a2[d], 3)
        push!(out, CenterInfo(id, t, Rx, Ry, steps, ctr))
    end
    return out
end

"Unique pairs (i ≤ j) with center distance ≤ cutoff. On-site pairs (i,i) included."
function enumerate_pairs(centers::Vector{CenterInfo}, cutoff::Real)
    byid = sort(centers; by = c -> c.id)
    pairs = Tuple{Int,Int}[]
    for (m, ci) in enumerate(byid)
        for cj in byid[m:end]
            d = sqrt(sum(abs2, ci.center .- cj.center))
            d <= cutoff && push!(pairs, (ci.id, cj.id))
        end
    end
    return pairs
end

"""
    build_batches(pairs, n_centers_per_batch) -> Vector{BatchSpec}

Each pair belongs to its anchor (= min id); consecutive anchors are merged
`n_centers_per_batch` at a time. Every pair lands in exactly one batch.
"""
function build_batches(pairs::Vector{Tuple{Int,Int}}, n_centers_per_batch::Int)
    by_anchor = Dict{Int,Vector{Tuple{Int,Int}}}()
    for p in pairs
        push!(get!(by_anchor, min(p[1], p[2]), Tuple{Int,Int}[]), p)
    end
    anchors = sort(collect(keys(by_anchor)))
    out = BatchSpec[]
    bid = 0
    for grp in Iterators.partition(anchors, n_centers_per_batch)
        bid += 1
        ps = reduce(vcat, (sort(by_anchor[a]) for a in grp))
        push!(out, BatchSpec(bid, collect(grp), ps))
    end
    return out
end

# ---- TSV IO (human-greppable; deterministic) ----
function write_centers(path::AbstractString, centers::Vector{CenterInfo})
    d = dirname(path); isempty(d) || mkpath(d)
    open(path, "w") do io
        println(io, "id\ttemplate\tRx\tRy\tsx\tsy\tsz\tcx\tcy\tcz")
        for c in sort(centers; by = c -> c.id)
            # repr(x::Float64) always includes a decimal point (e.g. "1.0"), which
            # distinguishes float columns from int columns on read-back.
            flt = join(repr.(c.center), '\t')
            int = join([c.id, c.template_id, c.Rx, c.Ry, c.steps...], '\t')
            println(io, int, '\t', flt)
        end
    end
end

function read_centers(path::AbstractString)
    out = CenterInfo[]
    for (n, line) in enumerate(eachline(path))
        n == 1 && continue
        f = split(line, '\t')
        push!(out, CenterInfo(parse(Int, f[1]), parse(Int, f[2]),
            parse(Int, f[3]), parse(Int, f[4]),
            (parse(Int, f[5]), parse(Int, f[6]), parse(Int, f[7])),
            (parse(Float64, f[8]), parse(Float64, f[9]), parse(Float64, f[10]))))
    end
    return out
end

function write_manifest(path::AbstractString, batches::Vector{BatchSpec})
    d = dirname(path); isempty(d) || mkpath(d)
    open(path, "w") do io
        println(io, "batch_id\tanchors\tK\tpairs")
        for b in sort(batches; by = b -> b.batch_id)
            ps = join(("$(i):$(j)" for (i, j) in b.pairs), ';')
            println(io, join([b.batch_id, join(b.anchors, ','), length(b.pairs), ps], '\t'))
        end
    end
end

function read_manifest(path::AbstractString)
    out = BatchSpec[]
    for (n, line) in enumerate(eachline(path))
        n == 1 && continue
        f = split(line, '\t')
        pairs = isempty(strip(f[4])) ? Tuple{Int,Int}[] :
            [(parse(Int, split(p, ':')[1]), parse(Int, split(p, ':')[2]))
             for p in split(f[4], ';')]
        push!(out, BatchSpec(parse(Int, f[1]), parse.(Int, split(f[2], ',')), pairs))
    end
    return out
end
