#=
Derive rerun campaign TOMLs from the published ones by retargeting `name` and
`root` only.

Why not regenerate with gen_benchmark_*.jl: those rebuild the orbital list from
the .xsf templates, and any change in that list (ordering, cutoff arithmetic)
would confound the comparison we are trying to make. The geometry, orbitals,
pairing and every solver tolerance must be byte-identical to the published run
so that the ONLY difference is correct_edges and precondition, which now come
from BI's defaults rather than the TOML.

Usage:
    julia --project=codes/lattice_scale codes/rerun_2026_07/scripts/gen_rerun_campaigns.jl [TAG]

TAG defaults to "_v2". Writes codes/lattice_scale/campaigns/<name><TAG>.toml and
refuses to touch an existing ceph root.
=#

using Dates

const TAG = isempty(ARGS) ? "_v2" : ARGS[1]
const CAMPAIGNS = normpath(joinpath(@__DIR__, "..", "..", "lattice_scale", "campaigns"))
const SOURCES = ["lattice_10x10_het3x", "lattice_conv_l3"]

for base in SOURCES
    src = joinpath(CAMPAIGNS, base * ".toml")
    dst = joinpath(CAMPAIGNS, base * TAG * ".toml")
    isfile(src) || error("missing source campaign $src")

    lines = readlines(src)
    nname = nroot = 0
    for (i, ln) in enumerate(lines)
        if startswith(ln, "name = ")
            lines[i] = "name = \"$(base)$(TAG)\""; nname += 1
        elseif startswith(ln, "root = ")
            lines[i] = "root = \"/mnt/ceph/users/xgao1/four_index/$(base)$(TAG)\""; nroot += 1
        end
    end
    nname == 1 || error("expected exactly one `name =` in $src, found $nname")
    nroot == 1 || error("expected exactly one `root =` in $src, found $nroot")

    newroot = "/mnt/ceph/users/xgao1/four_index/$(base)$(TAG)"
    isdir(newroot) && error("""
        rerun root already exists: $newroot
        Refusing to overwrite. Remove it or pick another TAG.""")

    # Provenance header, so the rerun TOML is self-describing.
    hdr = [
        "# RERUN of $(base) -- generated $(Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")) by rerun_2026_07/scripts/gen_rerun_campaigns.jl",
        "# Identical to $(base).toml except `name` and `root`. The behavioural change is in",
        "# BoundaryIntegral: solve_dielectric_lattice_batch now defaults",
        "# correct_edges = true and precondition = true (design note 2026-07-31).",
        "# The published run had BOTH false.",
    ]
    open(dst, "w") do io
        for h in hdr; println(io, h); end
        for ln in lines; println(io, ln); end
    end
    println("wrote $dst")
    println("  root -> $newroot")
end
