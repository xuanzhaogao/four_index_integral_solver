# Cross-check (spec §8.1): the campaign's WITHIN-BATCH V entries for the batch anchored
# at center 1 must match the existing .bie path (four_index_integrals with circshift
# LATTICE images) to solver tolerance. Run on a compute node AFTER demo_2x2 finished
# solve+consolidate+eval:   julia --project scripts/compare_anchor.jl
#
# CAVEAT: the .bie LATTICE path wraps on the 5x5 template grid; the campaign translates.
# Within a 2x2 flake all shifts are <= 1 cell with the orbital mid-grid, so wrap effects
# sit below support_rtol = 1e-4. If the comparison fails, FIRST check whether the .bie
# group/centroid conventions place center 1's partners identically (print both center lists).
using CampaignLib, BoundaryIntegral, Printf

function main()
    c = load_campaign(joinpath(@__DIR__, "..", "campaigns", "demo_2x2.toml"))
    centers = read_centers(centers_path(c))
    spec = first(read_manifest(manifest_path(c)))            # batch anchored at center 1
    vr = CampaignLib.load_v_rows(v_path(c, spec.batch_id))

    # build the equivalent .bie (center 1 + its batch partners as LATTICE images of the
    # SAME templates) in a temp dir, run the reference path
    byid = Dict(ct.id => ct for ct in centers)
    partners = sort(unique(reduce(vcat, [[p[1], p[2]] for p in spec.pairs])))
    bie = tempname() * ".bie"
    open(bie, "w") do io
        println(io, "UNITS bohr\n\nBEGIN_DIELECTRICS\nEPS_OUT 1.0")
        b = c.boxes[1]
        @printf(io, "  %.3f %.3f %.3f    %.3f %.3f %.3f    %.3f\n",
            b.center..., b.Lx, b.Ly, b.Lz, c.epses[1])
        println(io, "END_DIELECTRICS\n\nBEGIN_ORBITALS")
        for id in partners
            ct = byid[id]
            println(io, "  $id   $(c.xsf[ct.template_id])   LATTICE $(ct.Rx) $(ct.Ry) 0")
        end
        println(io, "END_ORBITALS\n\nBEGIN_GROUPING")
        println(io, "  1 : $(join([p[2] for p in spec.pairs], ' '))")
        println(io, "END_GROUPING\n\nBEGIN_SOLVE")
        for (k, v) in [("N_QUAD", 6), ("EDGE_REFINE_LEVEL", 2), ("RHS_TOL", 1e-3),
                       ("LHS_TOL", 1e-5), ("GMRES_RTOL", 1e-5), ("SUPPORT_RTOL", 1e-4),
                       ("VOLUME_TOL", 1e-5), ("MAX_ORDER", 8), ("MAX_DEPTH", 128)]
            println(io, "  $k $v")
        end
        println(io, "END_SOLVE")
    end
    ref = four_index_integrals(bie, 1)

    # compare the shared entries (campaign rows for this batch's own pairs)
    rowof = Dict(p => i for (i, p) in enumerate(vr.target_pairs))
    worst = 0.0
    for (a, pa) in enumerate(spec.pairs)
        for (b2, pb) in enumerate(spec.pairs)
            v_c = vr.V[rowof[pa], b2]
            v_r = ref.V[a, b2]
            rel = abs(v_c - v_r) / max(abs(v_r), 1e-300)
            worst = max(worst, rel)
            @printf("%-14s %-14s  campaign % .6e   ref % .6e   rel %.2e\n",
                "$(pa)", "$(pb)", v_c, v_r, rel)
        end
    end
    @printf("\nworst relative difference: %.3e  (expect <~ 1e-2 at these tolerances)\n", worst)
end

main()
