using HDF5
using CSV
using DataFrames
using Printf

const H5 = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/bilayer/k_161601_nb_288_c_15/crpa.mbpt_clean.h5"
const OUT_CSV = joinpath(@__DIR__, "..", "data", "coqui_bare_vloc_abcd.csv")
const OURS_CSV = joinpath(@__DIR__, "..", "data", "bare_bilayer_graphene.csv")

# Read Vloc_abcd. HDF5.jl returns the array with axis order reversed relative
# to the file's row-major (C) layout, so the on-disk shape {4, 4, 4, 4, 2}
# (with last axis = re/im) becomes Julia size (2, 4, 4, 4, 4).
function load_vloc()
    raw = h5read(H5, "/scf/iter0/downfolded_model/Vloc_abcd")
    @assert size(raw) == (2, 4, 4, 4, 4) "unexpected Vloc shape: $(size(raw))"
    re = @view raw[1, :, :, :, :]
    im = @view raw[2, :, :, :, :]
    return (re = Array(re), im = Array(im))
end

# Index convention check: CoQui-style four-index U_{abcd} = (ab|cd) with
# physics chemist convention U_{abcd} = ∫∫ φ_a*(r) φ_b(r) (1/|r-r'|) φ_c*(r')
# φ_d(r') dr dr'. For real orbitals:
#   onsite U_i = V_{i i i i}
#   density-density V_{ij} = V_{i i j j}
#   Hund spin-flip  J^sf_{ij} = V_{i j j i}
#   Hund pair-hop   J^ph_{ij} = V_{i j i j}
# We verify the channel-average matches the text-output number to confirm
# that this ordering matches what CoQui's "orbital-average" summary uses.

function channel_averages(V)
    n = 4
    onsite = [V[i, i, i, i] for i in 1:n]
    intra_avg = sum(onsite) / n

    inter_pairs = [(i, j) for i in 1:n for j in 1:n if i != j]
    inter_dd = [V[i, i, j, j] for (i, j) in inter_pairs]
    inter_avg = sum(inter_dd) / length(inter_dd)

    hund_sf = [V[i, j, j, i] for (i, j) in inter_pairs]
    hund_ph = [V[i, j, i, j] for (i, j) in inter_pairs]
    hund_sf_avg = sum(hund_sf) / length(hund_sf)
    hund_ph_avg = sum(hund_ph) / length(hund_ph)

    return (
        intra = intra_avg,
        inter = inter_avg,
        hund_sf = hund_sf_avg,
        hund_ph = hund_ph_avg,
    )
end

# Bare tensor in CoQui's `crpa.mbpt_clean.h5` is stored in Hartree (atomic units).
const HARTREE_TO_EV = 27.211386245988

function unique_pair_rows(V)
    rows = NamedTuple[]
    n = 4
    for i in 1:n
        push!(rows, (
            i = i, j = i, channel = "onsite",
            u_ev_coqui = V[i, i, i, i],
        ))
    end
    for i in 1:n, j in (i + 1):n
        push!(rows, (i = i, j = j, channel = "nn",      u_ev_coqui = V[i, i, j, j]))
        push!(rows, (i = i, j = j, channel = "hund_sf", u_ev_coqui = V[i, j, j, i]))
        push!(rows, (i = i, j = j, channel = "hund_ph", u_ev_coqui = V[i, j, i, j]))
    end
    return rows
end

function main()
    println("Reading $H5")
    vloc = load_vloc()
    Vre = vloc.re

    println("Max |Im(V)| (Hartree) = ", maximum(abs.(vloc.im)))

    # First check whether values look like eV or Hartree by averaging the
    # diagonal and comparing to the text-output 17.669390 eV.
    raw_intra = sum(Vre[i, i, i, i] for i in 1:4) / 4
    @printf("Raw intra-orbital average (h5 units) = %.6f\n", raw_intra)
    @printf("  → if Hartree, in eV = %.6f\n", raw_intra * HARTREE_TO_EV)
    @printf("  → if eV, in eV      = %.6f\n", raw_intra)
    println("  CoQui text output    = 17.669390 eV")

    use_hartree = abs(raw_intra * HARTREE_TO_EV - 17.669390) < abs(raw_intra - 17.669390)
    scale = use_hartree ? HARTREE_TO_EV : 1.0
    println(use_hartree ? "Detected Hartree storage; converting to eV." : "Detected eV storage; no conversion.")

    V = Vre .* scale

    avgs = channel_averages(V)
    println()
    println("Recomputed orbital averages from h5 (eV):")
    @printf("  intra   = %10.6f   (text: 17.669390)\n", avgs.intra)
    @printf("  inter   = %10.6f   (text:  5.785749)\n", avgs.inter)
    @printf("  hund_sf = %10.6f   (text:  0.041867)\n", avgs.hund_sf)
    @printf("  hund_ph = %10.6f   (text:  0.041867)\n", avgs.hund_ph)

    rows = unique_pair_rows(V)
    coqui_df = DataFrame(rows)
    CSV.write(OUT_CSV, coqui_df)
    println()
    println("Wrote $(nrow(coqui_df)) reference rows to $OUT_CSV")
    show(coqui_df; allrows = true, allcols = true)
    println()

    # Merge with ours
    ours = CSV.read(OURS_CSV, DataFrame)
    cmp = innerjoin(
        select(ours, [:i, :j, :channel, :u_ev, :pair_class]),
        coqui_df,
        on = [:i, :j, :channel],
    )
    cmp.diff_ev  = cmp.u_ev .- cmp.u_ev_coqui
    cmp.rel_err  = 100 .* cmp.diff_ev ./ cmp.u_ev_coqui
    sort!(cmp, [:channel, :i, :j])
    println()
    println("Per-pair comparison (ours vs CoQui Vloc_abcd):")
    show(cmp; allrows = true, allcols = true)
    println()
    CSV.write(joinpath(@__DIR__, "..", "data", "bare_bilayer_vs_coqui.csv"), cmp)
end

main()
