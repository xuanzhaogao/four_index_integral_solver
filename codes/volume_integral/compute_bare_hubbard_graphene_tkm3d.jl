using BoundaryIntegral
using CSV
using DataFrames
using LinearAlgebra
using TKM3D

import BoundaryIntegral as BI

include(joinpath(@__DIR__, "tkm3d_hubbard_utils.jl"))

const DEFAULT_TKM_TOLERANCES = [1.0e-2, 1.0e-3, 1.0e-4]
const DEFAULT_OUTPUT_FILE = joinpath(@__DIR__, "data", "hubbard_graphene_tkm3d.csv")
const GRAPHENE_1 = joinpath(@__DIR__, "../../density_data/graphene_00001_5x5x1_shifted.xsf")
const GRAPHENE_2 = joinpath(@__DIR__, "../../density_data/graphene_00002_5x5x1_shifted.xsf")
const Z_SHIFT = -7.920155482424242
const A1_PRIM = (2.465, 0.0, 0.0)
const E2_4PIEPS0 = 14.3996

function to_eV(raw, Na, Nb)
    return raw * 4π * E2_4PIEPS0 / (Na * Nb)
end

function parse_tolerances(args)
    return isempty(args) ? normalize_tolerances(DEFAULT_TKM_TOLERANCES) : normalize_tolerances(parse.(Float64, args))
end

function load_graphene_volume_sources()
    isfile(GRAPHENE_1) || error("Missing XSF input: $(GRAPHENE_1)")
    isfile(GRAPHENE_2) || error("Missing XSF input: $(GRAPHENE_2)")

    _, datagrid_1 = BI.read_xsf(GRAPHENE_1)
    datagrid_1.values .*= datagrid_1.values

    _, datagrid_2 = BI.read_xsf(GRAPHENE_2)
    datagrid_2.values .*= datagrid_2.values

    vs1 = BoundaryIntegral.VolumeSource(datagrid_1, shift = (0.0, 0.0, Z_SHIFT), tol = 1.0e-6)
    vs2 = BoundaryIntegral.VolumeSource(datagrid_2, shift = (0.0, 0.0, Z_SHIFT), tol = 1.0e-6)
    vs1_a1 = BoundaryIntegral.VolumeSource(datagrid_1, shift = (A1_PRIM[1], A1_PRIM[2], Z_SHIFT), tol = 1.0e-6)
    vs2_a1 = BoundaryIntegral.VolumeSource(datagrid_2, shift = (A1_PRIM[1], A1_PRIM[2], Z_SHIFT), tol = 1.0e-6)

    return (
        vs1 = vs1,
        vs2 = vs2,
        vs1_a1 = vs1_a1,
        vs2_a1 = vs2_a1,
    )
end

function geometry_kmax_metadata(source_positions)
    _, spacings, _, Δk, axis_nyquist = TKM3D._ltkm3dc_source_nyquist_geometry(source_positions)
    return (
        kmax = Float64(TKM3D._ltkm3dc_estimate_kmax(source_positions)),
        spacing_x = spacings[1],
        spacing_y = spacings[2],
        spacing_z = spacings[3],
        delta_k_x = Δk[1],
        delta_k_y = Δk[2],
        delta_k_z = Δk[3],
        axis_nyquist_x = axis_nyquist[1],
        axis_nyquist_y = axis_nyquist[2],
        axis_nyquist_z = axis_nyquist[3],
    )
end

function evaluate_tkm_pair(source_positions, source_charges, target_vs, tol, kmax_metadata)
    values = TKM3D.ltkm3dc(
        tol,
        source_positions;
        charges = source_charges,
        targets = target_vs.positions,
        pgt = 1,
        kmax = kmax_metadata.kmax,
    )
    values.ier == 0 || error("TKM3D.ltkm3dc failed with ier=$(values.ier)")

    target_weights = target_vs.density .* target_vs.weights
    potential = real.(values.pottarg)

    return (
        U_raw = sum(potential .* target_weights),
        kmax = kmax_metadata.kmax,
        spacing_x = kmax_metadata.spacing_x,
        spacing_y = kmax_metadata.spacing_y,
        spacing_z = kmax_metadata.spacing_z,
        delta_k_x = kmax_metadata.delta_k_x,
        delta_k_y = kmax_metadata.delta_k_y,
        delta_k_z = kmax_metadata.delta_k_z,
        axis_nyquist_x = kmax_metadata.axis_nyquist_x,
        axis_nyquist_y = kmax_metadata.axis_nyquist_y,
        axis_nyquist_z = kmax_metadata.axis_nyquist_z,
    )
end

function main(args = ARGS)
    tolerances = parse_tolerances(args)
    sources = load_graphene_volume_sources()

    N1 = sum(sources.vs1.density .* sources.vs1.weights)
    N2 = sum(sources.vs2.density .* sources.vs2.weights)
    source_charges = sources.vs1.density .* sources.vs1.weights
    pair_specs = [
        (pair = "U_00", target_vs = sources.vs1, Na = N1, Nb = N1),
        (pair = "U_01", target_vs = sources.vs2, Na = N1, Nb = N2),
        (pair = "U_02", target_vs = sources.vs1_a1, Na = N1, Nb = N1),
        (pair = "U_03", target_vs = sources.vs2_a1, Na = N1, Nb = N2),
    ]
    pair_samples = Dict(spec.pair => NamedTuple[] for spec in pair_specs)
    kmax_metadata = geometry_kmax_metadata(sources.vs1.positions)

    println("Using fixed geometry-derived kmax=$(kmax_metadata.kmax)")

    for tol in tolerances
        println("Running ltkm3dc for tol=$(tol)")

        for spec in pair_specs
            println("  Computing $(spec.pair)")
            values = evaluate_tkm_pair(sources.vs1.positions, source_charges, spec.target_vs, tol, kmax_metadata)
            push!(pair_samples[spec.pair], (
                tol = tol,
                U_raw = values.U_raw,
                U_ev = to_eV(values.U_raw, spec.Na, spec.Nb),
                kmax = values.kmax,
                spacing_x = values.spacing_x,
                spacing_y = values.spacing_y,
                spacing_z = values.spacing_z,
                delta_k_x = values.delta_k_x,
                delta_k_y = values.delta_k_y,
                delta_k_z = values.delta_k_z,
                axis_nyquist_x = values.axis_nyquist_x,
                axis_nyquist_y = values.axis_nyquist_y,
                axis_nyquist_z = values.axis_nyquist_z,
            ))
        end
    end

    rows = vcat([annotate_reference_errors(spec.pair, pair_samples[spec.pair]) for spec in pair_specs]...)

    table = DataFrame(rows)
    CSV.write(DEFAULT_OUTPUT_FILE, table)

    println("Wrote $(nrow(table)) rows to $(DEFAULT_OUTPUT_FILE)")
    show(table; allrows = true, allcols = true)
    println()
end

main()
