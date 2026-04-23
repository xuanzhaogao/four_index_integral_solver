module ScreenedOrbitalSolve

using BoundaryIntegral
using Krylov
using LinearAlgebra
import BoundaryIntegral as BI

export DEFAULT_A1,
    DEFAULT_A2,
    DEFAULT_BOND_VECTOR,
    DEFAULT_PAIR_LAYOUT,
    PAPER_SHELL_LAYOUT,
    E2_4PIEPS0,
    centered_graphene_sources,
    evaluate_scattered_potential,
    evaluate_screened_pair_interaction,
    evaluate_volume_potential,
    integrate_target_potential,
    pair_specs,
    pair_targets,
    shift_volume_source,
    solve_screened_mode,
    to_eV

const DEFAULT_A1 = (2.465, 0.0, 0.0)
const DEFAULT_BOND_VECTOR = (0.0, 1.4233656765561458, 0.0)
const DEFAULT_A2 = (DEFAULT_A1[1] / 2, 1.5 * DEFAULT_BOND_VECTOR[2], 0.0)
const PAPER_SHELL_LAYOUT = (
    (pair = :U_00, orbital = :vs1, shift = (0.0, 0.0, 0.0)),
    (pair = :U_01, orbital = :vs2, shift = (0.0, 0.0, 0.0)),
    (pair = :U_02, orbital = :vs1, shift = DEFAULT_A1),
    (pair = :U_03, orbital = :vs2, shift = DEFAULT_A1),
    (pair = :U_04, orbital = :vs2, shift = DEFAULT_A2),
    (pair = :U_05, orbital = :vs1, shift = (DEFAULT_A1[1] + DEFAULT_A2[1], DEFAULT_A1[2] + DEFAULT_A2[2], 0.0)),
)
const DEFAULT_PAIR_LAYOUT = PAPER_SHELL_LAYOUT
const E2_4PIEPS0 = 14.3996
const DEFAULT_ORBITAL_1 = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "density_data", "graphene_00001_5x5x1_shifted.xsf"))
const DEFAULT_ORBITAL_2 = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "density_data", "graphene_00002_5x5x1_shifted.xsf"))

function _load_squared_datagrid(path::AbstractString)
    _, raw_datagrid = BI.read_xsf(path)
    squared_values = copy(raw_datagrid.values)
    squared_values .*= squared_values
    return merge(raw_datagrid, (; values = squared_values))
end

function _shift_datagrid(datagrid, shift::NTuple{3, <:Real})
    origin = ntuple(i -> Float64(datagrid.origin[i]) + Float64(shift[i]), 3)
    return merge(datagrid, (; origin = origin))
end

function _density_centroid(vs)
    weights = vs.weights .* vs.density
    total_weight = sum(weights)
    iszero(total_weight) && throw(ArgumentError("cannot compute centroid for zero total density"))
    return tuple(
        (
            sum(vs.positions[dim, :] .* weights) / total_weight
            for dim in 1:size(vs.positions, 1)
        )...,
    )
end

function _centering_shift(datagrid; tol::Real = 1e-4)
    vs = BI.VolumeSource(datagrid, tol = tol)
    centroid = _density_centroid(vs)
    return ntuple(i -> -Float64(centroid[i]), 3)
end

function _centered_datagrid(path::AbstractString; tol::Real = 1e-4)
    datagrid = _load_squared_datagrid(path)
    return _shift_datagrid(datagrid, _centering_shift(datagrid; tol = tol))
end

function centered_graphene_sources(;
    orbital_1::AbstractString = DEFAULT_ORBITAL_1,
    orbital_2::AbstractString = DEFAULT_ORBITAL_2,
    tol::Real = 1e-3,
    z_center::Real = 0.0,
)
    datagrid_1_raw = _load_squared_datagrid(orbital_1)
    datagrid_2_raw = _load_squared_datagrid(orbital_2)
    centering_shift = _centering_shift(datagrid_1_raw; tol = tol)
    shared_shift = (
        centering_shift[1],
        centering_shift[2],
        centering_shift[3] + Float64(z_center),
    )
    datagrid_1 = _shift_datagrid(datagrid_1_raw, shared_shift)
    datagrid_2 = _shift_datagrid(datagrid_2_raw, shared_shift)
    vs1 = BI.VolumeSource(datagrid_1, tol = tol)
    vs2 = BI.VolumeSource(datagrid_2, tol = tol)
    return (
        datagrid_1 = datagrid_1,
        datagrid_2 = datagrid_2,
        vs1 = vs1,
        vs2 = vs2,
        shared_shift = shared_shift,
    )
end

function shift_volume_source(vs::BI.VolumeSource{T, 3}, shift::NTuple{3, <:Real}) where {T}
    shifted_positions = copy(vs.positions)
    for dim in 1:3
        shifted_positions[dim, :] .+= Float64(shift[dim])
    end
    return BI.VolumeSource(shifted_positions, copy(vs.weights), copy(vs.density))
end

function _pair_target(base_sources, layout)
    named_pairs = [
        Pair(spec.pair, shift_volume_source(getproperty(base_sources, spec.orbital), spec.shift))
        for spec in layout
    ]
    return (; named_pairs...)
end

function _normalization(vs::BI.VolumeSource)
    return sum(vs.weights .* vs.density)
end

function pair_targets(
    vs1::BI.VolumeSource{T, 3},
    vs2::BI.VolumeSource{T, 3};
    layout = DEFAULT_PAIR_LAYOUT,
) where {T}
    return _pair_target((; vs1, vs2), layout)
end

function pair_specs(
    vs1::BI.VolumeSource{T, 3},
    vs2::BI.VolumeSource{T, 3};
    layout = DEFAULT_PAIR_LAYOUT,
) where {T}
    targets = pair_targets(vs1, vs2; layout = layout)
    n1 = _normalization(vs1)
    n2 = _normalization(vs2)
    norms = (; vs1 = n1, vs2 = n2)
    return [
        (
            pair = String(spec.pair),
            target_vs = getproperty(targets, spec.pair),
            Na = n1,
            Nb = getproperty(norms, spec.orbital),
        )
        for spec in layout
    ]
end

function integrate_target_potential(vs::BI.VolumeSource, potential::AbstractVector{<:Real})
    length(potential) == length(vs.density) || throw(ArgumentError("potential length must match target density length"))
    return dot(vs.weights .* vs.density, potential)
end

function _target_matrix(targets::AbstractVector{<:NTuple{3, <:Real}})
    mat = Matrix{Float64}(undef, 3, length(targets))
    for i in eachindex(targets)
        mat[1, i] = Float64(targets[i][1])
        mat[2, i] = Float64(targets[i][2])
        mat[3, i] = Float64(targets[i][3])
    end
    return mat
end

_target_matrix(targets::AbstractMatrix{<:Real}) = Matrix{Float64}(targets)

function evaluate_volume_potential(
    source_vs::BI.VolumeSource{Float64, 3},
    targets;
    tol::Real = 1e-6,
    kmax = nothing,
)
    target_matrix = _target_matrix(targets)
    size(target_matrix, 1) == 3 || throw(ArgumentError("targets must have shape (3, n)"))
    charges = source_vs.weights .* source_vs.density
    resolved_kmax = isnothing(kmax) ? BI._estimate_tkm3dc_kmax(source_vs) : Float64(kmax)
    values = BoundaryIntegral.TKM3D.ltkm3dc(
        Float64(tol),
        source_vs.positions;
        charges = charges,
        targets = target_matrix,
        pgt = 1,
        kmax = resolved_kmax,
    )
    values.ier == 0 || error("TKM3D.ltkm3dc failed with ier=$(values.ier)")
    return real.(values.pottarg)
end

function evaluate_scattered_potential(
    interface,
    sigma::AbstractVector{Float64},
    targets;
    fmm_tol::Real = 1e-6,
    hcubature_atol::Real = 1e-6,
    range_factor::Real = 5.0,
    include_edges_src::Bool = false,
)
    target_matrix = _target_matrix(targets)
    op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
        interface,
        target_matrix,
        Float64(fmm_tol),
        Float64(hcubature_atol),
        Float64(range_factor);
        include_edges_src = include_edges_src,
    )
    return op * sigma
end

to_eV(raw, Na, Nb) = raw * 4π * E2_4PIEPS0 / (Na * Nb)

function evaluate_screened_pair_interaction(
    source_vs::BI.VolumeSource{Float64, 3},
    target_vs::BI.VolumeSource{Float64, 3},
    interface,
    sigma::AbstractVector{Float64};
    volume_tol::Real = 1e-6,
    scatter_fmm_tol::Real = 1e-6,
    scatter_atol::Real = 1e-6,
    direct_kmax = nothing,
    include_edges_src::Bool = false,
    volume_evaluator::Function = evaluate_volume_potential,
    scatter_evaluator::Function = evaluate_scattered_potential,
)
    target_positions = target_vs.positions
    target_weights = target_vs.weights .* target_vs.density
    u_int = volume_evaluator(source_vs, target_positions; tol = volume_tol, kmax = direct_kmax)
    u_scatter = scatter_evaluator(
        interface,
        sigma,
        target_positions;
        fmm_tol = scatter_fmm_tol,
        hcubature_atol = scatter_atol,
        include_edges_src = include_edges_src,
    )

    u_int_raw = dot(u_int, target_weights)
    u_scatter_raw = dot(u_scatter, target_weights)

    return (
        u_int_raw = u_int_raw,
        u_scatter_raw = u_scatter_raw,
        u_total_raw = u_int_raw + u_scatter_raw,
    )
end

function solve_screened_mode(
    source_vs::BI.VolumeSource{Float64, 3},
    target_specs,
    Lx::Real,
    Ly::Real,
    Lz::Real,
    eps_in::Real,
    eps_out::Real,
    mode;
    n_quad::Integer = 6,
    edge_refine_level::Integer = 4,
    rhs_tol::Real = 1e-4,
    lhs_tol::Real = 1e-6,
    gmres_atol::Real = 1e-6,
    gmres_rtol::Real = 1e-6,
    max_order::Integer = 128,
    max_depth::Integer = 100,
    include_edges_src::Bool = false,
    include_edges_trg::Bool = false,
    volume_tol::Real = 1e-6,
    scatter_range_factor::Real = 5.0,
    tkm_kmax = nothing,
)
    l_ec = Float64(Lz) / 2.0^edge_refine_level * 1.01
    screened_vs = BI.screened_volume_source(Float64(Lx), Float64(Ly), Float64(Lz), source_vs, Float64(eps_in), Float64(eps_out), mode; tol = Float64(rhs_tol))
    resolved_kmax = isnothing(tkm_kmax) ? BI._estimate_tkm3dc_kmax(screened_vs) : Float64(tkm_kmax)

    interface = BI.single_dielectric_box3d_rhs_adaptive(
        Float64(Lx),
        Float64(Ly),
        Float64(Lz),
        Int(n_quad),
        screened_vs,
        1.0,
        l_ec,
        Float64(rhs_tol),
        Float64(eps_in),
        Float64(eps_out),
        Float64;
        max_depth = Int(max_depth),
        tkm_kmax = resolved_kmax,
    )

    rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, screened_vs, 1.0, Float64(rhs_tol))
    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(
        interface,
        Float64(lhs_tol),
        Float64(lhs_tol),
        Int(max_order);
        include_edges_src = include_edges_src,
        include_edges_trg = include_edges_trg,
    )
    sigma, status = Krylov.gmres(lhs, rhs, atol = Float64(gmres_atol), rtol = Float64(gmres_rtol), verbose = 0)
    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))

    pair_results = NamedTuple[]
    for spec in target_specs
        raw = evaluate_screened_pair_interaction(
            screened_vs,
            spec.target_vs,
            interface,
            sigma;
            volume_tol = volume_tol,
            scatter_fmm_tol = lhs_tol,
            scatter_atol = lhs_tol,
            direct_kmax = resolved_kmax,
            include_edges_src = include_edges_src,
        )
        push!(pair_results, merge(spec, raw, (
            u_int_ev = to_eV(raw.u_int_raw, spec.Na, spec.Nb),
            u_scatter_ev = to_eV(raw.u_scatter_raw, spec.Na, spec.Nb),
            u_total_ev = to_eV(raw.u_total_raw, spec.Na, spec.Nb),
        )))
    end

    return (
        interface = interface,
        screened_vs = screened_vs,
        sigma = sigma,
        gmres_status = status,
        residual = residual,
        tkm_kmax = resolved_kmax,
        n_interface_points = BI.num_points(interface),
        pair_results = pair_results,
    )
end

end
