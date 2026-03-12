using LinearAlgebra
using Printf

using BoundaryIntegral
import BoundaryIntegral as BI

const ORBITAL_FILE = joinpath(@__DIR__, "../../../density_data/graphene_00001_5x5x1_shifted.xsf")
const TOL = 1.0e-4
const NUFFT_EPS = 1.0e-12

function _centered_mode_indices(n::Int)
    n >= 1 || throw(ArgumentError("mode count must be positive"))
    if iseven(n)
        return collect((-n ÷ 2):(n ÷ 2 - 1))
    end
    half = n ÷ 2
    return collect((-half):half)
end

function _smallest_cutoff_from_tail(
    radii::AbstractVector{<:Real},
    magnitudes::AbstractVector{<:Real},
    tol::Real,
)
    length(radii) == length(magnitudes) || throw(ArgumentError("radii and magnitudes must have the same length"))
    isempty(radii) && throw(ArgumentError("spectrum cannot be empty"))
    0.0 < tol < 1.0 || throw(ArgumentError("tol must satisfy 0 < tol < 1"))

    max_mag = maximum(magnitudes)
    max_mag > 0 || return zero(Float64)

    order = sortperm(radii)
    sorted_radii = Float64.(radii[order])
    sorted_magnitudes = Float64.(magnitudes[order])
    suffix_exclusive = Vector{Float64}(undef, length(sorted_magnitudes))
    suffix_exclusive[end] = 0.0

    @inbounds for i in (length(sorted_magnitudes) - 1):-1:1
        suffix_exclusive[i] = max(suffix_exclusive[i + 1], sorted_magnitudes[i + 1])
    end

    i = 1
    while i <= length(sorted_radii)
        j = i
        while j < length(sorted_radii) && sorted_radii[j + 1] == sorted_radii[i]
            j += 1
        end
        if suffix_exclusive[i] / max_mag < tol
            return sorted_radii[i]
        end
        i = j + 1
    end

    return sorted_radii[end]
end

function _tail_ratio(radii::AbstractVector{<:Real}, magnitudes::AbstractVector{<:Real}, cutoff::Real)
    max_mag = maximum(magnitudes)
    max_mag > 0 || return 0.0
    tail_max = 0.0
    @inbounds for i in eachindex(radii)
        if radii[i] > cutoff && magnitudes[i] > tail_max
            tail_max = Float64(magnitudes[i])
        end
    end
    return tail_max / max_mag
end

function _load_squared_datagrid(path::AbstractString)
    structure, datagrid = BI.read_xsf(path)
    datagrid.values .*= datagrid.values
    return structure, datagrid
end

function _datagrid_geometry(datagrid)
    A, B, C = BI.true_cell_vectors(datagrid)
    cell = hcat(A, B, C)
    reciprocal = 2π .* inv(cell)'
    volume = abs(det(cell))
    return cell, reciprocal, volume
end

function _datagrid_nufft_inputs(datagrid, cell_volume::Real)
    nx, ny, nz = datagrid.nx, datagrid.ny, datagrid.nz
    npts = nx * ny * nz
    weight = cell_volume / npts

    x = Vector{Float64}(undef, npts)
    y = Vector{Float64}(undef, npts)
    z = Vector{Float64}(undef, npts)
    charges = Vector{ComplexF64}(undef, npts)

    idx = 1
    @inbounds for kz in 1:nz, jy in 1:ny, ix in 1:nx
        x[idx] = 2π * (ix - 1) / nx
        y[idx] = 2π * (jy - 1) / ny
        z[idx] = 2π * (kz - 1) / nz
        charges[idx] = ComplexF64(weight * datagrid.values[ix, jy, kz])
        idx += 1
    end

    return x, y, z, charges
end

function _spectral_radii_and_magnitudes(coeff::AbstractArray{<:Complex, 3}, reciprocal::AbstractMatrix{<:Real})
    nx, ny, nz = size(coeff)
    mx = _centered_mode_indices(nx)
    my = _centered_mode_indices(ny)
    mz = _centered_mode_indices(nz)

    npts = length(coeff)
    radii = Vector{Float64}(undef, npts)
    magnitudes = Vector{Float64}(undef, npts)

    b1x, b1y, b1z = reciprocal[1, 1], reciprocal[2, 1], reciprocal[3, 1]
    b2x, b2y, b2z = reciprocal[1, 2], reciprocal[2, 2], reciprocal[3, 2]
    b3x, b3y, b3z = reciprocal[1, 3], reciprocal[2, 3], reciprocal[3, 3]

    idx = 1
    @inbounds for kz in eachindex(mz), jy in eachindex(my), ix in eachindex(mx)
        kx = b1x * mx[ix] + b2x * my[jy] + b3x * mz[kz]
        ky = b1y * mx[ix] + b2y * my[jy] + b3y * mz[kz]
        kz_phys = b1z * mx[ix] + b2z * my[jy] + b3z * mz[kz]
        radii[idx] = sqrt(kx^2 + ky^2 + kz_phys^2)
        magnitudes[idx] = abs(coeff[ix, jy, kz])
        idx += 1
    end

    return radii, magnitudes
end

function _tail_checkpoints(max_radius::Real, cutoff::Real)
    checkpoints = Float64[0.0, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0] .* Float64(max_radius)
    push!(checkpoints, Float64(cutoff))
    return sort!(unique(checkpoints))
end

function orbital_spectrum_summary(path::AbstractString; tol::Real = TOL, eps::Real = NUFFT_EPS)
    structure, datagrid = _load_squared_datagrid(path)
    nx, ny, nz = datagrid.nx, datagrid.ny, datagrid.nz

    cell, reciprocal, cell_volume = _datagrid_geometry(datagrid)
    x, y, z, charges = _datagrid_nufft_inputs(datagrid, cell_volume)
    coeff = BI.TKM3D.nufft3d1(x, y, z, charges, -1, eps, nx, ny, nz)
    if ndims(coeff) == 4 && size(coeff, 4) == 1
        coeff = dropdims(coeff; dims = 4)
    end

    radii, magnitudes = _spectral_radii_and_magnitudes(coeff, reciprocal)
    cutoff = _smallest_cutoff_from_tail(radii, magnitudes, tol)
    max_radius = maximum(radii)
    max_magnitude = maximum(magnitudes)
    cutoff_tail_ratio = _tail_ratio(radii, magnitudes, cutoff)

    spacing_vectors = (cell[:, 1] ./ nx, cell[:, 2] ./ ny, cell[:, 3] ./ nz)
    axis_spacings = map(norm, spacing_vectors)
    axis_lengths = map(norm, eachcol(cell))
    reciprocal_norms = map(norm, eachcol(reciprocal))
    axis_nyquist = π ./ axis_spacings

    checkpoints = _tail_checkpoints(max_radius, cutoff)
    checkpoint_ratios = [(_tail_ratio(radii, magnitudes, k), k) for k in checkpoints]

    return (
        structure = structure,
        datagrid = datagrid,
        cell = cell,
        reciprocal = reciprocal,
        cell_volume = cell_volume,
        axis_lengths = collect(axis_lengths),
        axis_spacings = collect(axis_spacings),
        reciprocal_norms = collect(reciprocal_norms),
        axis_nyquist = collect(axis_nyquist),
        max_radius = max_radius,
        max_magnitude = max_magnitude,
        cutoff = cutoff,
        cutoff_tail_ratio = cutoff_tail_ratio,
        tol = tol,
        checkpoints = checkpoint_ratios,
    )
end

function _print_vec(label::AbstractString, v::AbstractVector{<:Real})
    @printf("%s = [%.8f, %.8f, %.8f]\n", label, v[1], v[2], v[3])
end

function main(; tol::Real = TOL, eps::Real = NUFFT_EPS)
    summary = orbital_spectrum_summary(ORBITAL_FILE; tol = tol, eps = eps)
    datagrid = summary.datagrid

    println("Orbital spectral cutoff analysis")
    println("file = $(ORBITAL_FILE)")
    @printf("grid = (%d, %d, %d)\n", datagrid.nx, datagrid.ny, datagrid.nz)
    _print_vec("|A|, |B|, |C|", summary.axis_lengths)
    _print_vec("grid spacing", summary.axis_spacings)
    _print_vec("|b1|, |b2|, |b3|", summary.reciprocal_norms)
    _print_vec("axis Nyquist", summary.axis_nyquist)
    @printf("cell volume = %.8f\n", summary.cell_volume)
    @printf("max |k| in box = %.8f\n", summary.max_radius)
    @printf("max |F(k)| = %.12e\n", summary.max_magnitude)
    @printf("tol = %.3e\n", summary.tol)
    @printf("estimated k_cut = %.8f\n", summary.cutoff)
    @printf("tail ratio at k_cut = %.12e\n", summary.cutoff_tail_ratio)

    println("tail checkpoints:")
    for (ratio, cutoff) in summary.checkpoints
        @printf("  k = %.8f -> max_{|q|>k}|F(q)| / max|F| = %.12e\n", cutoff, ratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
