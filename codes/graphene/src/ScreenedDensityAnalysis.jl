module ScreenedDensityAnalysis

using BoundaryIntegral
using FFTW
using FINUFFT
import BoundaryIntegral as BI

export box_bounds,
    datagrid_bounds,
    default_orbital_file,
    default_volume_source,
    density_centroid,
    centering_shift,
    kz_zero_mode_spectrum,
    linear_upsample_uniform,
    log10_clamped,
    normalized_global_nufft_qz_decay,
    nufft_upsample_periodic,
    nufft_upsample_global,
    normalized_kz_decay,
    normalized_refined_qz_decay,
    plane_charge_density,
    refined_kz_zero_mode_spectrum,
    screened_density_at_points,
    screened_density_vector,
    spectral_upsample_periodic,
    shift_datagrid,
    yz_slice_at_x

function default_orbital_file()
    return normpath(joinpath(@__DIR__, "..", "scripts", "../../../density_data/graphene_00001_5x5x1_shifted.xsf"))
end

function box_bounds(Lx::T, Ly::T, Lz::T) where {T}
    return ((-Lx / 2, -Ly / 2, -Lz / 2), (Lx / 2, Ly / 2, Lz / 2))
end

function datagrid_bounds(datagrid)
    o, M, _, u_max, v_max, w_max = BI._datagrid_affine(datagrid)
    corners = NTuple{3, Float64}[]
    for u in (0.0, u_max), v in (0.0, v_max), w in (0.0, w_max)
        push!(
            corners,
            (
                o[1] + M[1, 1] * u + M[1, 2] * v + M[1, 3] * w,
                o[2] + M[2, 1] * u + M[2, 2] * v + M[2, 3] * w,
                o[3] + M[3, 1] * u + M[3, 2] * v + M[3, 3] * w,
            ),
        )
    end
    xs = first.(corners)
    ys = getindex.(corners, 2)
    zs = last.(corners)
    return ((minimum(xs), minimum(ys), minimum(zs)), (maximum(xs), maximum(ys), maximum(zs)))
end

function shift_datagrid(datagrid, shift::NTuple{3, <:Real})
    origin = ntuple(i -> Float64(datagrid.origin[i]) + Float64(shift[i]), 3)
    return merge(datagrid, (; origin = origin))
end

function density_centroid(vs)
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

function centering_shift(datagrid; tol::Real = 1e-4)
    vs = BI.VolumeSource(datagrid, tol = tol)
    centroid = density_centroid(vs)
    return ntuple(i -> -Float64(centroid[i]), 3)
end

function default_volume_source(; orbital_file = default_orbital_file(), shift = nothing, tol = 1e-4)
    _, raw_datagrid = BI.read_xsf(orbital_file)
    squared_values = copy(raw_datagrid.values)
    squared_values .*= squared_values
    datagrid = merge(raw_datagrid, (; values = squared_values))
    applied_shift = isnothing(shift) ? centering_shift(datagrid; tol = tol) : shift
    shifted_datagrid = shift_datagrid(datagrid, applied_shift)
    vs = BI.VolumeSource(shifted_datagrid, tol = tol)
    return shifted_datagrid, vs
end

function screened_density_vector(
    vs,
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    mode;
    tol::Real = 1e-10,
)
    min_corner, max_corner = bounds
    T = eltype(vs.density)
    return BI._screened_volume_density(
        vs,
        ntuple(i -> T(min_corner[i]), 3),
        ntuple(i -> T(max_corner[i]), 3),
        T(eps_in),
        T(eps_out),
        mode,
        T(tol),
    )
end

function screened_density_at_points(
    datagrid,
    points::AbstractVector{<:NTuple{3, <:Real}},
    mode;
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    tol::Real = 1e-10,
)
    min_corner, max_corner = bounds
    _, _, Minv, _, _, _ = BI._datagrid_affine(datagrid)
    values = Vector{Float64}(undef, length(points))
    for i in eachindex(points)
        point = points[i]
        raw = BI._datagrid_trilinear_value(datagrid, Minv, point; tol = tol)
        if isnan(raw)
            values[i] = NaN
            continue
        end
        eps_local = BI._screened_permittivity(
            mode,
            point,
            ntuple(j -> Float64(min_corner[j]), 3),
            ntuple(j -> Float64(max_corner[j]), 3),
            Float64(eps_in),
            Float64(eps_out),
            Float64(tol),
        )
        values[i] = raw / eps_local
    end
    return values
end

function yz_slice_at_x(
    datagrid,
    x0::Real,
    mode;
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    ny::Integer = 250,
    nz::Integer = 250,
    y_limits = nothing,
    z_limits = nothing,
    tol::Real = 1e-10,
)
    (grid_min, grid_max) = datagrid_bounds(datagrid)
    y_min, y_max = isnothing(y_limits) ? (grid_min[2], grid_max[2]) : y_limits
    z_min, z_max = isnothing(z_limits) ? (grid_min[3], grid_max[3]) : z_limits
    ys = collect(range(y_min, y_max; length = ny))
    zs = collect(range(z_min, z_max; length = nz))

    points = NTuple{3, Float64}[]
    sizehint!(points, ny * nz)
    for y in ys, z in zs
        push!(points, (Float64(x0), Float64(y), Float64(z)))
    end

    values = screened_density_at_points(
        datagrid,
        points,
        mode;
        bounds = bounds,
        eps_in = eps_in,
        eps_out = eps_out,
        tol = tol,
    )

    slice = reshape(values, nz, ny)'
    return ys, zs, slice
end

function log10_clamped(values::AbstractArray{<:Real}; floor_value::Real)
    floor_value > 0 || throw(ArgumentError("floor_value must be positive"))
    transformed = Matrix{Float64}(undef, size(values)...)
    floor_f = Float64(floor_value)
    for i in eachindex(values)
        value = values[i]
        if isnan(value)
            transformed[i] = NaN
        else
            transformed[i] = log10(max(Float64(value), floor_f))
        end
    end
    return transformed
end

function plane_charge_density(vs, rho::AbstractVector{<:Real})
    length(rho) == length(vs.density) || throw(ArgumentError("rho must match vs.density length"))
    zs = sort(unique(vs.positions[3, :]))
    index_by_z = Dict(z => i for (i, z) in enumerate(zs))
    masses = zeros(Float64, length(zs))
    for i in eachindex(rho)
        masses[index_by_z[vs.positions[3, i]]] += vs.weights[i] * rho[i]
    end
    length(zs) >= 2 || throw(ArgumentError("at least two z planes are required"))
    dz = minimum(diff(zs))
    return zs, masses ./ dz, dz
end

function spectral_upsample_periodic(values::AbstractVector{<:Real}; factor::Integer)
    factor >= 1 || throw(ArgumentError("factor must be positive"))
    factor == 1 && return collect(Float64, values)

    n = length(values)
    m = n * factor
    coeffs = FFTW.fftshift(fft(values))
    padded = zeros(ComplexF64, m)
    start = div(m - n, 2) + 1
    padded[start:(start + n - 1)] .= coeffs
    refined = ifft(FFTW.ifftshift(padded)) .* factor
    return real.(refined)
end

function linear_upsample_uniform(values::AbstractVector{<:Real}; factor::Integer)
    factor >= 1 || throw(ArgumentError("factor must be positive"))
    factor == 1 && return collect(Float64, values)

    n = length(values)
    refined = Vector{Float64}(undef, (n - 1) * factor + 1)
    for i in eachindex(refined)
        t = (i - 1) / factor
        j = clamp(floor(Int, t) + 1, 1, n)
        if j == n
            refined[i] = Float64(values[end])
        else
            α = t - (j - 1)
            refined[i] = (1 - α) * Float64(values[j]) + α * Float64(values[j + 1])
        end
    end
    return refined
end

function nufft_upsample_periodic(values::AbstractVector{<:Real}; factor::Integer, eps::Real = 1e-12)
    factor >= 1 || throw(ArgumentError("factor must be positive"))
    factor == 1 && return collect(Float64, values)

    n = length(values)
    x = 2π .* collect(0:(n - 1)) ./ n
    coeffs = FINUFFT.nufft1d1(x, ComplexF64.(values), -1, eps, n) ./ n
    x_refined = 2π .* collect(0:(n * factor - 1)) ./ (n * factor)
    refined = FINUFFT.nufft1d2(x_refined, 1, eps, coeffs)
    return real.(refined)
end

function nufft_upsample_global(
    values::AbstractVector{<:Real};
    factor::Integer,
    pad_factor::Integer = 4,
    eps::Real = 1e-12,
)
    factor >= 1 || throw(ArgumentError("factor must be positive"))
    pad_factor >= 1 || throw(ArgumentError("pad_factor must be positive"))
    factor == 1 && pad_factor == 1 && return collect(Float64, values)

    n = length(values)
    n_padded = n * pad_factor
    left_pad = fld(n_padded - n, 2)
    padded_values = zeros(ComplexF64, n_padded)
    padded_values[(left_pad + 1):(left_pad + n)] .= ComplexF64.(values)

    x = 2π .* collect(0:(n_padded - 1)) ./ n_padded
    coeffs = FINUFFT.nufft1d1(x, padded_values, -1, eps, n_padded) ./ n_padded
    x_refined = 2π .* collect(0:(n_padded * factor - 1)) ./ (n_padded * factor)
    refined_full = real.(FINUFFT.nufft1d2(x_refined, 1, eps, coeffs))

    left_pad_refined = left_pad * factor
    return refined_full[(left_pad_refined + 1):(left_pad_refined + n * factor)]
end

function z_plane_charges(vs, rho::AbstractVector{<:Real})
    length(rho) == length(vs.density) || throw(ArgumentError("rho must match vs.density length"))
    zs = unique(vs.positions[3, :])
    index_by_z = Dict(z => i for (i, z) in enumerate(zs))
    charges = zeros(ComplexF64, length(zs))
    for i in eachindex(rho)
        charges[index_by_z[vs.positions[3, i]]] += vs.weights[i] * rho[i]
    end
    return zs, charges
end

function kz_zero_mode_spectrum(vs, rho::AbstractVector{<:Real}, kz_values::AbstractVector{<:Real})
    zs, charges = z_plane_charges(vs, rho)
    return ComplexF64[
        sum(charges .* exp.(-im * kz .* zs))
        for kz in kz_values
    ]
end

function refined_kz_zero_mode_spectrum(
    vs,
    rho::AbstractVector{<:Real},
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    mode,
    kz_values::AbstractVector{<:Real};
    upsample_factor::Integer = 8,
    tol::Real = 1e-10,
)
    min_corner, max_corner = bounds
    zs, qz, dz = plane_charge_density(vs, rho)
    qz_refined = nufft_upsample_periodic(qz; factor = upsample_factor, eps = tol)
    dz_refined = dz / upsample_factor
    z0 = first(zs)
    z_refined = z0 .+ dz_refined .* collect(0:(length(qz_refined) - 1))
    screened_masses = ComplexF64[]
    sizehint!(screened_masses, length(z_refined))
    for (z, q) in zip(z_refined, qz_refined)
        eps_local = BI._screened_permittivity(
            mode,
            (0.0, 0.0, z),
            ntuple(i -> Float64(min_corner[i]), 3),
            ntuple(i -> Float64(max_corner[i]), 3),
            Float64(eps_in),
            Float64(eps_out),
            Float64(tol),
        )
        push!(screened_masses, (q / eps_local) * dz_refined)
    end

    return ComplexF64[
        sum(screened_masses .* exp.(-im * kz .* z_refined))
        for kz in kz_values
    ]
end

function global_nufft_kz_zero_mode_spectrum(
    vs,
    rho::AbstractVector{<:Real},
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    mode,
    kz_values::AbstractVector{<:Real};
    upsample_factor::Integer = 8,
    pad_factor::Integer = 4,
    tol::Real = 1e-10,
)
    min_corner, max_corner = bounds
    zs, qz, dz = plane_charge_density(vs, rho)
    qz_refined = nufft_upsample_global(qz; factor = upsample_factor, pad_factor = pad_factor, eps = tol)
    dz_refined = dz / upsample_factor
    z0 = first(zs)
    z_refined = z0 .+ dz_refined .* collect(0:(length(qz_refined) - 1))
    screened_masses = ComplexF64[]
    sizehint!(screened_masses, length(z_refined))
    for (z, q) in zip(z_refined, qz_refined)
        eps_local = BI._screened_permittivity(
            mode,
            (0.0, 0.0, z),
            ntuple(i -> Float64(min_corner[i]), 3),
            ntuple(i -> Float64(max_corner[i]), 3),
            Float64(eps_in),
            Float64(eps_out),
            Float64(tol),
        )
        push!(screened_masses, (q / eps_local) * dz_refined)
    end

    return ComplexF64[
        sum(screened_masses .* exp.(-im * kz .* z_refined))
        for kz in kz_values
    ]
end

function normalized_kz_decay(vs, rho::AbstractVector{<:Real}; num_points::Integer = 400, kmax = nothing)
    zs = sort(unique(vs.positions[3, :]))
    dz = minimum(diff(zs))
    kmax_value = isnothing(kmax) ? (π / dz) : kmax
    kz_values = collect(range(0.0, Float64(kmax_value); length = num_points))
    spectrum = kz_zero_mode_spectrum(vs, rho, kz_values)
    scale = abs(spectrum[1])
    if iszero(scale)
        scale = 1.0
    end
    return kz_values, abs.(spectrum) ./ scale
end

function normalized_refined_qz_decay(
    vs,
    rho::AbstractVector{<:Real},
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    mode;
    upsample_factor::Integer = 8,
    num_points::Integer = 400,
    kmax = nothing,
    tol::Real = 1e-10,
)
    zs = sort(unique(vs.positions[3, :]))
    dz = minimum(diff(zs))
    kmax_value = isnothing(kmax) ? (π / dz) : kmax
    kz_values = collect(range(0.0, Float64(kmax_value); length = num_points))
    spectrum = refined_kz_zero_mode_spectrum(
        vs,
        rho,
        bounds,
        eps_in,
        eps_out,
        mode,
        kz_values;
        upsample_factor = upsample_factor,
        tol = tol,
    )
    scale = abs(spectrum[1])
    if iszero(scale)
        scale = 1.0
    end
    return kz_values, abs.(spectrum) ./ scale
end

function normalized_global_nufft_qz_decay(
    vs,
    rho::AbstractVector{<:Real},
    bounds::Tuple{<:NTuple{3, <:Real}, <:NTuple{3, <:Real}},
    eps_in::Real,
    eps_out::Real,
    mode;
    upsample_factor::Integer = 8,
    pad_factor::Integer = 4,
    num_points::Integer = 400,
    kmax = nothing,
    tol::Real = 1e-10,
)
    zs = sort(unique(vs.positions[3, :]))
    dz = minimum(diff(zs))
    kmax_value = isnothing(kmax) ? (π / dz) : kmax
    kz_values = collect(range(0.0, Float64(kmax_value); length = num_points))
    spectrum = global_nufft_kz_zero_mode_spectrum(
        vs,
        rho,
        bounds,
        eps_in,
        eps_out,
        mode,
        kz_values;
        upsample_factor = upsample_factor,
        pad_factor = pad_factor,
        tol = tol,
    )
    scale = abs(spectrum[1])
    if iszero(scale)
        scale = 1.0
    end
    return kz_values, abs.(spectrum) ./ scale
end

end
