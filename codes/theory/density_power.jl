using LinearAlgebra
using Roots

function epsilon_from_contrast(alpha::Real; eps_0::Real = 1.0)
    if !(0 <= alpha < 1)
        throw(ArgumentError("contrast alpha must satisfy 0 <= alpha < 1"))
    end
    return eps_0 * (1 + alpha) / (1 - alpha)
end

"""
    theta_ODE_det(a, e, gamma)

Transmission determinant for the dielectric angular ODE on a multi-material
wedge. The vector `a` contains all but the final material angle; the final
angle is `2π - sum(a)`. The vector `e` contains relative permittivities.
Zeros in `gamma` are wedge potential powers.
"""
function theta_ODE_det(a::AbstractVector, e::AbstractVector, gamma::Real)
    M = Matrix(1.0I, 2, 2)
    aa = [a; 2pi - sum(a)]
    @assert aa[end] >= 0 "angles sum to more than 2π"
    @assert length(aa) == length(e) "angle and permittivity vectors have inconsistent lengths"

    next_e = circshift(e, -1)
    for (j, eps_j) in enumerate(e)
        if !isinf(eps_j)
            c = cos(gamma * aa[j])
            s = sin(gamma * aa[j])
            s_over_gamma = gamma == 0.0 ? aa[j] : s / gamma
            M = [1.0 0.0; 0.0 eps_j / next_e[j]] *
                [c s_over_gamma; -gamma * s c] * M
        end
    end
    return isinf(e[end]) ? M[1, 2] : det(M - I)
end

"""
    density_power_right_angle_edge(alpha)

Return the leading surface-density power for a dielectric right-angle edge
with contrast `alpha = (eps_d - eps_0) / (eps_d + eps_0)`.

The cube edge cross-section is modeled as a dielectric wedge of angle π/2
embedded in vacuum. The SLP density scales like `r^(gamma - 1)`, so this
function returns `gamma - 1`.
"""
function density_power_right_angle_edge(alpha::Real)
    eps_d = epsilon_from_contrast(alpha)
    gamma = fzero(g -> theta_ODE_det([3pi / 2], [1.0, eps_d], g), 0.8)
    return gamma - 1.0
end
