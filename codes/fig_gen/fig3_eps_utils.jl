using LinearAlgebra

function epsilon_from_contrast(alpha::Real; eps_0::Real = 1.0)
    if !(0 <= alpha < 1)
        throw(ArgumentError("contrast alpha must satisfy 0 <= alpha < 1"))
    end
    return eps_0 * (1 + alpha) / (1 - alpha)
end

function charge_integral(weights::AbstractVector, sigma::AbstractVector)
    length(weights) == length(sigma) ||
        throw(DimensionMismatch("weights and sigma must have the same length"))
    return dot(weights, sigma)
end

function charge_theory_interior_source(eps_d::Real; eps_0::Real = 1.0, source_charge::Real = 1.0)
    return source_charge * (1 - eps_0 / eps_d)
end

@inline function _dist3(a, b)
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]
    return sqrt(dx * dx + dy * dy + dz * dz)
end

function panel_min_edge_length(interface)
    l_min = Inf
    for panel in interface.panels
        a, b, c, d = panel.corners
        l_min = min(l_min, _dist3(a, b), _dist3(b, c), _dist3(c, d), _dist3(d, a))
    end
    return l_min
end
