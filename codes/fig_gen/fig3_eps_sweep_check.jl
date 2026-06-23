using LinearAlgebra

include("fig3_eps_utils.jl")

function _assert_close(actual, expected; rtol = 1e-12, atol = 1e-12)
    if !isapprox(actual, expected; rtol = rtol, atol = atol)
        error("expected $expected, got $actual")
    end
end

alphas = (10 / 33, 20 / 33, 30 / 33)
epsvals = epsilon_from_contrast.(alphas)

_assert_close(epsvals[1], 43 / 23)
_assert_close(epsvals[2], 53 / 13)
_assert_close(epsvals[3], 21.0)

for (alpha, eps_d) in zip(alphas, epsvals)
    _assert_close((eps_d - 1.0) / (eps_d + 1.0), alpha)
end

theory_charges = charge_theory_interior_source.(epsvals)
_assert_close(theory_charges[1], 20 / 43)
_assert_close(theory_charges[2], 40 / 53)
_assert_close(theory_charges[3], 20 / 21)

weights = [0.25, 0.25, 0.5]
sigma = [2.0, -4.0, 1.0]
_assert_close(charge_integral(weights, sigma), 0.0)

println("fig3 epsilon sweep checks passed")
