include(joinpath(@__DIR__, "..", "theory", "density_power.jl"))

function _assert_close(actual, expected; rtol = 1e-10, atol = 1e-12)
    if !isapprox(actual, expected; rtol = rtol, atol = atol)
        error("expected $expected, got $actual")
    end
end

alphas = (10 / 33, 20 / 33, 30 / 33)
epsvals = epsilon_from_contrast.(alphas)
density_powers = density_power_right_angle_edge.(alphas)

_assert_close(epsvals[1], 43 / 23)
_assert_close(epsvals[2], 53 / 13)
_assert_close(epsvals[3], 21.0)

_assert_close(density_powers[1], -0.09683046687353458)
_assert_close(density_powers[2], -0.19599668215817823)
_assert_close(density_powers[3], -0.3003965754379143)

println("fig3 theory density checks passed")
