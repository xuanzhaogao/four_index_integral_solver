using Serialization

include(joinpath(@__DIR__, "..", "theory", "density_power.jl"))

const contrast_specs = [(num = 10, den = 33), (num = 20, den = 33), (num = 30, den = 33)]

theory_results = [
    begin
        alpha = spec.num / spec.den
        eps_d = epsilon_from_contrast(alpha)
        power = density_power_right_angle_edge(alpha)
        (
            contrast_num = spec.num,
            contrast_den = spec.den,
            alpha = alpha,
            eps_d = eps_d,
            density_power = power,
        )
    end
    for spec in contrast_specs
]

out = (
    model = "right-angle dielectric edge",
    density_scaling = "sigma ~ r^density_power",
    theory_results = theory_results,
)

path = joinpath(@__DIR__, "fig3_theory_density.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path

for r in theory_results
    println("alpha=$(r.contrast_num)/$(r.contrast_den) eps=$(r.eps_d) density_power=$(r.density_power)")
end
