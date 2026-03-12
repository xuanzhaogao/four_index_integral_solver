using Test

include(joinpath(@__DIR__, "..", "scripts", "max_frequency.jl"))

@testset "smallest cutoff from relative pointwise tail" begin
    radii = [0.0, 1.0, 1.0, 2.0, 3.0, 3.0]
    magnitudes = [10.0, 5.0, 4.0, 0.8, 0.09, 0.08]

    cutoff = _smallest_cutoff_from_tail(radii, magnitudes, 0.1)

    @test cutoff == 2.0
end
