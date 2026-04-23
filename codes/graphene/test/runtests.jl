using Test

@testset "bilayer_slab" begin
    include(joinpath(@__DIR__, "..", "bilayer_slab", "test", "runtests.jl"))
end

@testset "monolayer" begin
    monolayer_runtests = joinpath(@__DIR__, "..", "monolayer", "test", "runtests.jl")
    if isfile(monolayer_runtests)
        include(monolayer_runtests)
    else
        @info "monolayer test suite not yet present; skipping"
    end
end
