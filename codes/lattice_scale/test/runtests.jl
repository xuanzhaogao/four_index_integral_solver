using CampaignLib
using BoundaryIntegral
using LinearAlgebra
using Serialization
using Test

@testset "CampaignLib" begin
    include("config.jl")
    include("manifest.jl")
    include("pipeline.jl")
end
