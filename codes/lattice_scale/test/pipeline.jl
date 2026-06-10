include("fixture_campaign.jl")

@testset "prepare + status" begin
    mktempdir() do dir
        c = load_campaign(write_fixture_campaign(dir))
        prepare(c)
        @test isfile(manifest_path(c)) && isfile(centers_path(c))
        centers = read_centers(centers_path(c))
        @test length(centers) == 2                       # 2x1 lattice, 1 sublattice
        @test centers[2].steps == (2, 0, 0)              # 2 steps/cell along a1
        batches = read_manifest(manifest_path(c))
        @test length(batches) == 2                       # anchors 1 and 2
        @test sort(reduce(vcat, b.pairs for b in batches)) == [(1,1), (1,2), (2,2)]

        @test sort(pending_batches(c, :solve)) == [1, 2]
        # a completed batch file flips status
        br = BatchResult(BoundaryIntegral.BATCH_FORMAT_VERSION, 1, [(1,1),(1,2)],
            [(1,1,1)], [0.1], ones(1, 2), nothing, ones(3, 2), Dict{String,Any}())
        save_batch_result(batch_path(c, 1), br)
        @test sort(pending_batches(c, :solve)) == [2]
        @test sort(pending_batches(c, :eval)) == [1, 2]  # no V files yet

        # re-running prepare with same params must NOT throw (idempotent)
        @test prepare(c) isa Vector{BatchSpec}

        # pending_batches without prepare errors out
        mktempdir() do dir2
            c2 = load_campaign(write_fixture_campaign(dir2))
            @test_throws ErrorException pending_batches(c2, :solve)
        end
    end
end
