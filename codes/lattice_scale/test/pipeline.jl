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

@testset "solve_batch" begin
    mktempdir() do dir
        c = load_campaign(write_fixture_campaign(dir))
        prepare(c)
        t = @elapsed solve_batch(c, 1)
        @test isfile(batch_path(c, 1))
        br = load_batch_result(batch_path(c, 1))
        @test br.batch_id == 1
        @test br.pair_ids == [(1, 1), (1, 2)]
        @test size(br.sigma, 2) == 2
        @test size(br.densities, 2) == 2 && size(br.densities, 1) == length(br.gidx)
        @test all(haskey(br.stats, k) for k in
                  ("t_setup", "t_assemble", "t_solve", "t_total", "niter", "dof", "n_support", "K", "hostname"))
        @test br.stats["dof"] > 0
        solve_batch(c, 1)                                # idempotent: skips, no error
        @test sort(pending_batches(c, :solve)) == [2]
        solve_batch(c, 2)
        @test isempty(pending_batches(c, :solve))
    end
end

@testset "consolidate" begin
    mktempdir() do dir
        c = load_campaign(write_fixture_campaign(dir))
        prepare(c); solve_batch(c, 1); solve_batch(c, 2)
        consolidate(c)
        @test isfile(targets_path(c)) && isfile(rho_store_path(c))

        T = open(Serialization.deserialize, targets_path(c))
        store = open(Serialization.deserialize, rho_store_path(c))
        br1 = load_batch_result(batch_path(c, 1))
        br2 = load_batch_result(batch_path(c, 2))

        @test T.gidx == sort(union(br1.gidx, br2.gidx))          # exact union, sorted
        @test size(T.positions) == (3, length(T.gidx))
        @test store.pair_ids == vcat(br1.pair_ids, br2.pair_ids) # batch order
        # contraction vectors: tw = w .* rho on the pair's own support rows
        k = 1                                                    # pair (1,1) from batch 1
        @test store.tw[k] ≈ (br1.weights .* br1.densities[:, 1])
        # t_idx maps the pair's support into T
        @test T.gidx[store.t_idx[k]] == br1.gidx
        # consolidate is idempotent
        consolidate(c)
        @test store.pair_ids == open(Serialization.deserialize, rho_store_path(c)).pair_ids
    end
end
