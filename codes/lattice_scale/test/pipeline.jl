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
        k2 = length(br1.pair_ids) + 1                  # first pair of batch 2
        @test store.pair_ids[k2] == br2.pair_ids[1]
        @test T.gidx[store.t_idx[k2]] == br2.gidx
        # consolidate is idempotent
        consolidate(c)
        @test store.pair_ids == open(Serialization.deserialize, rho_store_path(c)).pair_ids
    end
end

@testset "eval_batch: V file + plumbing + four_index sanity" begin
    mktempdir() do dir
        c = load_campaign(write_fixture_campaign(dir))
        prepare(c); solve_batch(c, 1); solve_batch(c, 2); consolidate(c)
        eval_batch(c, 1)
        @test CampaignLib._is_complete_v(v_path(c, 1))
        vr = CampaignLib.load_v_rows(v_path(c, 1))
        @test vr.source_pairs == [(1, 1), (1, 2)]
        @test vr.target_pairs == [(1, 1), (1, 2), (2, 2)]   # ALL pairs, manifest order
        @test size(vr.V) == (3, 2)                          # n_target_pairs × K_sources

        br = load_batch_result(batch_path(c, 1))
        T = open(Serialization.deserialize, targets_path(c))
        store = open(Serialization.deserialize, rho_store_path(c))
        temps = CampaignLib.load_templates!(c); dg = temps[1][2]
        pos = Matrix{Float64}(undef, 3, length(br.gidx))
        for (r, g) in enumerate(br.gidx)
            p = BoundaryIntegral.grid_point(dg, g[1], g[2], g[3]); pos[:, r] .= p
        end
        srcs = [BoundaryIntegral.VolumeSource(copy(pos), copy(br.weights), br.densities[:, k])
                for k in 1:2]
        At, Bt, Ct = BoundaryIntegral.true_cell_vectors(dg)
        far_pad = c.far_pad_steps * maximum((LinearAlgebra.norm(At)/dg.nx,
            LinearAlgebra.norm(Bt)/dg.ny, LinearAlgebra.norm(Ct)/dg.nz))

        # (1) PLUMBING (tight): reproduce eval_batch's OWN computation — evaluate Φ at the
        # shared target set T and contract via store.t_idx/store.tw. Pins source
        # reconstruction, the target-set Φ evaluation, the contraction indexing, and the
        # V-file round-trip, independent of TKM's absolute accuracy.
        ΦT = evaluate_batch_potential(br.interface, br.sigma, srcs, T.positions;
            lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], far_pad = far_pad)
        V_plumb = [LinearAlgebra.dot(store.tw[kl], view(ΦT, store.t_idx[kl], a))
                   for kl in 1:length(store.pair_ids), a in 1:2]
        @test maximum(abs.(vr.V .- V_plumb)) < 1e-10 * max(maximum(abs.(V_plumb)), eps())

        # (2) PHYSICAL SANITY (loose): within-batch block vs four_index_matrix-style
        # contraction at the batch's OWN grid. The fixture is deliberately under-resolved,
        # and TKM's k-grid depends on the source+target bbox (Task 4 finding), so
        # eval-at-T differs from eval-at-own-grid by ~1e-3 here; on production-resolution
        # data this collapses to ~1e-8. The tight accuracy anchor is Task 16
        # (scripts/compare_anchor.jl on real Wannier data via the .bie path).
        Φg = evaluate_batch_potential(br.interface, br.sigma, srcs, pos;
            lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], far_pad = far_pad)
        V_fi = [LinearAlgebra.dot(br.weights .* br.densities[:, a], Φg[:, bb])
                for a in 1:2, bb in 1:2]
        # vr.V rows = target pairs, cols = source pairs; the first 2 target rows are this
        # batch's pairs, so vr.V[1:2, :] aligns with V_fi[target a, source bb].
        @test maximum(abs.(vr.V[1:2, :] .- V_fi)) < 1e-2 * maximum(abs.(V_fi))

        eval_batch(c, 2)
        @test isempty(pending_batches(c, :eval))
    end
end
