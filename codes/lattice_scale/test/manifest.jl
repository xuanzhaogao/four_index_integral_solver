@testset "manifest" begin
    # synthetic geometry: 2x2 cells, 2 sublattices, square lattice a1=(1,0,0) a2=(0,1,0)
    primvec = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 10.0]
    centroids = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]      # per template
    steps_per_cell = ((2, 0, 0), (0, 2, 0))             # grid steps for a1, a2

    centers = enumerate_centers(2, 2, primvec, centroids, steps_per_cell)
    @test length(centers) == 8
    @test allunique(getfield.(centers, :id))
    c1 = centers[findfirst(c -> c.Rx == 0 && c.Ry == 0 && c.template_id == 1, centers)]
    c5 = centers[findfirst(c -> c.Rx == 1 && c.Ry == 1 && c.template_id == 2, centers)]
    @test c1.steps == (0, 0, 0)
    @test c5.steps == (2, 2, 0)
    @test collect(c5.center) ≈ [1.5, 1.5, 0.0]

    # dedup: i <= j, each unordered pair once
    pairs = enumerate_pairs(centers, 1.2)                # nn cutoff: dist <= 1.2
    @test all(p -> p[1] <= p[2], pairs)
    @test allunique(pairs)
    @test length(pairs) == length(unique(pairs))
    onsite = count(p -> p[1] == p[2], pairs)
    @test onsite == 8                                    # every center pairs with itself

    batches = build_batches(pairs, 1)
    @test sum(b -> length(b.pairs), batches) == length(pairs)   # every pair exactly once
    @test allunique(getfield.(batches, :batch_id))
    # anchored: every pair's min id belongs to the batch's anchor set
    for b in batches
        @test all(p -> min(p[1], p[2]) in b.anchors, b.pairs)
    end
    b2 = build_batches(pairs, 2)
    @test sum(b -> length(b.pairs), b2) == length(pairs)
    @test length(b2) < length(batches)

    # TSV round-trip
    mktempdir() do dir
        cp = joinpath(dir, "centers.tsv"); mp = joinpath(dir, "manifest.tsv")
        write_centers(cp, centers); write_manifest(mp, batches)
        @test read_centers(cp) == centers
        @test read_manifest(mp) == batches
    end
end
