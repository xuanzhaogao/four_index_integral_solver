@testset "config" begin
    mktempdir() do dir
        toml = joinpath(dir, "c.toml")
        write(toml, """
        [campaign]
        name = "t"
        root = "$dir/out"

        [orbitals]
        xsf = ["/a/one.xsf", "/a/two.xsf"]

        [lattice]
        nx = 2
        ny = 3
        neighbor_cutoff = 5.0

        [dielectrics]
        eps_out = 1.0
        boxes = [[0.0, 0.0, 7.5, 90.0, 90.0, 3.35, 3.5]]

        [solve]
        n_quad = 6
        edge_refine_level = 2
        rhs_tol = 1e-3
        lhs_tol = 1e-5
        gmres_rtol = 1e-5
        support_rtol = 1e-4
        volume_tol = 1e-5
        max_order = 8
        max_depth = 128

        [batching]
        n_centers_per_batch = 1

        [eval]
        far_pad_steps = 2.0
        """)
        c = load_campaign(toml)
        @test c.name == "t" && c.nx == 2 && c.ny == 3
        @test length(c.boxes) == 1 && c.epses == [3.5]
        @test c.boxes[1].Lz == 3.35
        @test c.solve["rhs_tol"] == 1e-3
        @test endswith(batch_path(c, 7), "batches/batch_0007.jls")
        @test endswith(v_path(c, 12), "V/V_0012.jls")
        @test CampaignLib.campaign_l_ec(c) ≈ 3.35 / 4 * 1.01
    end
end
