include(joinpath(@__DIR__, "single_box3d_utils.jl"))

Ls = [5.0, 10.0, 20.0]
as = [0.5]
# ps = [1, 2, 4]
ps_min = [2, 4]
ps_max = [8, 16]
rs = 0:2:6
eps = 4.0
src = (0.2, 0.3, 0.4)


df = joinpath(@__DIR__, "data/single_thin_box_reduce_convergence.csv")
CSV.write(df, DataFrame(L = [], nxy = [], nz = [], p_max = [], r = [], gi = [], n_val = [], n_iter = []))

for L in Ls
    for a in as
        nxy = ceil(Int, L / a)
        nz = ceil(Int, 1.0 / a)
        for p_min in ps_min
            for p_max in ps_max
                for r in rs
                    println("L = $L, nxy = $nxy, nz = $nz, p_min = $p_min, p_max = $p_max, r = $r")
                    tbox, sigma, gi, n_val, n_iter = solve_single_thin_box3d(eps, L, L, 1.0, nxy, nxy, nz, p_min, p_max, r, r, src)
                    CSV.write(df, DataFrame(L = [L], nxy = [nxy], nz = [nz], p_min = [p_min], p_max = [p_max], r = [r], gi = [gi], n_val = [n_val], n_iter = [n_iter]), append = true)
                    res = Dict("tbox" => tbox, "sigma" => sigma, "gi" => gi)
                    name = "single_thin_box_reduce_convergence_L$(Int(L))_r$(Int(r))_p_min$(p_min)_p_max$(p_max)_r$(Int(r)).jld2"
                    save(joinpath(@__DIR__, "cache/$name"), res)
                end
            end
        end
    end
end
