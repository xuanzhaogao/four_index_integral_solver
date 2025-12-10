using FastGaussQuadrature
using BoundaryIntegral
import BoundaryIntegral as BI
using CSV, DataFrames

# target point is at (1.0, 0.0, δ), norm is (1.0, 0.0, 0.0)
# surface is at (-1.0, -1.0, 0.0) to (1.0, 1.0, 0.0)

function near_int_surface(δ, p)

    ns, ws = gausslegendre(p)
    t = 0.0

    for i in 1:length(ns)
        for j in 1:length(ns)
            r = sqrt((1.0 - ns[i])^2 + (ns[j])^2 + δ^2)
            Ex = (1.0 - ns[i]) / r^3
            t += Ex * ws[i] * ws[j]
        end
    end

    return t
end

function main()

    df = joinpath(@__DIR__, "data/near_int_surface.csv")
    CSV.write(df, DataFrame(δ = [], p = [], t = []))

    for δ in 0.01:0.01:0.1
        for p in 2:2:256
            t = near_int_surface(δ, p)
            println("δ = $δ, p = $p, t = $t")
            CSV.write(df, DataFrame(δ = [δ], p = [p], t = [t]), append = true)
        end
    end
end

function near_int_surface_error_gl(p)
    ns, ws = gausslegendre(p)
    t = 0.0
    d = 1 + ns[1]

    for i in 1:length(ns)
        for j in 1:length(ns)
            r = sqrt((1.0 - ns[i])^2 + (ns[j])^2 + d^2)
            Ex = (1.0 - ns[i]) / r^3
            t += Ex * ws[i] * ws[j]
        end
    end

    return t, d
end

function main_2()
    ps = []
    ds = []
    ts = []
    ref_ts = []
    self_errs = []

    for p in 2:2:256
        t, d = near_int_surface_error_gl(p)
        push!(ps, p)
        push!(ds, d)
        push!(ts, t)
        ref_t = near_int_surface(d, 2048)
        ref_t2 = near_int_surface(d, 4096)
        push!(ref_ts, ref_t)
        push!(self_errs, abs.(ref_t - ref_t2))
        println("p = $p, d = $d, t = $t, ref_t = $ref_t, self_err = $(abs.(ref_t - ref_t2))")
    end

    df = DataFrame(p = ps, d = ds, t = ts, ref_t = ref_ts, self_err = self_errs)
    CSV.write(joinpath(@__DIR__, "data/near_int_surface_error_gl.csv"), df)

end

main()
main_2()