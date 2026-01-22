# check the rhs-based refinement by integrating the total total_flux

using BoundaryIntegral
import BoundaryIntegral as BI
using CSV, DataFrames

function total_flux(box, source)
    flux = 0.0
    for trg in BI.eachpoint(box)
        p = trg.panel_point
        flux += BI.laplace3d_grad(source.point, p.point, p.normal) * p.weight
    end
    return flux
end

df = joinpath(@__DIR__, "data/total_flux.csv")
CSV.write(df, DataFrame(id = Int[], L = Float64[], r = Int[], p = Int[], tol = Float64[], total_flux = Float64[]))

for L in [5.0, 10.0, 20.0]
    for r in [2, 4, 6]
        for p in [2, 4, 6]
            for tol in [1e-2, 1e-4, 1e-6, 1e-8]
                for (id, ps) in enumerate([PointSource((0.1, 0.2, 1.0), 1.0), PointSource((L - 0.4, L - 0.5, 1.0), 1.0)])
                    Lx = L
                    Ly = L
                    Lz = 1.0

                    l_ec = 1.0 / 2^r * 1.01

                    eps_in = 4.0
                    eps_out = 1.0

                    # tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, ps, 1.0, l_ec, 1e-6, eps_in, eps_out, max_depth = 100)
                    tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, ps, 1.0, l_ec, tol, eps_in, eps_out)

                    tf = total_flux(tbox, ps)

                    CSV.write(df, DataFrame(id = id, L = L, r = r, p = p, tol = tol, total_flux = tf), append = true)
                end
            end
        end
    end
end