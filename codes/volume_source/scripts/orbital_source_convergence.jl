using BoundaryIntegral
import BoundaryIntegral as BI
using CSV, DataFrames
using Krylov
using Random

Random.seed!(1234)

L = 90.0
Lx = L
Ly = L
Lz = 9.0

orbital_file = joinpath(@__DIR__, "../../density_data/graphene_00002_cropped.xsf")
structure, datagrid = BI.read_xsf(orbital_file)
datagrid.values .*= datagrid.values
vs = BoundaryIntegral.VolumeSource(datagrid, shift = (0.0, 0.0, - 7.920155482424242))

eps_out = 1.0
eps_in = 6.0
max_order = 128

df = joinpath(@__DIR__, "data/orbital_source_convergence.csv")
CSV.write(df, DataFrame(p = Int[], r = Int[], loc_id = [], potential = []))

ps = [4, 6]
rs = collect(0:2:6)

c = (11.0, 12.0, 0.0)

locs = [c .+ (randn(), randn(), randn()) for _ in 1:10]

for (p, r) in [(x, y) for x in ps for y in rs] ∪ [(6, 7)]
    l_ec = 9.0 / 2^r * 1.01

    interface = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p, vs, eps_in, l_ec, 1e-6, eps_in, eps_out, Float64)

    rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, vs, eps_in, 1e-6)
    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(interface, 1e-6, 1e-6, 128, include_edges_src = false, include_edges_trg = false)
    sigma, status = Krylov.gmres(lhs, rhs, atol=1e-6, verbose = 1)

    for (loc_id, loc) in enumerate(locs)
        potential = BI.laplace3d_pottrg_near(interface, loc, sigma, 1e-8)
        CSV.write(df, DataFrame(p = p, r = r, loc_id = loc_id, potential = potential), append=true)
        @info "p = $p, r = $r, loc_id = $loc_id, potential = $potential"
    end

end