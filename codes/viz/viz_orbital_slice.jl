using BoundaryIntegral
import BoundaryIntegral as BI
using GLMakie

for i in [1, 2]

    orbital_file = joinpath(@__DIR__, "../../density_data/graphene_0000$(i)_5x5x1.xsf")
    structure, datagrid = BI.read_xsf(orbital_file)

    for z in 0.0:0.1:15.8
        fig = BI.viz_3d_zslice(
            datagrid;
            z = z,
            interpolation = :trilinear,
            nx_sample = 200,
            ny_sample = 200,
            add_colorbar = true,
            log_density = false,
        )

        save(joinpath(@__DIR__, "zslices/g$(i)_orbital_slice_z$(z).png"), fig)

    end
end