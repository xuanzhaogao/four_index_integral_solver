using BoundaryIntegral
using GLMakie

const DENSITY_DIR = normpath(joinpath(@__DIR__, "../../density_data"))
const OUT_DIR = joinpath(@__DIR__, "figs", "zslice_verify")
mkpath(OUT_DIR)

files = sort(filter(f -> endswith(lowercase(f), ".xsf"), readdir(DENSITY_DIR; join = true)))
isempty(files) && error("No XSF files found in $(DENSITY_DIR)")

println("Verifying ", length(files), " XSF files from: ", DENSITY_DIR)

for f in files
    base = splitext(basename(f))[1]
    _, datagrid = BoundaryIntegral.read_xsf(f)

    # Verify datagrid slicing behavior
    s_mid = BoundaryIntegral.datagrid_zslice(datagrid)
    @assert size(s_mid.values) == (datagrid.nx, datagrid.ny)

    iz = clamp(cld(datagrid.nz, 2), 1, datagrid.nz)
    s_i = BoundaryIntegral.datagrid_zslice(datagrid; iz = iz)
    @assert s_i.iz == iz

    z_target = BoundaryIntegral.grid_point(datagrid, 1, 1, iz)[3]
    s_z = BoundaryIntegral.datagrid_zslice(datagrid; z = z_target)
    @assert s_z.iz == iz

    # Verify GLMakie rendering path
    fig = BoundaryIntegral.viz_3d_zslice(datagrid; iz = iz, add_colorbar = true)
    out = joinpath(OUT_DIR, base * "_zslice_k$(iz).png")
    save(out, fig)
    @assert isfile(out)

    println("OK: ", basename(f), " -> ", basename(out), " (nx=$(datagrid.nx), ny=$(datagrid.ny), nz=$(datagrid.nz), iz=$(iz))")
end

println("All XSF + GLMakie z-slice verifications passed")
