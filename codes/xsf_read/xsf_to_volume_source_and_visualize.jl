using BoundaryIntegral
using GLMakie

const XSF_PATH = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002_cropped.xsf"
# const OUT_PATH = "/Users/xgao/Works/four_indices_integral_solver/codes/xsf_read/graphene_00002_volume_source.png"

structure, datagrid = BoundaryIntegral.read_xsf(XSF_PATH)
# datagrid.values .*= datagrid.values

vs = BoundaryIntegral.VolumeSource(datagrid)

fig = BoundaryIntegral.viz_3d(
    sources = vs;
    size = (1200, 850),
    markersize = 2.5,
    alpha = 0.7,
    log_density = false,
    add_colorbar = true,
    max_points = 140_000,
)

ax = content(fig[1, 1])
ax.xlabel = "x"
ax.ylabel = "y"
ax.zlabel = "z"

atom_pos = structure.positions
scatter!(ax, atom_pos[:, 1], atom_pos[:, 2], atom_pos[:, 3]; color = :black, markersize = 8)

save(OUT_PATH, fig)
println("Saved: ", OUT_PATH)
println("Grid size: ", size(vs.density))
