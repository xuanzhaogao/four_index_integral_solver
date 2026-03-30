using BoundaryIntegral
import BoundaryIntegral as BI
using CairoMakie

const SUBSTRATE_CENTER = (0.0, 0.0, -15.0)
const SUBSTRATE_SIZE = (30.0, 30.0, 30.0)
const FILM_CENTER = (0.0, 0.0, 1.4)
const FILM_SIZE = (20.0, 20.0, 2.8)

const EPSES = [10.0, 2.4]
const EPS_OUT = 1.0

const BOXES = [
    (center = SUBSTRATE_CENTER, Lx = SUBSTRATE_SIZE[1], Ly = SUBSTRATE_SIZE[2], Lz = SUBSTRATE_SIZE[3]),
    (center = FILM_CENTER, Lx = FILM_SIZE[1], Ly = FILM_SIZE[2], Lz = FILM_SIZE[3]),
]

# Representative mesh config from convergence study
const N_QUAD = 4
const EDGE_REFINE = 3

function build_interface()
    l_ec = FILM_SIZE[3] / 2.0^EDGE_REFINE * 1.01
    return BI.multi_dielectric_box3d(N_QUAD, l_ec, BOXES, EPSES, EPS_OUT)
end

function main()
    println("Building interface (n_quad=$N_QUAD, edge_refine=$EDGE_REFINE)...")
    interface = build_interface()
    println("  n_pts = $(BI.num_points(interface))")

    # Source location: orbital at center of graphene film
    src_x, src_y, src_z = 0.0, 0.0, 1.4

    fig = BI.viz_3d(interface;
                    show_points = false,
                    highlight_edges = true,
                    base_color = :steelblue,
                    edge_color = :orange,
                    size = (800, 700))

    # Retrieve the Axis3 from the figure to add the source marker and labels
    ax = nothing
    for obj in fig.content
        if obj isa Axis3
            ax = obj
            break
        end
    end

    if ax !== nothing
        # Mark the orbital source at the graphene center
        scatter!(ax, [src_x], [src_y], [src_z];
                 color = :red, markersize = 16, marker = :circle)
        ax.xlabel = "x"
        ax.ylabel = "y"
        ax.zlabel = "z"
        ax.title = "Graphene slab mesh  (n_quad=$N_QUAD, edge_refine=$EDGE_REFINE)"
    end

    mkpath("figs")
    save("figs/graphene_slab_mesh.png", fig)
    save("figs/graphene_slab_mesh.svg", fig)
    println("Saved figs/graphene_slab_mesh.{png,svg}")
end

main()
