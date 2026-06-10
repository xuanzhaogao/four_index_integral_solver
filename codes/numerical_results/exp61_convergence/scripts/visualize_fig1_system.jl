# Visualize the Fig.-1 system geometry: 10x10x1 slab on two 10x10x10 cubes.
# Builds the COARSE panelization with the actual multi-box builder (so shared
# faces / contact regions are exactly what the solver will see), renders the
# panels in 3D colored by interface type, and prints a face census.
# Output: figs/fig61_system_fig1_geometry.png (+ .pdf)

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using CairoMakie, Printf
using GeometryBasics: Point3f, TriangleFace

const SYS = system_fig1()
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

# coarse, non-adaptive build: p = 2, one edge-refinement level
iface = BI.multi_dielectric_box3d(2, Harness.l_ec_of(SYS, 1), SYS.boxes, SYS.epses, SYS.eps_out)
panels = iface.panels
println("panels: ", length(panels), "   points: ", BI.num_points(iface))

# classify each panel by its (eps_in, eps_out) pair
pairkey(i) = (round(min(iface.eps_in[i], iface.eps_out[i]); digits = 6),
              round(max(iface.eps_in[i], iface.eps_out[i]); digits = 6))
LABEL = Dict(
    (1.0, 4.0) => "Ω₁ | vac",
    (1.0, 12.0) => "Ω₂ | vac",
    (4.0, 12.0) => "Ω₁ | Ω₂  (shared)",
    (1.0, 10.0) => "slab | vac",
    (4.0, 10.0) => "slab | Ω₁ (contact)",
    (10.0, 12.0) => "slab | Ω₂ (contact)",
)
COLOR = Dict(
    (1.0, 4.0) => "#E8C547",      # warm yellow
    (1.0, 12.0) => "#9aa0a6",     # gray
    (4.0, 12.0) => "#CC79A7",     # magenta (internal cube-cube face)
    (1.0, 10.0) => "#56B4E9",     # light blue slab
    (4.0, 10.0) => "#D55E00",     # orange contact
    (10.0, 12.0) => "#009E73",    # green contact
)

census = Dict{Tuple{Float64, Float64}, Int}()
for i in eachindex(panels)
    census[pairkey(i)] = get(census, pairkey(i), 0) + 1
end
println("\nface census (eps pair => panel count):")
for (k, v) in sort(collect(census))
    @printf("  %-22s  %5d panels\n", get(LABEL, k, string(k)), v)
end

# duplicate-center scan (double-counted shared faces would show here)
let seen = Dict{NTuple{3, Float64}, Int}(), dup = 0
    for p in panels
        c = (p.corners[1] .+ p.corners[2] .+ p.corners[3] .+ p.corners[4]) ./ 4
        key = (round(c[1]; digits = 9), round(c[2]; digits = 9), round(c[3]; digits = 9))
        dup += (seen[key] = get(seen, key, 0) + 1) > 1 ? 1 : 0
    end
    println("duplicate-center panels: ", dup, " (expected 0)")
end

# ---------------------------------------------------------------- 3D rendering
fig = Figure(size = (980, 560))
ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.72π, elevation = 0.16π,
    title = "Fig.-1 system: 10×10×1 slab (ε=10) on two 10×10×10 cubes (ε=4, ε=12)",
    xlabel = "x", ylabel = "y", zlabel = "z")

# y < 0 panels culled -> cutaway view exposing the junction and contact faces
function draw_group!(ax, ids; cutaway = true)
    isempty(ids) && return
    pts = Point3f[]; faces = TriangleFace{Int}[]
    for i in ids
        c = panels[i].corners
        cutaway && all(q -> q[2] > 1e-9, c) && continue  # keep y <= 0 half
        n0 = length(pts)
        append!(pts, [Point3f(q...) for q in c])
        push!(faces, TriangleFace(n0 + 1, n0 + 2, n0 + 3))
        push!(faces, TriangleFace(n0 + 1, n0 + 3, n0 + 4))
    end
    isempty(faces) && return
    mesh!(ax, pts, faces; color = COLOR[pairkey(first(ids))], shading = NoShading,
          transparency = false)
    # panel borders
    seg = Point3f[]
    for i in ids
        c = panels[i].corners
        cutaway && all(q -> q[2] > 1e-9, c) && continue
        for (a, b) in ((1, 2), (2, 3), (3, 4), (4, 1))
            push!(seg, Point3f(c[a]...)); push!(seg, Point3f(c[b]...))
        end
    end
    linesegments!(ax, seg; color = (:black, 0.25), linewidth = 0.4)
end

for k in keys(census)
    draw_group!(ax, [i for i in eachindex(panels) if pairkey(i) == k])
end
scatter!(ax, [Point3f(SYS.src_center...)]; color = :red, markersize = 14, label = "source")
scatter!(ax, [Point3f(SYS.tgt_center...)]; color = :black, markersize = 10, marker = :diamond)

# legend via proxies
for (k, v) in sort(collect(census))
    scatter!(ax, [Point3f(NaN, NaN, NaN)]; color = COLOR[k], marker = :rect,
             markersize = 14, label = get(LABEL, k, string(k)))
end
axislegend(ax; position = :rt, framevisible = false, labelsize = 11)

save(joinpath(FIGS, "fig61_system_fig1_geometry.png"), fig; px_per_unit = 2)
save(joinpath(FIGS, "fig61_system_fig1_geometry.pdf"), fig)
println("\nwrote figs/fig61_system_fig1_geometry.{png,pdf}")
