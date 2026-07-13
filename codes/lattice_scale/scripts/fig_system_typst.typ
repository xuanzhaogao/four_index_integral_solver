// System-geometry schematic (Typst + cetz), faithful to the benchmark:
//   left  — overview: two 270-Å Si|SiO2 cubes with the small (~43-Å) epsilon=10 host slab patch
//           sitting on the junction; an arrow points to the zoomed slab.
//   right — top-down view of the slab: the REAL 198 benchmark graphene atoms (square window),
//           A/B sublattices, split at Δx=0 by the buried Si|SiO2 junction.
// Compile: typst compile --root <lattice_scale> scripts/fig_system_typst.typ figs/fig_system_typst.pdf
#import "@preview/cetz:0.3.4"
#import "../figs/atoms_data.typ": ATOMS

#set page(width: auto, height: auto, margin: 6pt)
#set text(size: 10pt)

#let c_si = rgb("#8fa9cc"); #let c_ox = rgb("#e0cfa6"); #let c_slab = rgb("#8fc79a")
#let c_A = rgb("#4477AA"); #let c_B = rgb("#EE7733")
#let pr(x, y, z) = (x + 0.40 * y, z + 0.32 * y)     // oblique projection (x right, y back, z up)

#cetz.canvas(length: 15pt, {
  import cetz.draw: *

  // ============ left: overview (cubes + small slab patch) ============
  let box3(x0, x1, y0, y1, z0, z1, col, sw) = {
    line(pr(x0,y1,z1), pr(x1,y1,z1), pr(x1,y0,z1), pr(x0,y0,z1), close: true, fill: col.lighten(12%), stroke: sw + black)
    line(pr(x1,y0,z0), pr(x1,y1,z0), pr(x1,y1,z1), pr(x1,y0,z1), close: true, fill: col.darken(14%), stroke: sw + black)
    line(pr(x0,y0,z0), pr(x1,y0,z0), pr(x1,y0,z1), pr(x0,y0,z1), close: true, fill: col, stroke: sw + black)
  }
  // 270-Å cubes (side 4 units); slab is ~43 Å ≈ 0.64 units wide, 9 Å ≈ 0.13 thick → small patch
  box3(0, 4, 0, 4, 0, 4, c_si, 0.7pt)            // Si cube  (Δx<0)
  box3(4, 8, 0, 4, 0, 4, c_ox, 0.7pt)            // SiO2 cube (Δx>0)
  box3(3.68, 4.32, 1.7, 2.36, 4, 4.15, c_slab, 0.5pt)   // ε=10 slab patch on the junction
  line(pr(4,0,0), pr(4,0,4), stroke: (dash: "dashed", paint: black, thickness: 1pt))   // buried junction
  rect((pr(3.3,2.36,4.15).at(0), pr(3.3,1.7,4.0).at(1)),   // red callout around the slab patch
       (pr(4.5,2.36,4.15).at(0), pr(4.5,2.36,4.3).at(1)),
       stroke: (dash: "dashed", paint: rgb("#c0392b"), thickness: 0.9pt))
  content(pr(0.5,4,1.8), text(fill: black)[Si])
  content(pr(4.5,4,1.8), text(fill: black)[SiO#sub[2]])

  // ============ right: top-down view of the slab ============
  let CX = 15.5; let CY = 2.6; let s = 0.34; let YMID = 10.3
  let td(dx, y) = (CX + s * dx, CY + s * (y - YMID))
  // slab footprint (top-down rectangle), light green
  line(td(-12, -2.2), td(12, -2.2), td(12, 22.8), td(-12, 22.8), close: true,
       fill: c_slab.transparentize(72%), stroke: 0.8pt + c_slab.darken(20%))
  // C–C bonds + atoms
  for i in range(0, ATOMS.len()) {
    for j in range(i + 1, ATOMS.len()) {
      let dx = ATOMS.at(i).at(0) - ATOMS.at(j).at(0)
      let dy = ATOMS.at(i).at(1) - ATOMS.at(j).at(1)
      if dx * dx + dy * dy < 2.6 {
        line(td(ATOMS.at(i).at(0), ATOMS.at(i).at(1)), td(ATOMS.at(j).at(0), ATOMS.at(j).at(1)), stroke: 0.9pt + gray)
      }
    }
  }
  line(td(0, -2.2), td(0, 22.8), stroke: (dash: "dashed", paint: black, thickness: 1pt))   // junction Δx=0
  for a in ATOMS {
    circle(td(a.at(0), a.at(1)), radius: 0.42 * s, fill: if a.at(2) == 1 { c_A } else { c_B }, stroke: 0.35pt + black)
  }

  // ============ arrow: slab patch → top-down zoom ============
  line((5.4, 4.8), (11.0, 4.0), mark: (end: ">", fill: rgb("#c0392b")), stroke: 1.1pt + rgb("#c0392b"))
})
