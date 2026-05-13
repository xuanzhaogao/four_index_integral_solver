Write a JCP-style subsection describing a component-level numerical experiment that validates the near-field quadrature correction for close panel-panel interactions in the boundary integral discretization.

Context:
We discretize the dielectric-interface boundary integral equation

    (1/2) gamma^{-1} sigma + D^T sigma = - partial_n u_inc  on Gamma,

using tensor-product Gauss-Legendre Nyström quadrature on flat rectangular panels. For two well-separated panels P and Q, the standard p x p Gauss-Legendre rule accurately approximates the adjoint double-layer interaction block

    (D^T_{P Q} sigma_Q)(x)
    =
    (1 / 4 pi) int_Q ((x - y) dot n_x) / |x - y|^3 sigma_Q(y) dS_y,
    x in P.

However, when P and Q are geometrically close, the kernel is nearly singular on the target panel, and the standard panelwise rule may lose accuracy. The near-field correction replaces the standard block by an upsampled block for panel pairs satisfying

    d(P,Q) <= c min(h_P, h_Q),

where h_P and h_Q denote the local Gauss-Legendre grid spacings and c is typically 5.

The goal of this experiment is to isolate this quadrature issue at the panel-pair level and verify that the upsampled near-field correction restores accuracy without increasing the number of Nyström unknowns.

Experiment design:
Use two flat rectangular panels P and Q. Let both panels have side length 1. Take P to be the target panel in the plane z = 0 with normal n_P = e_z, and Q to be the source panel in the parallel plane z = d. For example,

    P = [-1/2, 1/2]^2 x {0},
    Q = [-1/2, 1/2]^2 x {d}.

Use a smooth source density on Q, such as

    sigma_Q(y_1,y_2) = exp(-20 ((y_1 - 0.1)^2 + (y_2 + 0.15)^2)),

or another nontrivial smooth test density. Evaluate the adjoint double-layer potential at the p x p Gauss-Legendre nodes on P.

Important orientation detail:
Choose the target normal n_P so that the kernel does not vanish. For the parallel-panel setup above, n_P = e_z gives

    (x - y) dot n_P = -d,

so the interaction is nonzero. Do not use coplanar panels, because the adjoint double-layer kernel vanishes identically for interactions within the same flat plane.

For a fixed quadrature order p, for example p = 8 or p = 10, vary the separation distance d. Express the separation using the dimensionless ratio d/h, where h is the local node spacing on the panels, for example h = 1/p or the minimum nearest-neighbor distance among the Gauss-Legendre nodes mapped to the physical panel. Test values of d/h ranging from well below 1 to about 10 or 20.

Compute three approximations:
1. Standard p x p Gauss-Legendre panel quadrature.
2. Near-field corrected quadrature using source-panel upsampling when d <= 5h.
3. Near-field corrected quadrature using source-panel upsampling when d <= 8h.

For the upsampled rule, interpolate sigma_Q from the original p x p Gauss-Legendre nodes to a finer tensor-product Gauss-Legendre grid on Q, or evaluate the analytic test density directly on the upsampled grid for this diagnostic. The first option is closer to the actual algorithm; the second option isolates pure quadrature error. State clearly which one is used.

Error metric:
For each separation d, compute the relative block-action error

    E_near(d)
    =
    || D^T_{P Q} sigma_Q - D^{T,ref}_{P Q} sigma_Q ||_2
    /
    || D^{T,ref}_{P Q} sigma_Q ||_2,

where the norm is over the p x p target nodes on P. Plot E_near versus d/h.

Construct a 1 x 3 figure:

Panel (a):
Show a schematic of the two-panel geometry. Label P, Q, the separation distance d, the target normal n_P, and the local grid spacing h. Indicate that target nodes lie on P and quadrature nodes lie on Q.

Panel (b):
Plot E_near(d) versus d/h on a log scale. Compare standard Gauss-Legendre quadrature, 5h near-field correction, and 8h near-field correction. The expected behavior is that the uncorrected rule deteriorates when d/h is small, while the corrected rule maintains the target accuracy over the near-field regime. The 8h rule should be at least as accurate as the 5h rule, at the cost of treating more panel pairs as near-field in full geometries.

Panel (c):
Plot the cost/accuracy tradeoff for varying near-field threshold c in

    d(P,Q) <= c h.

For the same panel-pair test, this can be represented by the largest error over the near-field interval d/h <= c, together with the relative cost of the upsampled block. Alternatively, if using a representative multi-panel slab geometry, plot the fraction of panel pairs classified as near-field as a function of c. The goal is to justify the choice c = 5 as a practical compromise between accuracy and overhead.

Manuscript points to emphasize:
1. The experiment isolates the quadrature error for close panel-panel interactions and does not test the full BIE solver.
2. Standard Gauss-Legendre quadrature loses accuracy when the source panel lies within a few local node spacings of the target panel.
3. Upsampling the source panel for near interactions restores the desired accuracy while leaving the number of BIE unknowns unchanged.
4. The threshold c = 5 is a heuristic accuracy/cost compromise; increasing it to c = 8 improves robustness for higher accuracy requirements but increases the number of corrected interactions in full geometries.
5. The reference values are computed by over-resolved quadrature or adaptive cubature, not by FMM.

Suggested figure caption:
“Near-field quadrature correction for close panel-panel interactions. (a) Two parallel rectangular panels separated by distance d; the adjoint double-layer potential is evaluated at Gauss-Legendre nodes on the target panel P due to a smooth density on the source panel Q. (b) Relative block-action error E_near as a function of d/h. Standard Gauss-Legendre quadrature deteriorates when the panels are separated by only a few local grid spacings, while the upsampled near-field correction maintains the target accuracy. (c) Accuracy/cost tradeoff for the near-field threshold c in d(P,Q) <= c h, illustrating why c = 5 is used as the default and c = 8 may be used for stricter tolerances.”

Reference solution:
Compute the reference value for each target node using HCubature.jl applied directly to the source panel Q. For each target point x in P, evaluate

    (D^{T,ref}_{P Q} sigma_Q)(x)
    =
    (1 / 4 pi) int_Q ((x - y) dot n_P) / |x - y|^3 sigma_Q(y) dS_y

by two-dimensional adaptive integration over the local coordinates of Q.

For the parallel-panel test with

    Q = [-1/2, 1/2]^2 x {d},

write y = (eta_1, eta_2, d), eta_1, eta_2 in [-1/2, 1/2]. Then for each target node x,

    I_ref(x)
    =
    (1 / 4 pi) int_{-1/2}^{1/2} int_{-1/2}^{1/2}
    ((x - y(eta)) dot n_P) / |x - y(eta)|^3
    sigma_Q(eta_1, eta_2)
    d eta_1 d eta_2.

Use HCubature.jl with absolute and relative tolerances at least one to two orders of magnitude smaller than the smallest error reported in the plot, for example rtol = 1e-12 and atol = 1e-14, if feasible. The reference should be recomputed independently for each target node and each separation d. Do not use FMM or the same Gauss-Legendre rule as the reference.

Mention in the manuscript that HCubature.jl is used only to generate reference values for this small panel-pair diagnostic and is not part of the production BIE solver.