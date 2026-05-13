Write a JCP-style subsection describing a numerical experiment on edge singularities and edge-local refinement for the boundary integral equation.

Context:
We solve a dielectric-interface Poisson problem using a second-kind boundary integral equation for the single-layer density sigma,

    (1/2) gamma^{-1} sigma + D^T sigma = - partial_n u_inc  on Gamma,

where Gamma is the union of dielectric interfaces. For polyhedral dielectric geometries, sigma develops singular behavior near edges and corners even when the incident field is smooth. The goal of this experiment is to demonstrate that edge-local dyadic refinement resolves the singular density sufficiently well for the BIE residual to converge.

Experiment design:
Use a single dielectric cube embedded in vacuum. Let the cube be Omega_d = [-1,1]^3, with permittivity epsilon_d and exterior permittivity epsilon_0 = 1. Choose a moderate dielectric contrast, for example epsilon_d / epsilon_0 = 10 or 20. Drive the system by a distant point source located outside the cube, for example x_s = (0,0,4), so that partial_n u_inc is smooth on the cube faces. This choice is deliberate: any non-smooth behavior in sigma near edges and corners is caused by the geometry, not by a localized or under-resolved right-hand side.

Use a fixed Gauss-Legendre order p, for example p = 8 or p = 10. Compare a sequence of edge-local dyadic refinements. The refinement is applied only near edges and corners until the minimum edge-panel width reaches h_min or, equivalently, until a prescribed edge-refinement depth ell is reached. Keep all other numerical parameters fixed. Solve the BIE to a tolerance well below the plotted discretization error.

Construct a 1 x 3 figure:

Panel (a):
Plot the computed surface density on one face of the cube, for example the top face z = 1. Plot either log10(|sigma|) or a signed logarithmic transform of sigma. The purpose is to show concentration and apparent algebraic growth of |sigma| near the edges and corners, even though the incident field is smooth. The plot should make the four edges of the face visibly more singular than the interior.

Panel (b):
Choose a line on the same face that approaches one edge normally. For example, on the top face z = 1, fix y = 0 and let

    r = 1 - x,   x = 1 - r,   y = 0,   z = 1,

so that r is the distance to the edge x = 1. For several edge-refinement depths ell, plot |sigma_ell(r)| versus r on a log-log scale. Interpolate sigma from the panel Gauss-Legendre nodes to the line samples using tensor-product barycentric interpolation on each panel. The curves should show self-convergence away from the unresolved innermost region and increasingly resolve the near-edge singular profile. Optionally include an inset showing the line-profile self-convergence error

    E_sigma,line(ell)
    =
    ( sum_m |sigma_ell(r_m) - sigma_ref(r_m)|^2 w_m )^{1/2}
    /
    ( sum_m |sigma_ref(r_m)|^2 w_m )^{1/2},

where sigma_ref is obtained from the deepest edge refinement. When computing this error, avoid comparing at distances r much smaller than the mesh scale of the coarser refinement level; restrict to a resolved interval such as r >= c h_min,ell.

Panel (c):
Plot the relative BIE residual

    R_BIE =
    || (1/2) gamma^{-1} sigma + D^T sigma + partial_n u_inc ||_{L^2(Gamma)}
    /
    || partial_n u_inc ||_{L^2(Gamma)}

as a function of the number of boundary unknowns N. Evaluate this residual on an independent oversampled check grid, not on the Nyström nodes used to solve the linear system. For example, use a q x q Gauss-Legendre check rule on each leaf panel, with q = 2p or q = 3p, and interpolate sigma to the check nodes before evaluating the residual. This is important: otherwise the residual mostly reflects the GMRES residual of the discrete system rather than the discretization error.

In panel (c), compare edge-local refinement against uniform refinement. The expected result is that edge-local refinement reduces R_BIE much more efficiently with respect to N, whereas uniform refinement requires many more degrees of freedom to achieve the same residual.

Manuscript points to emphasize:
1. The point source is placed far from the cube so that the incident Neumann data is smooth on Gamma.
2. The observed growth of |sigma| near the edges and corners is therefore a geometric edge singularity.
3. Edge refinement does not regularize or remove the singularity. It resolves the singular density sufficiently well that its contribution to the boundary integral equation is captured.
4. The residual R_BIE is evaluated on an independent oversampled check grid, so it measures discretization error rather than the linear solver residual.
5. The comparison with uniform refinement demonstrates that local dyadic refinement near edges and corners is essential for efficiency.

Suggested figure caption:
“Edge singularity and edge-refinement convergence for a dielectric cube driven by a distant point source. (a) Surface plot of log10(|sigma|) on the top face, showing concentration near edges and corners despite a smooth incident field. (b) Edge-normal profiles |sigma(r)|, with r = dist(x, edge), for increasing edge-refinement depth, demonstrating self-convergence of the resolved singular profile. (c) Relative BIE residual R_BIE, evaluated on an independent oversampled check grid, versus the number of boundary unknowns N. Edge-local refinement achieves a given residual with substantially fewer unknowns than uniform refinement.”