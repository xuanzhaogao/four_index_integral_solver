Write a JCP-style subsection describing an end-to-end convergence experiment for the complete screened four-index integral under edge refinement.

Context:
Earlier component-level experiments have already validated the RHS-driven adaptive refinement, the edge singularity structure of the BIE density, near-field panel-panel correction, the TKM volume solver, and post-refinement for repeated evaluation. This experiment is part of the final Numerical Results section and should therefore focus on the complete pipeline and the final target quantity.

We compute the screened Coulomb integral

    V = int rho_tar(x) phi_src(x) dx,

where phi_src is obtained by solving the dielectric-interface Poisson problem using the adaptive BIE solver. The goal is to test how the final integral converges as the edge-refinement depth is increased, and how the GMRES iteration count changes for the same discretizations.

Experiment design:
Fix one representative dielectric geometry, such as a dielectric cube or finite slab embedded in vacuum. Fix the source density rho_src and target density rho_tar. These may be smooth Gaussian densities or representative orbital-pair densities from the intended application. Fix all numerical parameters except the Gauss-Legendre order p and the edge-refinement depth r.

The fixed parameters should include:
- RHS-adaptive refinement tolerance;
- volume grid resolution;
- near-field threshold, for example the 5h rule;
- post-refinement threshold h0;
- GMRES relative tolerance;
- dielectric parameters;
- source and target densities.

Vary only:
- the tensor-product Gauss-Legendre order p on each panel, for example p = 2, 4, 6;
- the edge-refinement depth r.

The edge-refinement depth r controls the minimum edge/corner panel size. For example,

    h_ec(r) = 2^{-r} h_ec(0),

or equivalently r dyadic refinement levels are applied to panels near edges and corners.

Reference value:
Use an over-resolved computation as the reference,

    V_ref = V_{p_max, r_max},

with the largest p and r used in the study, or with an even finer discretization. Verify the reliability of V_ref by checking that it changes negligibly when compared with V_{p_max, r_max - 1} and/or V_{p_max - 2, r_max}. Mention this stability check in the text. The reference error should be smaller than the smallest plotted error by at least one order of magnitude if feasible.

Construct a 1 x 2 figure:

Panel (a): final integral convergence.
Plot the relative integral error

    E_V(p,r) = | V_{p,r} - V_ref | / | V_ref |

as a function of the edge-refinement depth r. Use a logarithmic scale for the vertical axis. Plot one curve for each Gauss-Legendre order p. The expected behavior is that, at small r, the error is dominated by under-resolution of the edge and corner singularities, so increasing p alone gives limited improvement. As r increases, the singular density is better resolved and E_V decreases until other fixed error sources dominate.

Panel (b): GMRES iteration count.
For the same discretizations, plot the number of GMRES iterations

    N_iter(p,r)

required to reach the prescribed solver tolerance as a function of r. Use the same set of p values and the same markers/colors as in panel (a). This panel characterizes the linear-solver cost associated with increasing p and edge-refinement depth. If the iteration count grows only mildly, state that the edge-refined second-kind discretization remains practical for the tested geometry and dielectric contrast. If the growth is substantial, report it directly and interpret it as part of the cost of resolving edge-localized density features.

Important implementation details:
1. Do not vary the RHS tolerance, volume grid, near-field threshold, post-refinement threshold, GMRES tolerance, dielectric geometry, or densities across the curves.
2. This experiment should measure convergence of the complete integral V, not merely the BIE residual or layer density.
3. The BIE residual does not need to be plotted here if it has already been shown in the component-level edge-refinement diagnostic.
4. Report the boundary unknown count

       N_Gamma(p,r) = p^2 M(r),

   where M(r) is the number of leaf panels at edge-refinement depth r. Since p changes the number of unknowns by p^2, the text or caption should state the range of N_Gamma values corresponding to the plotted curves. Timing and memory scaling with N_Gamma should be shown separately in a later figure.
5. Use a GMRES tolerance smaller than the target discretization errors so that solver error does not contaminate the convergence curves.
6. If the curves plateau, explain that the plateau is caused by other fixed error sources, such as the volume grid, near-field threshold, post-refinement threshold, reference accuracy, or floating-point precision.

Manuscript points to emphasize:
1. This is an end-to-end convergence test for the final screened integral, not a component-level diagnostic.
2. Varying r isolates the effect of resolving edge and corner singularities, because all other numerical parameters are fixed.
3. Comparing several p values shows the interaction between high-order panel quadrature and edge-local mesh refinement.
4. The GMRES iteration plot quantifies the linear-solver cost of the refined BIE discretization.
5. The corresponding CPU time and memory growth with N_Gamma are addressed separately in the scaling experiment.

Suggested figure caption:
“End-to-end convergence of the screened integral under edge refinement. The dielectric geometry, source and target densities, RHS-adaptive tolerance, volume grid, near-field threshold, post-refinement threshold, dielectric parameters, and GMRES tolerance are fixed. (a) Relative error E_V = |V_{p,r} - V_ref| / |V_ref| versus edge-refinement depth r for several Gauss-Legendre orders p. (b) Number of GMRES iterations required to reach the prescribed solver tolerance for the same discretizations. The corresponding boundary unknown counts N_Gamma = p^2 M(r) are reported in the text; timing and memory scaling with N_Gamma are shown separately.”