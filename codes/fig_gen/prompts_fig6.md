Write a JCP-style subsection describing a component-level numerical experiment that validates the post-refinement strategy for repeated evaluation of layer-potential contributions.

Context:
After solving the boundary integral equation, the single-layer density sigma is fixed on the Nyström discretization of the dielectric interface Gamma. For many target densities rho_j supported in a nearby volume region, we must repeatedly evaluate the layer-potential contribution

    V_j,layer
    =
    (1 / 4 pi)
    int rho_j(x)
        int_Gamma sigma(y) / |x - y| dS_y
    dx.

This evaluation stage is distinct from the BIE solve. The purpose of post-refinement is not to change the solved unknowns sigma or improve the BIE discretization. Instead, for a fixed sigma, the source panels near the target volume are dyadically refined only for evaluation. The density sigma is interpolated from the original Nyström grid to the refined panels. This reduces the number of genuinely near source-target interactions that require expensive direct or adaptive quadrature, and allows most interactions to be handled efficiently by FMM.

Post-refinement rule:
For each original source panel P, refine P if

    d(P, supp rho_j) < 5 h_P,

or, for a batch of target densities supported in a common region Omega_t,

    d(P, Omega_t) < 5 h_P.

Refine dyadically until the refined panel width satisfies

    |P'| <= h0,

where h0 is the evaluation-stage threshold. Interpolate sigma to the refined panels using tensor-product barycentric interpolation.

Experiment design:
Use a fixed dielectric geometry, such as a cube or slab, and first solve the BIE to obtain sigma on the original adaptive Nyström discretization. Then freeze sigma. Choose a target region Omega_t located near one interface, for example a small box or Gaussian-support region close to one face. Construct a family of smooth target densities rho_j supported in the same Omega_t, representing repeated orbital-pair density evaluations in the same material region.

The key point is that all target densities share approximately the same support Omega_t. Therefore, the post-refined source discretization can be constructed once and reused across many target densities. This is essential for the amortization experiment.

Construct a 1 x 3 figure:

Panel (a): post-refinement geometry.
Show the original interface panels, the target support Omega_t, and the source panels selected for post-refinement. Then show the refined child panels near Omega_t. Color panels by refinement level or distinguish original and post-refined panels. Emphasize that the post-refinement is applied only for evaluation after sigma has already been solved.

Panel (b): accuracy versus post-refinement threshold.
For one representative target density rho, compute

    E_V(h0)
    =
    | V_h0 - V_ref | / | V_ref |,

where V_h0 is the layer contribution evaluated using post-refinement threshold h0. Plot E_V versus h0, or versus the number of post-refined evaluation nodes N_eval. The reference V_ref should be computed using a much finer post-refinement level together with high-accuracy direct or adaptive quadrature for the remaining near interactions. The expected behavior is that the error decreases as h0 is reduced and eventually reaches a plateau determined by the target-volume quadrature, the reference accuracy, or the original BIE discretization error.

Panel (c): repeated-target acceleration.
Plot total evaluation time as a function of the number of target densities N_t. Compare at least two strategies:

1. No post-refinement:
   Use the original panel discretization of sigma. Far interactions are evaluated by FMM, but many near source-target interactions remain and must be evaluated by expensive direct or adaptive quadrature.

2. With post-refinement:
   Construct the post-refined source discretization once for the common target region Omega_t. Interpolate sigma once onto the refined panels. Then evaluate each target density using FMM plus a much smaller residual near-interaction correction.

The post-refinement curve should include the one-time setup cost. The expected result is a break-even point: for small N_t, post-refinement may not be advantageous because of the setup overhead; for larger N_t, the smaller per-target evaluation cost dominates and the post-refined strategy becomes faster.

Optional additional quantities for panel (c) or a small table:
Report the number of near source-target interactions before and after post-refinement, the number of evaluation nodes after post-refinement, the setup time for interpolation/refinement, and the per-target evaluation time.

Important implementation details:
1. Use the same fixed solved density sigma for all post-refinement tests.
2. Do not rerun the BIE solve when h0 changes. Only the evaluation representation changes.
3. The post-refined source density should be obtained by interpolation from the solved Nyström density, not by resolving the BIE on the refined panels.
4. The target densities in the timing experiment should share a common support region Omega_t; otherwise the post-refinement setup cannot be cleanly amortized.
5. The reference for panel (b) should be more accurate than the plotted errors and should not reuse the same h0-level approximation.
6. Separate the one-time post-refinement setup cost from the per-target evaluation cost, but include both in the total-time curve.

Manuscript points to emphasize:
1. Post-refinement is an evaluation-stage device, not a modification of the BIE discretization.
2. It resolves near source-target geometry after sigma has been computed, so it does not increase the number of unknowns in the GMRES solve.
3. Refining near the target support converts most interactions into the far-field regime and reduces the expensive residual near-interaction set.
4. The cost is amortized when many target densities share the same spatial support region.
5. The accuracy is controlled by h0 until other error sources dominate.

Suggested figure caption:
“Post-refinement for repeated layer-potential evaluation. (a) Original interface panels and evaluation-stage post-refined panels near a common target support region Omega_t. The BIE density sigma is solved on the original discretization and interpolated to the refined panels only after the solve. (b) Relative error in the layer contribution V_layer versus post-refinement threshold h0 or evaluation nodes N_eval. (c) Total evaluation time versus the number of target densities N_t. The post-refined strategy includes its one-time setup cost and becomes advantageous once the setup is amortized over sufficiently many targets.”