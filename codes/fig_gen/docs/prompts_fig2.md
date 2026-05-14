Write a JCP-style subsection describing a numerical experiment that validates the adaptive refinement of the boundary-integral right-hand side for a dielectric-interface Poisson problem.

Context:
We solve the dielectric-interface Poisson problem using the second-kind boundary integral equation

    (1/2) gamma^{-1} sigma + D^T sigma = - partial_n u_inc  on Gamma,

where Gamma is the union of dielectric interfaces and

    f = - partial_n u_inc

is the right-hand side. When a source density is close to an interface, f can vary rapidly over a small region of Gamma. A uniform panelization is inefficient because it must globally refine the entire interface to resolve this localized variation. The purpose of this experiment is to demonstrate that RHS-driven adaptive dyadic refinement resolves f efficiently.

Experiment design:
Use a simple dielectric box or slab embedded in vacuum. The geometry should be simple enough that the refinement pattern is easy to interpret. For example, use a rectangular dielectric slab or cube with exterior permittivity epsilon_0 = 1 and interior permittivity epsilon_d = 10. Place a localized smooth source density near one face of the interface, for example a normalized Gaussian

    rho(x) = C exp(-|x - x_0|^2 / (2 s^2)),

with x_0 close to the selected face and s small compared with the box size. The source should be close enough that the induced right-hand side f = - partial_n u_inc has a sharp localized feature on the nearby face, but not so close that it is numerically singular.

Use a fixed tensor-product Gauss-Legendre order p, for example p = 8 or p = 10. Starting from a coarse panelization of Gamma, apply RHS-driven dyadic refinement. On each panel P, approximate f by a tensor-product polynomial interpolant I_P f of degree p - 1 in each local coordinate. Refine P until the local interpolation error satisfies

    || f - I_P f ||_{L^2(P)} < epsilon,

or, preferably, a relative local/global criterion such as

    || f - I_P f ||_{L^2(P)} < epsilon ||f||_{L^2(Gamma)}

with an appropriate normalization.

The interpolation error must be measured on an oversampled check rule, not on the original p x p Gauss-Legendre interpolation nodes. For example, use a q x q Gauss-Legendre rule with q = 2p or q = 3p on each panel. This avoids falsely reporting zero error at the interpolation nodes.

Construct a figure with two or three panels, depending on available space:

Panel (a):
Show the final adaptive dyadic panelization of the interface. Color each panel by refinement level. Mark or indicate the location of the source center x_0. The refinement should concentrate near the projection of the source onto the closest face, while remaining coarse far from the localized feature.

Panel (b):
Plot the magnitude of the right-hand side |f| or log10(|f|) on the same interface face, together with the local interpolation error |f - I_P f| evaluated on an oversampled check grid. If using a single panel, show these as two side-by-side surface/heat maps; if space is limited, show only the interpolation error heat map. This panel should demonstrate that the refined mesh tracks the localized variation of f.

Panel (c), optional but recommended:
Plot the global relative interpolation error

    E_f =
    ( sum_{P subset Gamma} || f - I_P f ||_{L^2(P)}^2 )^{1/2}
    /
    || f ||_{L^2(Gamma)}

as a function of the total number of boundary unknowns N. Compare RHS-driven adaptive refinement with uniform dyadic refinement. For the adaptive curve, vary the tolerance epsilon. For the uniform curve, vary the global refinement level. The expected result is that adaptive refinement reaches a given E_f with far fewer unknowns because it refines only near the localized feature.

Important implementation details:
1. Keep p fixed while varying the refinement tolerance epsilon or uniform refinement level. Otherwise the experiment mixes p-convergence with mesh refinement.
2. Use the same source density and geometry for all curves.
3. Compute E_f on an independent oversampled check grid.
4. If the source is very close to the interface, verify that the TKM or direct evaluation of u_inc and partial_n u_inc is accurate enough so that the plotted error reflects interpolation/refinement error rather than error in evaluating f.
5. This figure should test only the RHS-resolution component. Do not include the solved layer density sigma or final four-index integral error in this figure.

Manuscript points to emphasize:
1. The right-hand side f is smooth but highly localized when the source density approaches an interface.
2. Uniform refinement is inefficient because it resolves the entire interface at the scale required only near the source-proximal region.
3. The adaptive criterion based on the oversampled interpolation error refines only where f is under-resolved.
4. The observed decay of E_f with N demonstrates that the adaptive panelization gives an efficient discretization of the incident Neumann data.
5. This experiment is independent of edge singularities in sigma; it validates the RHS-driven part of the refinement strategy.

Suggested figure caption:
“RHS-driven adaptive refinement for a localized source near a dielectric interface. (a) Final dyadic panelization colored by refinement level; refinement concentrates near the projection of the source onto the closest interface. (b) Magnitude of the incident Neumann data |f| = |partial_n u_inc| and/or the oversampled interpolation error |f - I_P f| on the refined interface. (c) Relative interpolation error E_f versus the number of boundary unknowns N. Adaptive refinement reaches a given interpolation error with substantially fewer unknowns than uniform refinement.”