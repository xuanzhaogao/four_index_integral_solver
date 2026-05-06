"""Generate the TKM partition schematic (Fig. 4).

Shows the source support and its 5h-neighborhood as rectangles on the xOy
plane, and a Gaussian source density rho(x) plotted as a translucent surface
above the source region.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import art3d  # noqa: F401  (registers 3d projection)

OUT = Path(__file__).with_name("tkm_partition.png")


def gaussian_sum(X, Y, centers, sigma=0.18, amps=None):
    if amps is None:
        amps = [1.0] * len(centers)
    Z = np.zeros_like(X)
    for (cx, cy), a in zip(centers, amps):
        Z += a * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * sigma**2))
    return Z


def main():
    # source support: [src_lo, src_hi]^2; 5h-neighborhood: [nbh_lo, nbh_hi]^2
    src_lo, src_hi = -0.5, 0.5
    h5 = 0.22  # represents 5h
    nbh_lo, nbh_hi = src_lo - h5, src_hi + h5

    fig = plt.figure(figsize=(8.0, 6.4))
    ax = fig.add_subplot(111, projection="3d")

    # --- xOy plane rectangles ---
    # 5h-neighborhood: dashed blue
    nbh = Rectangle(
        (nbh_lo, nbh_lo),
        nbh_hi - nbh_lo,
        nbh_hi - nbh_lo,
        linewidth=1.6,
        linestyle="--",
        edgecolor="#1f4e79",
        facecolor="none",
        zorder=1,
    )
    ax.add_patch(nbh)
    art3d.pathpatch_2d_to_3d(nbh, z=0, zdir="z")

    # source support: filled orange (outline emphasized)
    src = Rectangle(
        (src_lo, src_lo),
        src_hi - src_lo,
        src_hi - src_lo,
        linewidth=2.0,
        edgecolor="#b35900",
        facecolor="#ffd9a8",
        alpha=0.55,
        zorder=2,
    )
    ax.add_patch(src)
    art3d.pathpatch_2d_to_3d(src, z=0, zdir="z")

    # --- Gaussian source density above the plane ---
    grid = np.linspace(src_lo, src_hi, 120)
    X, Y = np.meshgrid(grid, grid)
    centers = [(-0.18, -0.10), (0.18, 0.14)]
    amps = [1.0, 0.7]
    Z = gaussian_sum(X, Y, centers, sigma=0.13, amps=amps)
    # taper to zero at the source-box boundary so support is clearly inside
    taper = (
        np.cos(np.pi * X / (2 * src_hi)) ** 2
        * np.cos(np.pi * Y / (2 * src_hi)) ** 2
    )
    Z = Z * taper
    Z = 0.55 * Z / Z.max()  # rescale for visual height

    ax.plot_surface(
        X,
        Y,
        Z,
        rstride=2,
        cstride=2,
        cmap="viridis",
        alpha=0.88,
        linewidth=0,
        antialiased=True,
        zorder=3,
    )

    # --- annotations ---
    # 5h gap on the LEFT side of the figure (along x at y = src_hi midline)
    y_left = src_hi + (nbh_hi - src_hi) / 2
    ax.plot(
        [nbh_lo, src_lo],
        [y_left, y_left],
        [0.0, 0.0],
        color="#1f4e79",
        lw=1.5,
    )
    for xv in (nbh_lo, src_lo):
        ax.plot([xv, xv], [y_left - 0.025, y_left + 0.025], [0, 0],
                color="#1f4e79", lw=1.2)
    ax.text(
        (nbh_lo + src_lo) / 2,
        y_left + 0.06,
        0.0,
        r"$5h$",
        color="#1f4e79",
        fontsize=13,
        ha="center",
        va="bottom",
    )

    # L: diagonal of the source bounding box (drawn on the plane in dark orange)
    ax.plot(
        [src_lo, src_hi],
        [src_lo, src_hi],
        [0.001, 0.001],
        color="#7a3a00",
        lw=1.4,
        linestyle=(0, (4, 2)),
    )
    ax.text(
        0.05,
        -0.18,
        0.0,
        r"$L$",
        color="#7a3a00",
        fontsize=14,
        ha="left",
        va="top",
    )

    # supp rho label on the RIGHT side of the source square
    ax.text(
        src_hi + 0.03,
        src_lo + 0.05,
        0.0,
        r"$\mathrm{supp}\,\rho$",
        color="#7a3a00",
        fontsize=13,
        ha="left",
        va="bottom",
    )

    # rho(x) label hovering above the surface peak
    pk = np.unravel_index(np.argmax(Z), Z.shape)
    ax.text(
        X[pk],
        Y[pk],
        Z.max() + 0.06,
        r"$\rho(x)$",
        fontsize=14,
        ha="center",
        va="bottom",
    )

    # --- axes cosmetics ---
    span = nbh_hi - nbh_lo + 0.15
    ax.set_xlim(nbh_lo - 0.08, nbh_hi + 0.18)
    ax.set_ylim(nbh_lo - 0.12, nbh_hi + 0.08)
    ax.set_zlim(0.0, max(0.75, Z.max() * 1.4))
    ax.set_xlabel("$x$", fontsize=12)
    ax.set_ylabel("$y$", fontsize=12)
    ax.set_zlabel("$z$", fontsize=12)

    ax.view_init(elev=22, azim=-62)
    ax.set_box_aspect((span, span, 0.85 * span))

    # lighten panes
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor((0.85, 0.85, 0.85, 1.0))
        axis._axinfo["grid"]["color"] = (0.9, 0.9, 0.9, 1.0)

    plt.tight_layout()
    fig.savefig(OUT, dpi=220, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
