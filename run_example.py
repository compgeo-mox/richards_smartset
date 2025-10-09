"""Run a short example using the Richards1D solver and save a figure."""

import numpy as np
from richards1d import Richards1D
import matplotlib.pyplot as plt


def main():
    # domain
    L = 1.0
    N = 101
    z = np.linspace(0, L, N)
    dt = 60.0  # 1 minute
    t_max = 3600.0  # 1 hour

    vg = dict(alpha=1.5, n=1.35, theta_s=0.43, theta_r=0.078, Ks=1e-5)
    solver = Richards1D(z, dt, vg)

    # initial: hydrostatic profile from bottom head -1 m to top -1 m
    h0 = -1.0

    # top: atmospheric head (0 m) - Dirichlet; bottom: fixed head -1 m
    bc_top = ("dirichlet", 0.0)
    bc_bot = ("dirichlet", -1.0)

    times, h_hist, theta_hist = solver.solve(h0, t_max, bc_top=bc_top, bc_bottom=bc_bot)

    # plot final head and theta
    final_h = h_hist[-1]
    final_theta = theta_hist[-1]

    fig, ax1 = plt.subplots()
    ax1.plot(final_h, z, label="h (m)")
    ax1.set_xlabel("Pressure head h (m)")
    ax1.set_ylabel("z (m)")
    ax1.invert_yaxis()
    ax1.legend(loc="upper left")

    ax2 = ax1.twiny()
    ax2.plot(final_theta, z, "r", label="theta")
    ax2.set_xlabel("theta (m3/m3)")
    ax2.legend(loc="upper right")

    plt.title(f"Richards 1D final profiles at t={times[-1]:.0f}s")
    plt.tight_layout()
    plt.savefig("richards_1d_final.png")
    print("Saved figure richards_1d_final.png")


if __name__ == "__main__":
    main()
