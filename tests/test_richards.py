import numpy as np
from richards1d import Richards1D


def test_solver_runs_and_theta_bounds():
    z = np.linspace(0, 0.5, 51)
    dt = 60.0
    vg = dict(alpha=1.5, n=1.35, theta_s=0.43, theta_r=0.078)
    solver = Richards1D(z, dt, vg, Ks=1e-5)
    times, h_hist, theta_hist = solver.solve(
        h0=-0.5, t_max=600.0, bc_top=("dirichlet", 0.0), bc_bottom=("dirichlet", -1.0)
    )
    # final theta should be within [theta_r, theta_s]
    final_theta = theta_hist[-1]
    assert np.all(final_theta >= vg["theta_r"] - 1e-8)
    assert np.all(final_theta <= vg["theta_s"] + 1e-8)
    # shape checks
    assert h_hist.shape[0] == len(times)
    assert h_hist.shape[1] == len(z)
