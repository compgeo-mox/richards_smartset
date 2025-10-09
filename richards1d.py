"""
1D Richards equation solver (Picard iteration, implicit time stepping)
Implements van Genuchten soil model for retention and relative permeability.

Usage: import Richards1D from this module and call solve(...)

This implementation favors clarity and robustness over raw performance.
"""

import numpy as np


class Richards1D:
    def __init__(self, z, dt, van_genuchten_params):
        """
        z : 1D array of node coordinates (m)
        dt: time step (s)
        van_genuchten_params: dict with keys alpha (1/m), n, theta_s, theta_r
        Ks : saturated hydraulic conductivity (m/s)
        l  : pore-connectivity parameter (commonly 0.5)
        """
        self.z = np.asarray(z)
        self.N = len(z)
        self.dz = np.diff(self.z)[0]

        self.Se_min = 1e-12
        self.dt = dt

        self.Ks = van_genuchten_params.get("Ks", 1e-5)
        self.l = van_genuchten_params.get("l", 0.5)
        self.alpha = van_genuchten_params.get("alpha")
        self.n = van_genuchten_params.get("n")
        self.m = 1.0 - 1.0 / self.n
        self.theta_s = van_genuchten_params.get("theta_s")
        self.theta_r = van_genuchten_params.get("theta_r")

    def theta(self, h):
        """Volumetric water content as function of pressure head h (m).
        h can be array-like.
        Uses van Genuchten model (h negative in unsaturated zone).
        For h >= 0 we take theta = theta_s (saturated).
        """
        h = np.asarray(h)
        theta = np.empty_like(h, dtype=float)
        # unsaturated
        mask = h < 0
        if np.any(mask):
            ah = np.abs(self.alpha * h[mask])
            Se = (1.0 + ah**self.n) ** (-self.m)
            theta[mask] = self.theta_r + Se * (self.theta_s - self.theta_r)
        # saturated or near-surface
        if np.any(~mask):
            theta[~mask] = self.theta_s
        return theta

    def K(self, h):
        """Hydraulic conductivity (m/s) as function of h using Mualem-van Genuchten"""
        h = np.asarray(h)
        K = np.empty_like(h, dtype=float)
        mask = h < 0
        if np.any(mask):
            # reuse the theta() function to compute effective saturation Se
            theta_mask = self.theta(h[mask])
            Se = (theta_mask - self.theta_r) / (self.theta_s - self.theta_r)
            # Avoid negative or zero
            Se = np.clip(Se, self.Se_min, 1.0)
            kr = Se**self.l * (1 - (1 - Se ** (1 / self.m)) ** self.m) ** 2
            K[mask] = self.Ks * kr
        if np.any(~mask):
            K[~mask] = self.Ks
        return K

    def C(self, h):
        """Specific moisture capacity dtheta/dh. Using central finite difference for stability/robustness."""
        h = np.asarray(h, dtype=float)
        eps = 1e-6
        h_plus = h + eps
        h_minus = h - eps
        th_plus = self.theta(h_plus)
        th_minus = self.theta(h_minus)
        C = (th_plus - th_minus) / (2 * eps)
        # enforce non-negative
        C = np.maximum(C, 1e-12)
        return C

    def solve(
        self,
        h0,
        t_max,
        bc_top=("dirichlet", 0.0),
        bc_bottom=("dirichlet", -1.0),
        tol=1e-4,
        max_iter=50,
    ):
        """
        Solve Richards equation from t=0 to t=t_max (s).
        h0 : initial head array of length N or scalar (uniform)
            dz = self.dz
            For Neumann flux, positive flux means downward flux (m/s).
        tol: picard tolerance for max head change (m)
        max_iter: max picard iterations per time step

        Returns (z, h, theta, times)
        - h is array shape (nt, N)
        - theta same shape
        - times array of time points
        """
        if np.isscalar(h0):
            h_old = np.ones(self.N) * h0
        else:
            h_old = np.asarray(h0, dtype=float).copy()
        times = [0.0]
        h_hist = [h_old.copy()]
        theta_hist = [self.theta(h_old)]
        t = 0.0
        nt = int(np.ceil(t_max / self.dt))
        # prepare dz array per interface in uniform case
        dz = self.dz
        dz_arr = np.ones(self.N - 1) * dz

        for step in range(nt):
            t += self.dt
            # Picard iterations
            h_iter = h_old.copy()
            for it in range(max_iter):
                # compute K at nodes and K at interfaces (harmonic mean)
                K_nodes = self.K(h_iter)
                K_half = np.empty(self.N - 1)
                # harmonic mean to avoid zero conductivity issues
                for i in range(self.N - 1):
                    a = K_nodes[i]
                    b = K_nodes[i + 1]
                    if a <= 0 or b <= 0:
                        K_half[i] = 0.0
                    else:
                        K_half[i] = 2.0 * a * b / (a + b)

                C_nodes = self.C(h_iter)

                # Assemble tridiagonal system A h_new = rhs
                A_diag = np.zeros(self.N)
                A_lower = np.zeros(self.N - 1)  # sub-diagonal (i-1)
                A_upper = np.zeros(self.N - 1)  # super-diagonal (i+1)
                rhs = np.zeros(self.N)

                # interior nodes
                for i in range(1, self.N - 1):
                    kL = K_half[i - 1]
                    kR = K_half[i]
                    A_diag[i] = C_nodes[i] / self.dt + (kL + kR) / (
                        dz_arr[i - 1] * dz_arr[i - 1]
                    )
                    A_lower[i - 1] = -kL / (dz_arr[i - 1] * dz_arr[i - 1])
                    A_upper[i] = -kR / (dz_arr[i - 1] * dz_arr[i - 1])
                    rhs[i] = C_nodes[i] * h_old[i] / self.dt + (1.0 / dz_arr[i - 1]) * (
                        kL - kR
                    )

                # Top boundary i=0
                if bc_top[0] == "dirichlet":
                    A_diag[0] = 1.0
                    rhs[0] = bc_top[1]
                elif bc_top[0] == "neumann":
                    # specified flux q_top (positive downward). Implement as ghost node flux.
                    q_top = bc_top[1]
                    # for node 0: use interface flux at top = q_top
                    kR = K_half[0]
                    A_diag[0] = C_nodes[0] / self.dt + kR / (dz_arr[0] * dz_arr[0])
                    A_upper[0] = -kR / (dz_arr[0] * dz_arr[0])
                    rhs[0] = C_nodes[0] * h_old[0] / self.dt + (1.0 / dz_arr[0]) * (
                        q_top - kR
                    )
                else:
                    raise ValueError("Unknown bc_top")

                # Bottom boundary i=N-1
                if bc_bottom[0] == "dirichlet":
                    A_diag[-1] = 1.0
                    rhs[-1] = bc_bottom[1]
                elif bc_bottom[0] == "neumann":
                    q_bot = bc_bottom[1]
                    kL = K_half[-1]
                    A_diag[-1] = C_nodes[-1] / self.dt + kL / (dz_arr[-1] * dz_arr[-1])
                    A_lower[-1] = -kL / (dz_arr[-1] * dz_arr[-1])
                    rhs[-1] = C_nodes[-1] * h_old[-1] / self.dt + (1.0 / dz_arr[-1]) * (
                        kL - q_bot
                    )
                else:
                    raise ValueError("Unknown bc_bottom")

                # Solve tridiagonal system using Thomas algorithm
                h_new = self._solve_tridiag(A_lower, A_diag, A_upper, rhs)

                diff = np.max(np.abs(h_new - h_iter))
                h_iter = h_new
                if diff < tol:
                    break

            # done time step
            h_old = h_iter.copy()
            h_hist.append(h_old.copy())
            theta_hist.append(self.theta(h_old))
            times.append(t)

        return np.array(times), np.array(h_hist), np.array(theta_hist)

    @staticmethod
    def _solve_tridiag(a, b, c, d):
        """Solve tridiagonal system with vectors a (n-1), b (n), c (n-1), d (n)
        Returns x of length n.
        """
        n = len(b)
        ac = a.copy()
        bc = b.copy()
        cc = c.copy()
        dc = d.copy()
        # forward
        for i in range(1, n):
            if bc[i - 1] == 0:
                raise ZeroDivisionError("Zero pivot in tridiagonal")
            m = ac[i - 1] / bc[i - 1]
            bc[i] = bc[i] - m * cc[i - 1]
            dc[i] = dc[i] - m * dc[i - 1]
        # back substitution
        x = np.zeros(n)
        x[-1] = dc[-1] / bc[-1]
        for i in range(n - 2, -1, -1):
            x[i] = (dc[i] - cc[i] * x[i + 1]) / bc[i]
        return x
