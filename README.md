# Richards 1D solver

This small package implements an implicit Picard solver for the 1D Richards equation with the van Genuchten soil model.

Files
- `richards1d.py`: solver implementation (Richards1D class)
- `run_example.py`: example runner that saves a figure `richards_1d_final.png`
- `requirements.txt`: minimal dependencies
- `tests/test_richards.py`: simple test to check the solver runs


Equations and numerical scheme
------------------------------

Continuous PDE (mixed form)

The solver implements the 1D Richards equation in mixed form:
$$
\frac{\partial \theta(h)}{\partial t} = \frac{\partial}{\partial z}\Bigl[ K(h)\Bigl(\frac{\partial h}{\partial z} + 1\Bigr) \Bigr],
$$

where
- $h(z,t)$ is the pressure head ($m$),
- $\theta(h)$ is volumetric water content ($m^3/m^3$),
- $K(h)$ is the unsaturated hydraulic conductivity ($m/s$).

Constitutive relations (van Genuchten — Mualem)

The code uses the standard van Genuchten retention curve combined with the Mualem relative permeability model. Below we present the formulae and a LaTeX-style parameter table that explains every symbol and its units/typical ranges.

Effective saturation and retention (van Genuchten):

$$
S_e(h) = \begin{cases}
  (1 + |\alpha h|^n)^{-m}, & h<0,\\
  1, & h\ge 0,
\end{cases}\qquad m = 1 - \frac{1}{n},
$$

and

$$
\theta(h) = \theta_r + S_e(h)\, (\theta_s - \theta_r).
$$

Relative permeability (Mualem-van Genuchten):

$$
k_r(S_e) = S_e^{\ell} \Bigl[1 - \bigl(1 - S_e^{1/m}\bigr)^{m}\Bigr]^2,\qquad K(h)=K_s\,k_r(S_e(h)).
$$

List of symbols used.
<div align="center">

| Symbol | Meaning | Units / typical range |
|---|---:|---|
| $\alpha$ | van Genuchten inverse scale parameter (related to air-entry pressure) | $m^{-1}$ (e.g. 0.01–10) |
| $n$ | van Genuchten shape parameter (controls retention curve steepness) | dimensionless (>1, e.g. 1.1–3) |
| $m$ | auxiliary parameter, $m=1-1/n$ | dimensionless |
| $\theta_s$ | saturated volumetric water content (porosity) |$ m^3 / m^3$ (e.g. 0.3–0.5) |
| $\theta_r$ | residual volumetric water content | $m^{3}/m^{3}$ (e.g. 0.01–0.1) |
| $S_e$ | effective saturation, $S_e=(\theta-\theta_r)/(\theta_s-\theta_r)$ | dimensionless (0–1) |
| $K_s$ | saturated hydraulic conductivity | $m/s$ (e.g. 1e-8–1e-2) |
| $\ell$ | Mualem pore-connectivity / tortuosity parameter | dimensionless (often $\approx0.5$) |
| $k_r$ | relative hydraulic conductivity (function of $S_e$) | dimensionless (0–1) |

</div>

Notes:
- $\alpha$ is roughly the inverse of an entry pressure (large $\alpha$ $\,\rightarrow$ lower entry suction).
- $n$ controls the curvature/steepness of the retention curve; values close to 1 give gentler curves, large $n$ give steeper transitions.
- $\ell$ (ell) is often set to 0.5 in the literature (Mualem-Freundlich mixing parameter) but can be calibrated.
- Typical units and ranges above are rough guidelines; please consult soil-specific databases or laboratory fits for exact values.

Numerical discretisation used in this code
-----------------------------------------

- Time: backward Euler (implicit).
- Nonlinearity: Picard iteration — at each Picard iterate coefficients $C(h)=d\theta/dh$ and $K(h)$ are evaluated and a linear system for $h^{n+1}$ is solved.
- Space: second-order central finite differences on a uniform grid ($\Delta z$). Interface conductivity $K_{i+1/2}$ is computed by harmonic mean of neighbouring $K$ values.

Discrete algebraic form (interior node $i$)

The linear equation assembled for interior node $i$ (using $C_i$ and $K_{i\pm1/2}$ evaluated at the Picard iterate) is:
$$
\Bigl(\frac{C_i}{\Delta t} + \frac{K_{i-1/2}+K_{i+1/2}}{\Delta z^2}\Bigr) h_i^{n+1}
- \frac{K_{i-1/2}}{\Delta z^2} h_{i-1}^{n+1}
- \frac{K_{i+1/2}}{\Delta z^2} h_{i+1}^{n+1}
\;=\;
\frac{C_i}{\Delta t} h_i^{n} + \frac{K_{i-1/2}-K_{i+1/2}}{\Delta z},
$$
Boundary conditions
-------------------
- Dirichlet: prescribe head (enforced directly in the linear system).
- Neumann: prescribe flux q (code assumes positive downward). Treated via ghost-node/flux representation.

Notes
-----
- Code assumes uniform grid.
- C(h) is computed numerically (central finite difference) for robustness; an analytical derivative can be added for speed.

Quick start

1. Create a virtual environment and install dependencies:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

2. Run the example:

```bash
python run_example.py
```

3. Run tests:

```bash
pytest -q
```
