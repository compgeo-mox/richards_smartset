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

\[\]
\\frac{\\partial \\theta(h)}{\\partial t} = \\frac{\\partial}{\\partial z}\\Bigl[ K(h)\\Bigl(\\frac{\\partial h}{\\partial z} + 1\\Bigr) \\Bigr],
\[\]

where
- h(z,t) is the pressure head (m),
- \\theta(h) is volumetric water content (m^3/m^3),
- K(h) is the unsaturated hydraulic conductivity (m/s).

Constitutive relations (van Genuchten + Mualem)

- Effective saturation
	\[\]
	S_e(h) = \\begin{cases}
		(1 + |\\alpha h|^n)^{-m}, & h<0,\\\\
		1, & h\\ge 0, 
	\\end{cases}
	\quad m = 1 - \\frac{1}{n}.
	\[\]
- Retention (van Genuchten)
	\\theta(h) = \\theta_r + S_e(h)\\, (\\theta_s - \\theta_r).
- Relative permeability (Mualem-van Genuchten)
	\[\]
	k_r(S_e) = S_e^{\\ell} \\Bigl[1 - \\bigl(1 - S_e^{1/m}\\bigr)^{m}\\Bigr]^2,\\qquad K(h)=K_s\\,k_r(S_e(h)).
	\[\]

Numerical discretisation used in this code
-----------------------------------------

- Time: backward Euler (implicit).
- Nonlinearity: Picard iteration — at each Picard iterate coefficients C(h)=d\\theta/dh and K(h) are evaluated and a linear system for h^{n+1} is solved.
- Space: second-order central finite differences on a uniform grid (\\Delta z). Interface conductivity K_{i+1/2} is computed by harmonic mean of neighbouring K values.

Discrete algebraic form (interior node i)

The linear equation assembled for interior node i (using C_i and K_{i\pm1/2} evaluated at the Picard iterate) is:

\[\]
\\Bigl(\\frac{C_i}{\\Delta t} + \\frac{K_{i-1/2}+K_{i+1/2}}{\\Delta z^2}\\Bigr) h_i^{n+1}
- \\frac{K_{i-1/2}}{\\Delta z^2} h_{i-1}^{n+1}
- \\frac{K_{i+1/2}}{\\Delta z^2} h_{i+1}^{n+1}
\\;=\\;
\\frac{C_i}{\\Delta t} h_i^{n} + \\frac{K_{i-1/2}-K_{i+1/2}}{\\Delta z},
\[\]

Boundary conditions
-------------------
- Dirichlet: prescribe head (enforced directly in the linear system).
- Neumann: prescribe flux q (code assumes positive downward). Treated via ghost-node/flux representation.

Notes
-----
- Code assumes uniform grid. If you need non-uniform grids or a higher-order spatial scheme, I can extend it.
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
