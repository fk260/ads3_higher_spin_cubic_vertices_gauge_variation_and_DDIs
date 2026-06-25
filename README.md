# Accompanying code for the paper "Metric-like Cubic Vertices for Massless Bosonic Higher-Spin Fields in AdS3".

This repository contains the symbolic code used in the paper *Metric-like Cubic Vertices for Massless Bosonic Higher-Spin Fields in AdS(_3)*. It provides two independent checks of the dimension-dependent identities alongside an alternative verification for one of the two-derivative DDIs (see Section 3.3 and Appendix B), together with two semi-independent checks of the gauge invariance of the two- and three-derivative cubic vertices (see Sections 3.4.1, 3.4.2 and Appendix C).

## Repository contents
```text
Gauge_Variation.py
AdS_DDIs_Analysis.py
Gauge invariance via differential operators.nb
DDI Check Via Forward Reduction Matching.nb
AdSDDIReducer.wl
README.md
LICENSE
```

## Files

### Requirements

The code in files 1) and 2) were written for Python 3. A recent Python 3 installation is recommended. Only packages from the Python standard library are used:

```python
re
time
collections
fractions
itertools
```

The Mathematica notebooks 3), 4) and Wolfram Language package 5) require Wolfram Mathematica.

### 1) Gauge_Variation.py

This is the main file containing the infrastructure for automating the gauge variation algorithm set out in Section 3.2 of the paper. Note that the procedure is fully generalised and we do not use the specialised gauge variation operator (3.19) in the text - this is instead used in a separate file as a semi-independent check.

Below is a dictionary for the translating the notation used to the paper to the notation used in the code:

| Code notation | Paper notation | Meaning |
|---|---|---|
| `P1`, `P2`, `P3` | $\nabla_1$, $\nabla_2$, $\nabla_3$ | Covariant derivatives acting on fields 1, 2, and 3. In the code, index contractions are written explicitly, e.g. `P1_a`, `P1^a`. |
| `U1`, `U2`, `U3` | $\partial_{a_1}$, $\partial_{a_2}$, $\partial_{a_3}$ | Derivatives with respect to the auxiliary variables of fields 1, 2, and 3. |
| `a1`, `a2`, `a3` | $a_1$, $a_2$, $a_3$ | Auxiliary variables used in the generating-function description of the fields and gauge parameters, see equation (). |
| `Y1` | $y_1 = \partial_{a_1}\cdot\nabla_2$ | Implemented in the code by replacing `Y1` with e.g. `U1_a*P2^a`. |
| `Y2` | $y_2 = \partial_{a_2}\cdot\nabla_3$ | Implemented in the code by replacing `Y2` with e.g. `U2_a*P3^a`. |
| `Y3` | $y_3 = \partial_{a_3}\cdot\nabla_1$ | Implemented in the code by replacing `Y3` with e.g. `U3_a*P1^a`. |
| `Z1` | $z_1 = \partial_{a_2}\cdot\partial_{a_3}$ | Implemented in the code by replacing `Z1` with e.g. `U2_a*U3^a`. |
| `Z2` | $z_2 = \partial_{a_3}\cdot\partial_{a_1}$ | Implemented in the code by replacing `Z2` with e.g. `U1_a*U3^a`; the contraction is symmetric. |
| `Z3` | $z_3 = \partial_{a_1}\cdot\partial_{a_2}$ | Implemented in the code by replacing `Z3` with e.g. `U1_a*U2^a`. |
| `s1`, `s2`, `s3` | $s_1$, $s_2$, $s_3$ | Spins of the three fields. |
| `l` | $\ell$ | AdS radius. Curvature corrections appear with factors of `1/l^2`. |
| `n1`, `n2`, `n3` | $n_1$, $n_2$, $n_3$ | Formal powers of `Z1`, `Z2`, `Z3` in the two-derivative vertex sector. In the final spin-dependent expressions these are substituted using the spin labels. |
| `p1`, `p2`, `p3` | $p_1$, $p_2$, $p_3$ | Formal powers of `Z1`, `Z2`, `Z3` used for the three-derivative vertex sector. |
| `D_z1`, `D_z2`, `D_z3` | $\partial_{z_1}$, $\partial_{z_2}$, $\partial_{z_3}$ | Differential-operator notation used when rewriting powers such as `n1`, `n2`, `n3` as `Z1*D_z1`, `Z2*D_z2`, `Z3*D_z3`. |
| `m1`, `m2`, `m3` | $m_1$, $m_2$, $m_3$ | AdS Mass-shell/d'Alembertian placeholders. |
| `A`, `B`, `C` | Ansatz coefficients | Coefficients of lower-derivative AdS correction terms in some test expressions. |



### 2) AdS_DDIs_Analysis.py

### 3) Gauge invariance via differential operators.nb

### 4) DDI Check Via Forward Reduction Matching.nb

### 5) AdSDDIReducer.wl
