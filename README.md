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

### 1) `Gauge_Variation.py`

This is the main file containing the infrastructure for automating the gauge variation algorithm set out in Section 3.2 of the paper. Note that the procedure is fully generalised and we do not use the specialised gauge variation operator (3.19) in the text - this is instead used in a separate file as a semi-independent check. This file also generates the DDIs in accordance with equation (3.21) of the paper. That is, we compute all DDI forms generated from mulitplying the over-antisymmetrised expressions by additional factors of Yi, keeping all such expressions whose total number of derivatives (Yi's) lies between two and four. These are further analysed in 

Below is a dictionary for the translating the notation used to the paper to the notation used in the code:

| Code notation | Paper notation | Meaning |
|---|---|---|
| `P1`, `P2`, `P3` | $\nabla_1$, $\nabla_2$, $\nabla_3$ | Covariant derivatives acting on fields 1, 2, and 3. In the code, index contractions are written explicitly, e.g. `P1_a`, `P1^a`. |
| `U1`, `U2`, `U3` | $\partial_{a_1}$, $\partial_{a_2}$, $\partial_{a_3}$ | Derivatives with respect to the auxiliary variables of fields 1, 2, and 3. |
| `a1`, `a2`, `a3` | $a_1$, $a_2$, $a_3$ | Auxiliary variables used in the generating-function description of the fields and gauge parameters, see e.g. equation (2.12). |
| `Y1` | $y_1 = \partial_{a_1}\cdot\nabla_2$ | Implemented in the code by replacing `Y1` with e.g. `U1_a*P2^a`. |
| `Y2` | $y_2 = \partial_{a_2}\cdot\nabla_3$ | Implemented in the code by replacing `Y2` with e.g. `U2_a*P3^a`. |
| `Y3` | $y_3 = \partial_{a_3}\cdot\nabla_1$ | Implemented in the code by replacing `Y3` with e.g. `U3_a*P1^a`. |
| `Z1` | $z_1 = \partial_{a_2}\cdot\partial_{a_3}$ | Implemented in the code by replacing `Z1` with e.g. `U2_a*U3^a`. |
| `Z2` | $z_2 = \partial_{a_3}\cdot\partial_{a_1}$ | Implemented in the code by replacing `Z2` with e.g. `U1_a*U3^a` |
| `Z3` | $z_3 = \partial_{a_1}\cdot\partial_{a_2}$ | Implemented in the code by replacing `Z3` with e.g. `U1_a*U2^a`. |
| `s1`, `s2`, `s3` | $s_1$, $s_2$, $s_3$ | Spins of the three fields. |
| `l` | $\ell$ | AdS radius. Curvature corrections appear with factors of `1/l^2`. |
| `n1`, `n2`, `n3` | $n_1$, $n_2$, $n_3$ | Formal powers of `Z1`, `Z2`, `Z3` in the two-derivative vertex sector. In the final spin-dependent expressions these are substituted using the spin labels. |
| `p1`, `p2`, `p3` | $p_1$, $p_2$, $p_3$ | Formal powers of `Z1`, `Z2`, `Z3` used for the three-derivative vertex sector. |
| `D_z1`, `D_z2`, `D_z3` | $\partial_{z_1}$, $\partial_{z_2}$, $\partial_{z_3}$ | Differential-operator notation used when rewriting powers such as `n1`, `n2`, `n3` as `Z1*D_z1`, `Z2*D_z2`, `Z3*D_z3`. |
| `m1`, `m2`, `m3` | $m_1$, $m_2$, $m_3$ | AdS Mass-shell/d'Alembertian placeholders. |
| `A`, `B`, `C` | $\alpha$, $\alpha_1$, $\alpha_2$, $\alpha_3$ | Ansatz coefficients of lower-derivative AdS correction terms in some test expressions. |

To run a particular check:

1. Uncomment the relevant gauge-variation commutator expression assigned to `commutator`.
2. Uncomment the corresponding lower-derivative ansatz assigned to `z_term_comm`.
3. Set the `no_deriv` argument of `fully_process_three_deriv_gauge_variation_eq` (currently defined around line 6260) to the number of derivatives, i.e. the number of Y_i factors, in the expression assigned to commutator.
4. The fully evaluated gauge variation in the {Y_i, Z_i} basis is stored in `full_gauge_variation_equation`.
5. The gauge variation reduced modulo the DDIs is stored in `final_gauge_variation_equation`.

Note that the gauge-variation procedure is intended to be fully general within the symbolic setup used here. The DDI reduction procedure, however, is specialised to the two- and three-derivative vertices considered in the paper.

The code also contains several helper functions for basic symbolic manipulation. These are mostly self-explanatory and are not described individually here. Instead, we list below the core functions that implement the main steps of the gauge-variation and DDI-reduction algorithms, since these are the functions most relevant for reproducing and checking the results of the paper.


| Function | Role in the reduction |
|---|---|
| `move_all_a_left(equation)` | Moves explicit auxiliary variables `a1`, `a2`, and `a3` to the left using the relevant commutation rules. Terms with leftmost auxiliary variables are then removed, corresponding to setting the auxiliary variables to zero after differentiation. |
| `remove_traces(equation)` | Removes trace terms, i.e. terms containing contractions that vanish in the transverse-traceless sector. |
| `pull_all_non_canon_UP_left(equation)` | Rewrites non-canonical contractions of the form `U_i*P_j`, with `i \neq j`, into the canonical `Y_i` structures and divergence terms using integrations by parts and commutator relations. |
| `pull_all_Divs_right(equation)` | Moves divergence terms of the form `U_i*P_i` to the right, where they are discarded using the transverse condition. |
| `pull_all_PiPjs_left(equation)` | Reorders non-canonical derivative contractions `P_i*P_j` and rewrites them in terms of the standard derivative-contraction structures used in the reduction. |
| `pull_all_Bs_right(equation)` | Moves d'Alembertian/mass-shell terms, represented by `B_i`-type placeholders, to the right so that the on-shell equations can be imposed. |


### 2) `AdS_DDIs_Analysis.py`

This file contains supplementary analysis of the AdS3 DDIs generated by Gauge_Variation.py. It includes some intermediate checks and exploratory simplifications, retained for completeness. Its main purpose is to rewrite the DDIs into a compact differential-operator form closer to the presentation used in the paper.

To verify the output:

1. Ensure that Gauge_Variation.py is in the same directory.
2. Run `python AdS_DDIs_Analysis.py`
3. The simplified forms of the DDIs are stored in simplified_DDIs.

Since this file imports objects from Gauge_Variation.py, running it may also execute parts of the gauge-variation script.

### 3) `Gauge invariance via differential operators.nb`

### 4) `DDI Check Via Forward Reduction Matching.nb`

### 5) `AdSDDIReducer.wl`

This package can be found in the folder `ads_ddi_reducer`, whose goal is to perform AdS reduction in Mathemtica. The folder itself contains a Wolfram Language implementation of an ordered-operator reducer for three-dimensional AdS dimension-dependent identities (DDIs) acting on three higher-spin fields.

The core reducer works with ordered noncommuting derivative operators `P[i, idx]`, `U[i, idx]`, and curvature-generated `A[i, idx]`. It reduces raw antisymmetrised DDI expressions through the pipeline

```text
A cleanup
trace removal
Y' reductions
A cleanup
Div reductions
A cleanup
P.P reductions
A cleanup
Box reductions
A cleanup
trace removal
repeat until canonical
```

as illustrated in the paper.

This repository folder contains an AdS DDI reducer, `ads_ddi_reducer`; DDI-output scripts for a user to run the reducer locally and reproduce the results seen in the paper, and also a series of tests/examples to illustrate how each stage of the algorithm works.

### Folder layout

```text
ads_ddi_reducer/
  package/
    AdSDDIReducer.wl

  ddis_v1/
    run_ddis_v1.wl
    run_ddis_v1.sh

  ddis_vz/
    run_ddis_vz.wl
    run_ddis_vz.sh

  terms_and_stage_tests/
    test_terms.wl
    test_terms.nb
    run_terms_and_stage_tests.sh

  helpers/
    run_ddis_v1.sh
    run_ddis_vz.sh
    run_terms_and_stage_tests.sh

  run_ddis_v1.sh
  run_ddis_vz.sh
  run_terms_and_stage_tests.sh
```

### Requirements

- WolframScript or Mathematica/Wolfram Language.
- The scripts were developed against Wolfram Language 14.x.

Check your installation with:

```bash
wolframscript -code 'Print[$Version]'
```

### Running the reducer outputs

From this folder:

```bash
bash run_ddis_v1.sh
```

This computes the 10 DDIs for `V = 1` using the reducer pipeline. Outputs are written to the folder:

```text
ddi_v1_outputs/
```

For the vertex

```text
V(z) = z1^n1 z2^n2 z3^n3
```

run:

```bash
bash run_ddis_vz.sh
```

This uses symbolic `ZPow[i, ni]` implementation. Outputs are written to:

```text
ddi_vz_outputs/
```

The reducer uses the following `ZPow` when performing clean-up of the A_i operators:

```text
ZPow[j,n] A_i -> A_i ZPow[j,n] + n ZPow[j,n-1] partner-U
```

only when `z_j` contains `U_i`. This is a compressed Leibniz rule for symbolic powers of `z_j`.

For cyclic outputs over fields `i = 1,2,3`, run:

```bash
RUN_CYCLIC_DDIS=True bash run_ddis_v1.sh
RUN_CYCLIC_SYMBOLIC=True bash run_ddis_vz.sh
```

### User terms and stage tests

Run:

```bash
bash run_terms_and_stage_tests.sh
```

Outputs are written to:

```text
term_stage_outputs/
```

This evaluates indexed input words such as

```text
U3^dU3_dU2^cU1^aP1_aP1_cP1_bP2^b
```

and runs stage-level tests for A-cleanup, trace removal, Y' movement, divergences, P.P movement, boxes, grouped stages, full-cycle reduction, symbolic `ZPow` A-cleanup orientation, and Lichnerowicz commutator checks.

You can also open the notebook:

```text
terms_and_stage_tests/test_terms.nb
```

Example interactive use:

```wl
Get["package/AdSDDIReducer.wl"];

opts = {
  "MaxCycles" -> 30,
  "MaxSubsteps" -> 6000,
  "Lambda" -> λAdS,
  "Dim" -> 3,
  "Mass" -> (m[#] &)
};

TermReduce["U3^dU3_dU2^cU1^aP1_aP1_cP1_bP2^b", Sequence @@ opts]
TermPipeline["U3^dU2_dU2^cU1^aP1_aP1_cP1_bP2^b", Sequence @@ opts]
```

### Main package entry points

```wl
DDIRawByKey["DDI1", i, {0,0,0}]
DDIRawByKey["DDI1", i, SymbolicZPowers[]]
ReduceAdS[raw]
ToYZPolynomial[reduced]
NonCanonicalReport[reduced]
CanonicalQ[reduced]

TermOps["U3^dU2_d"]
TermReduce["U3^dU2_d", opts]
TermPipeline["U3^dU2_d", opts]

StageReduce["ACleanup", word, opts]
StagePipelineRun[{"StartCleanup", "YpBlock", "DivBlock"}, word, opts]
ReducerStageTestSuite[opts]
```

`CanonicalQ` follows Wolfram naming conventions for Boolean predicates (`NumberQ`, `MatrixQ`, etc.). It means the result has no unresolved `A`, `Comm`, `RawPComm`, `Pmarked`, traces, divergences, boxes, noncanonical `Y'`, or mixed `P.P` structures, and projects only to terminal `y[i]` and `z[i]` contractions.
