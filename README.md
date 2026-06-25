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

The Mathematica notebooks 3), 4) require Wolfram Mathematica. The Wolfram Language package 5) was developed in Wolfram Language 14.x. Check your installation with:

```bash
wolframscript -code 'Print[$Version]'
```

### 1) `Gauge_Variation.py`

This is the main file containing the infrastructure for automating the gauge variation algorithm set out in Section 3.2 of the paper. Note that the procedure is fully generalised and we do not use the specialised gauge variation operator (3.19) in the text - this is instead used in `Gauge invariance via differential operators.nb` as a semi-independent check. This file also generates the DDI expressions in accordance with equation (3.21) of the paper. That is, we compute all DDI forms generated from mulitplying the over-antisymmetrised expressions by additional factors of Yi, keeping all such expressions whose total number of derivatives (Yi's) lies between two and four. These are further analysed in `AdS_DDIs_Analysis.py`

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
4. Run `Gauge_Variation.py`.
5. The fully evaluated gauge variation in the {Y_i, Z_i} basis is stored in `full_gauge_variation_equation`.
6. The gauge variation reduced modulo the DDIs is stored in `final_gauge_variation_equation`.

Note that the gauge-variation procedure is intended to be fully general within the symbolic setup used here. The DDI reduction procedure, however, is specialised to the two- and three-derivative vertices considered in the paper.

The code contains a large number of helper functions for basic symbolic manipulation and replacement rules. These are mostly self-explanatory and are not described individually here. Instead, we list below the core functions that implement the main steps of the {Yi,Zi} basis algorithm, since these are the functions most relevant for reproducing and checking the results of the paper.

| Function | Role in the algorithm |
|---|---|
| `move_all_a_left(equation)` | Moves explicit auxiliary variables `a1`, `a2`, and `a3` to the left using the relevant commutation rules. Terms with leftmost auxiliary variables are then removed, corresponding to setting the auxiliary variables to zero after the vertex has acted on the generating functions. |
| `remove_traces(equation)` | Removes trace terms, i.e. terms containing contractions that vanish in the transverse-traceless sector. |
| `pull_all_non_canon_UP_left(equation)` | Rewrites non-canonical contractions of the form `Ui*Pi-1` into the canonical `Y_i` structures and divergence terms using integrations by parts and AdS commutator relations. |
| `pull_all_Divs_right(equation)` | Moves divergence terms of the form `U_i*P_i` to the right, where they are discarded using the transverse condition. |
| `pull_all_PiPjs_left(equation)` | Reorders non-canonical derivative contractions `P_i*P_j` with i$\neq$j and rewrites them in terms of d'Alembertians using integration by parts. |
| `pull_all_Bs_right(equation)` | Moves d'Alembertian's, i.e. terms of the form `Pi_a*Pi^a`, to the right so that the on-shell conditions can be imposed. |

These operations are combined in the wrapper function `perform_full_operation(equation)`.

A few example checks are included here to illustrate the use of these functions. Readers are encouraged to run additional tests if they wish to further verify the intermediate reductions and final outputs.

### Example checks:

Example 1:

```python
move_all_a_left('U1_w*P2^w*U1_x*P2^x*a1_a*P1^a')
```

returns

```python
'U1_w*P2^w*P2_a*P1^a + P2_a*P2^x*U1_x*P1^a'
```

Example 2:

```python
pull_all_non_canon_UP_left('U3_b*P3^b*U1_a*P3^a*U2_c*P3^c')
```

returns

```python
'(1/l^2)*U2^a*U1_a*U3_b*P3^b + -(1/l^2)*U2^b*U1_a*U3_b*P3^a + -d*(1/l^2)*U2_c*U1_a*U3^a*P3^c + (1/l^2)*U2_c*U1_a*U3^a*P3^c + -U2_c*U1_a*U3_b*P2^a*P3^b*P3^c + -U2_c*U1_a*U3_b*P1^a*P3^b*P3^c'
```

Example 3:

```python
pull_all_PiPjs_left('U2_c*U3^c*U2_a*P3^a*P3_b*P2^b')
```

returns

```python
'-(1/l^2)*U2_a*U2^a*U3_b*P2^b + (1/l^2)*U2_a*U2_b*U3^a*P2^b + (1/2)*U2_a*U3^c*U2_c*P1_b*P1^b*P3^a + -(1/2)*U2_a*U3^c*U2_c*P2_b*P2^b*P3^a + -(1/2)*U2_a*U3^c*U2_c*P3_b*P3^b*P3^a'
```

One can verify the accuracy of these outputs using the commutator relations listed in Appendix A of the paper.

### 2) `AdS_DDIs_Analysis.py`

This file contains supplementary analysis of the AdS3 DDIs generated by Gauge_Variation.py. It includes some intermediate checks and exploratory simplifications, retained for completeness. Its main purpose is to rewrite the DDIs into a compact differential-operator form closer to the presentation used in the paper.

To verify the output:

1. Ensure that Gauge_Variation.py is in the same directory.
2. Run `python AdS_DDIs_Analysis.py`
3. The simplified forms of the DDIs are stored in simplified_DDIs.

Since this file imports objects from Gauge_Variation.py, running it may also execute parts of the gauge-variation script.

### 3) `Gauge invariance via differential operators.nb`

This Mathematica notebook provides a semi-independent check of the gauge-invariance conditions using the specialised differential-operator form of the gauge variation given in equation (3.19) of the paper, as well as the differential-operator forms of the DDIs (B.1-B.10). Unlike Gauge_Variation.py, this notebook does not build the gauge variation by explicitly commuting the underlying Ui, Pi, and ai operators.

The notebook makes use of the following:

| Object/function     | Role                                                                                                       |
| ------------------- | ---------------------------------------------------------------------------------------------------------- |
| `dy1`, `dy2`, `dy3` | Derivative operators with respect to `y1`, `y2`, and `y3`.                                                 |
| `dz1`, `dz2`, `dz3` | Derivative operators with respect to `z1`, `z2`, and `z3`.                                                 |
| `AdSGaugeOp1`       | The specialised AdS gauge-variation operator corresponding to equation (3.19) of the paper.                |
| `V2flat`            | The flat-space part of the two-derivative vertex.                                                          |
| `V2ads`             | The lower-derivative AdS correction ansatz for the two-derivative vertex.                                  |
| `V2`                | The full two-derivative ansatz, including both `V2flat` and `V2ads`.                                       |
| `dV2`               | The gauge variation of the two-derivative ansatz.                                                          |
| `miRules`           | The mass-shell substitution rules for the `mi` placeholders.                                              |
| `nRules`            | The substitution rules relating the formal powers `ni` to the spin labels in the two-derivative sector.   |
| `dV23yPart`         | The three-derivative part of the two-derivative gauge variation.                                           |
| `dV21yPart`         | The one-derivative part of the two-derivative gauge variation.                                             |
| `ruleYcubedAll`     | The DDI reduction rules used for cubic expressions in the `yi`.                                           |
| `ruleYsqY`          | The DDI reduction rules used for the remaining two-derivative-sector terms.                                |
| `dV2DDIs`           | The two-derivative gauge variation after the relevant DDI substitutions have been prepared.                |
| `dV2DDIreduced`     | The final DDI-reduced form of the two-derivative gauge variation.                                          |
| `V3`                | The full three-derivative ansatz, including its lower-derivative AdS correction terms.                     |
| `dV3`               | The gauge variation of the three-derivative ansatz.                                                        |
| `pRules`            | The substitution rules relating the formal powers `pi` to the spin labels in the three-derivative sector. |
| `dV34yPart`         | The four-derivative part of the three-derivative gauge variation.                                          |
| `dV32ySimplif`      | The simplified two-derivative part of the three-derivative gauge variation.                                |
| `dV30ySimplif`      | The simplified zero-derivative part of the three-derivative gauge variation.                               |
| `ruleYsqYYAll`      | The DDI reduction rules used for the three-derivative-sector gauge variation.                              |
| `dV3DDIReduced`     | The final DDI-reduced form of the three-derivative gauge variation.                                        |

The notebook first applies `AdSGaugeOp1` to the two-derivative ansatz `V2`, producing `dV2`. After applying `miRules` and `nRules`, the result is separated into `dV23yPart` and `dV21yPart`, reduced using `ruleYcubedAll` and `ruleYsqY`, and stored in the final form `dV2DDIreduced`.

The same procedure is then repeated for the three-derivative ansatz `V3`. The gauge variation `dV3` is simplified using `miRules` and `pRules`, separated into `dV34yPart`, `dV32ySimplif`, and `dV30ySimplif`, and then reduced using `ruleYsqYYAll`. The final reduced expression is stored in `dV3DDIReduced`.

These reduced expressions are then used either to verify the lower-derivative coefficients or to solve for them directly. The resulting coefficients agree with those obtained in `Gauge_Variation.py`.

### 4) `DDI Check Via Forward Reduction Matching.nb`

This Mathematica notebook provides an alternative check of a two-derivative DDI used in the paper. The purpose is to track the terms that are dropped or modified when the over-antisymmetrised expression is reduced in the transverse-traceless, on-shell basis, and then to systematically reintroduce those terms in order to make the vanishing of the original DDI expression manifest.

In other words, the notebook checks that the compact DDI used in the paper is genuinely equivalent to the original over-antisymmetrised expression, up to the reductions allowed in the on-shell transverse-traceless sector. It does this by first recording the trace, divergence, box, and integration-by-parts terms removed or modified during the forward reduction, and then reversing these steps starting from the reduced DDI expression.

The notebook is organised around the following core objects and functions.

| Object/function(s) | Role |
|---|---|
| `antisymPart`, `antiExpr` | Define the original over-antisymmetrised expression and the DDI expression to be reduced. |
| `writeComponents` | Expands the expression into explicit component form, allowing the vanishing of the over-antisymmetrised expression to be checked directly. |
| `splitTraceTerms`, `traceDropped`, `ddiAfterTrace` | Separate and store trace terms, which are dropped in the transverse-traceless reduction. |
| `reduceAdjacentBoxPairsStep`, `boxChangedTermsStep1`, `ddiAfterBoxStep1` | Apply and record the adjacent-box/mass-shell reductions. |
| `reduceNonCanonicalUPStep`, `upChangedTermsStep1`, `ddiAfterUPStep1` | Apply and record the integration-by-parts reductions of non-canonical `U_i P_{i-1}` contractions. |
| `splitDivTerms`, `divDropped`, `ddiAfterDiv` | Separate and store divergence terms, which are dropped using the transverse condition. |
| `canonicalYZForm`, `ddiCanonicalYZ` | Rewrite the remaining expression into the compact canonical `{y[i], z[i]}` basis. |
| `ddiLHScomm`, `ddiRHScomm` | Define the compact reduced DDI expression and its AdS correction in commuting `{y[i], z[i]}` notation. |
| `ddiLHSop`, `ncExpand`, `removeNCExceptSameP` | Convert the compact DDI into an ordered operator expression while retaining the non-commutativity of same-field covariant derivatives. |
| `reverseReplacementData`, `restoredUPIBP`, `restoredBoxIBP` | Reintroduce the terms that were changed during the integration-by-parts and box-reduction steps. |
| `orderNCFixed` | Orders the remaining same-field derivative products and introduces explicit commutator terms where necessary. |
| `replaceCommP` | Replaces the derivative commutators by the corresponding AdS curvature terms. |
| `eliminateA` | Commutes to the left the explicit auxiliary-variable factors generated by the derivative commutator. |
| `DerivCommFullDDILHS` | Stores the fully restored and commutator-corrected DDI left-hand side. |
| `DDILHSminusRHS` | Stores the difference between the restored left-hand side and the compact AdS correction term. |
| `reduceTr`, `LHSminusRHSTrReduce` | Applies the final trace reduction. The final expression `LHSminusRHSTrReduce` should vanish. |

To run the check, open the notebook in Mathematica and evaluate all cells. The final result is the vanishing of `LHSminusRHSTrReduce`, which verifies that the differential-operator form of the DDI is equivalent to the original over-antisymmetrised expression after systematically restoring the terms dropped or modified during the forward reduction.


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

## Referencing

If you use this repository, please cite both the accompanying paper and this code repository. Until a DOI is available, the repository may be cited using the GitHub URL.

```bibtex
@software{KingPomfretCode2026,
  author       = {King, Freddie and Pomfret,Taylor},
  title        = {Accompanying code for the paper "Metric-like Cubic Vertices for Massless Bosonic Higher-Spin Fields in AdS3"},
  year         = {2026},
  version      = {1.0.0},
  url          = {https://github.com/fk260/ads3_higher_spin_cubic_vertices_gauge_variation_and_DDIs}
}
```

## Contact

Questions about the code, reproducibility, or the calculations in the accompanying paper can be directed to:

**Freddie King** `fk260@sussex.ac.uk`

**Taylor Pomfret** `pomfret@mpa-garching.mpg.de`

Please include enough detail to identify the relevant part of the repository, such as the file name, function name, notebook cell, or expression being checked.
