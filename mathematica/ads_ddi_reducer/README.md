# Ordered AdS3 DDI Reducer

This folder contains a Wolfram Language implementation of an ordered-operator reducer for three-dimensional AdS dimension-dependent identities (DDIs) acting on three higher-spin fields.

The core reducer works with ordered noncommuting operators `U[i, idx]`, `P[i, idx]`, and curvature-generated `A[i, idx]`. It reduces raw antisymmetrised DDI expressions through the pipeline

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

as outlined.

This repository folder contains an AdS DDI reducer, DDI-output scripts for a user to run, and tests to inspect how each stage of the algorithm operates.

## Folder layout

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

## Requirements

- WolframScript or Mathematica/Wolfram Language.
- The scripts were developed against Wolfram Language 14.x.

Check your installation with:

```bash
wolframscript -code 'Print[$Version]'
```

## Running the reducer outputs

From this folder:

```bash
bash run_ddis_v1.sh
```

This computes `DDI1`--`DDI10` for `V = 1` using only the sThis package can be found in the folder `ads_ddi_reducer`, whose goal is to perform AdS reduction in Mathemtica. The folder itself contains a Wolfram Language implementation of an ordered-operator reducer for three-dimensional AdS dimension-dependent identities (DDIs) acting on three higher-spin fields.

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

## Folder layout

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

## Requirements

- WolframScript or Mathematica/Wolfram Language.
- The scripts were developed against Wolfram Language 14.x.

Check your installation with:

```bash
wolframscript -code 'Print[$Version]'
```

## Running the reducer outputs

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

## User terms and stage tests

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

## Main package entry points

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
