(* ::Package:: *)
(* User-term reducer examples and stage-by-stage tests. *)
If[$InputFileName =!= "", SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
Get[FileNameJoin[{Directory[], "package", "AdSDDIReducer.wl"}]];

opts = {"MaxCycles" -> 30, "MaxSubsteps" -> 6000, "Lambda" -> \[Lambda]AdS, "Dim" -> 3, "Mass" -> (m[#] &)};
outDir = FileNameJoin[{Directory[], "term_stage_outputs"}];
If[DirectoryQ[outDir], DeleteDirectory[outDir, DeleteContents -> True]]; CreateDirectory[outDir];
write[name_, expr_] := Export[FileNameJoin[{outDir, name}], ToString[expr, InputForm], "Text"];
polyOf[result_] := If[AssociationQ[result] && KeyExistsQ[result, "Polynomial"], result["Polynomial"], $Failed[Polynomial]];
reportOf[result_] := If[AssociationQ[result] && KeyExistsQ[result, "NonCanonicalReport"], result["NonCanonicalReport"], $Failed[Report]];

Print["Package version: ", $AdSDDIReducerVersion];
Print["Running user indexed terms and reducer-stage tests..."];

ClearAll[runTerm];
runTerm[name_String, word_String] := Module[{ops, pipe, base},
  base = name;
  Print["  term: ", name, " = ", word];
  ops = Quiet@Check[TermOps[word], $Failed];
  write[base <> "_parsed_ops.txt", ops];
  pipe = If[ops === $Failed, $Failed, Quiet@Check[TermPipeline[word, Sequence @@ opts], $Failed]];
  write[base <> "_pipeline.wl", pipe];
  write[base <> "_polynomial.txt", polyOf[pipe]];
  write[base <> "_report.txt", reportOf[pipe]];
  <|"Name" -> name, "Word" -> word, "ParsedOps" -> ops, "Pipeline" -> pipe, "Polynomial" -> polyOf[pipe], "Report" -> reportOf[pipe]|>
];

terms = <|
  "requested_trace_word_P2b" -> "U3^dU3_dU2^cU1^aP1_aP1_cP1_bP2^b",
  "older_trace_word_P1b" -> "U3^dU3_dU2^cU1_aP1_aP1_cP1_bP1^b",
  "trace_free_sibling" -> "U3^dU2_dU2^cU1^aP1_aP1_cP1_bP2^b",
  "box_mixed_probe" -> "U2^cU1^aP1_aP1_cP1_bP2^b",
  "yp_div_pp_probe" -> "U3_aP2^aU1_bP1^bP1_cP2^c",
  "a_cleanup_probe" -> "U1^aA1_aU2^bU3_b",
  "pp_tail_probe" -> "P1_aP1_bP1^bU2^a",
  "div_box_z_probe" -> "U1^aP1_aP1_bP1^bU2^cU3_c"
|>;
termResults = AssociationMap[runTerm[#, terms[#]] &, Keys[terms]];
write["indexed_term_results.wl", termResults];

(* Symbolic ZPow/A-cleanup orientation tests.  These check that A_i moves left through a ZPow block to its left, and does not act on a ZPow block to its right. *)
ClearAll[zpowLeftExpr, zpowWrongSideExpr, zpowExplicitExpr, zpowCleanupTests];
zpowLeftExpr[power_] := {Term[1, {ZPow[2, power], A[3, Cov["a"]], P[2, Con["a"]]}]};
zpowWrongSideExpr[power_] := {Term[1, {A[3, Cov["a"]], ZPow[2, power], P[2, Con["a"]]}]};
zpowExplicitExpr[k_Integer] := Block[{$labelCounter = 0}, {Term[1, Join[ZVertex[{0, k, 0}], {A[3, Cov["a"]], P[2, Con["a"]]}]]}];
zpowCleanupTests[] := Module[{n2 = Global`n2, sym, wrong, expected, explicitRows},
  sym = ToYZPolynomial[CleanA[zpowLeftExpr[n2], Sequence @@ opts]];
  wrong = ToYZPolynomial[CleanA[zpowWrongSideExpr[n2], Sequence @@ opts]];
  expected = ZPowerSimplify[n2*y[1]*z[2]^(n2 - 1)];
  explicitRows = Table[
    With[{symK = ZPowerSimplify[sym /. n2 -> k], expK = ToYZPolynomial[CleanA[zpowExplicitExpr[k], Sequence @@ opts]]},
      <|"k" -> k, "SymbolicAtK" -> symK, "ExplicitExpansion" -> expK, "Difference" -> ZPowerSimplify[symK - expK]|>
    ],
    {k, 1, 4}
  ];
  <|
    "LeftCrossingInput" -> zpowLeftExpr[n2],
    "LeftCrossingPolynomial" -> sym,
    "LeftCrossingExpected" -> expected,
    "LeftCrossingDifference" -> ZPowerSimplify[sym - expected],
    "WrongSideInput" -> zpowWrongSideExpr[n2],
    "WrongSidePolynomial" -> wrong,
    "WrongSideExpected" -> 0,
    "WrongSideDifference" -> ZPowerSimplify[wrong],
    "ExplicitExpansionComparisons" -> explicitRows
  |>
];
zpowTests = zpowCleanupTests[];
write["zpow_a_cleanup_orientation_tests.wl", zpowTests];
write["zpow_a_cleanup_left_crossing_difference.txt", zpowTests["LeftCrossingDifference"]];
write["zpow_a_cleanup_wrong_side_difference.txt", zpowTests["WrongSideDifference"]];

(* Built-in single-stage and block-stage tests. *)
stageSuite = ReducerStageTestSuite[Sequence @@ opts];
write["reducer_stage_test_suite.wl", stageSuite];

stageWalks = <|
  "A_cleanup_only" -> StagePipelineRun[{"ACleanup"}, "U1^aA1_a", Sequence @@ opts],
  "Trace_removal_only" -> StagePipelineRun[{"TraceRemoval"}, "U3^aU3_a", Sequence @@ opts],
  "Yp_block_then_cleanup" -> StagePipelineRun[{"StartCleanup", "YpBlock", "ACleanup"}, "U3_aP2^aU2_bP3^b", Sequence @@ opts],
  "Div_block_then_cleanup" -> StagePipelineRun[{"StartCleanup", "DivBlock", "ACleanup"}, "U1_aP1^aU2_bP3^b", Sequence @@ opts],
  "PP_block_then_cleanup" -> StagePipelineRun[{"StartCleanup", "PPBlock", "ACleanup"}, "P1_aP2^aU1_bP2^b", Sequence @@ opts],
  "Box_block_then_cleanup" -> StagePipelineRun[{"StartCleanup", "BoxBlock", "ACleanup"}, "P1_aP1^aU2_bU3^b", Sequence @@ opts],
  "One_full_cycle_on_trace_free_sibling" -> StagePipelineRun[{"FullCycle"}, terms["trace_free_sibling"], Sequence @@ opts],
  "Grouped_pipeline_on_trace_free_sibling" -> StagePipelineRun[{"StartCleanup", "YpBlock", "ACleanup", "DivBlock", "ACleanup", "PPBlock", "ACleanup", "BoxBlock", "EndCleanup"}, terms["trace_free_sibling"], Sequence @@ opts],
  "Full_reduce_on_trace_free_sibling" -> StagePipelineRun[{"FullReduce"}, terms["trace_free_sibling"], Sequence @@ opts]
|>;
write["stage_walkthroughs.wl", stageWalks];

lich = LichnerowiczCommutatorTests[Sequence @@ opts];
lichDetails = LichnerowiczCommutatorDetails[Sequence @@ opts];
write["lichnerowicz_commutator_tests.wl", lich];
write["lichnerowicz_commutator_details.wl", lichDetails];

write["README.txt", "This folder contains user-indexed U/P/A term reductions, stage-by-stage reducer examples, ZPow/A-cleanup orientation tests, and local Lichnerowicz commutator checks.\n"];
Print["ZPow left-crossing difference: ", zpowTests["LeftCrossingDifference"]];
Print["ZPow wrong-side difference: ", zpowTests["WrongSideDifference"]];
Print["Lichnerowicz matches: ", lich[[All, "Matches"]]];
Print["Outputs written to folder: ", outDir];
