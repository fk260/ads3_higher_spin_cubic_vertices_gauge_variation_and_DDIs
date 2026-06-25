(* ::Package:: *)
(* Produce strict reducer outputs for DDI1--DDI10 at V=1.  No reference formulas are loaded here. *)
If[$InputFileName =!= "", SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
Get[FileNameJoin[{Directory[], "package", "AdSDDIReducer.wl"}]];

opts = {"MaxCycles" -> 36, "MaxSubsteps" -> 8000, "Lambda" -> \[Lambda]AdS, "Dim" -> 3, "Mass" -> (m[#] &)};
fields = If[TrueQ[ToExpression[Environment["RUN_CYCLIC_DDIS"] /. $Failed -> "False"]], {1,2,3}, {1}];
keys = Table["DDI" <> ToString[r], {r, 1, 10}];
outDir = FileNameJoin[{Directory[], "ddi_v1_outputs"}];
If[DirectoryQ[outDir], DeleteDirectory[outDir, DeleteContents -> True]]; CreateDirectory[outDir];
write[name_, expr_] := Export[FileNameJoin[{outDir, name}], ToString[expr, InputForm], "Text"];
runOne[k_, i_] := Module[{raw, red, poly, report, base},
  base = k <> "_i" <> ToString[i] <> "_V1";
  Print["Running ", base];
  raw = DDIRawByKey[k, i, {0,0,0}];
  red = ReduceAdS[raw, Sequence @@ opts];
  poly = ToYZPolynomial[red];
  report = NonCanonicalReport[red];
  write[base <> "_polynomial.txt", poly];
  write[base <> "_report.txt", report];
  write[base <> "_canonical.txt", CanonicalQ[red]];
  write[base <> "_raw_term_count.txt", TermCount[raw]];
  write[base <> "_reduced_term_count.txt", TermCount[red]];
  <|"Kind"->k,"i"->i,"RawTermCount"->TermCount[raw],"ReducedTermCount"->TermCount[red],"CanonicalQ"->CanonicalQ[red],"Polynomial"->poly,"Report"->report|>
];
results = Association@Flatten[Table[(k<>"_i"<>ToString[i])->runOne[k,i], {i,fields},{k,keys}],1];
write["summary.wl", results];
write["README.txt", "Strict reducer outputs for DDI1--DDI10 at V=1. No reference formulas are used by this script.\n"];
Print["Outputs written to folder: ", outDir];
