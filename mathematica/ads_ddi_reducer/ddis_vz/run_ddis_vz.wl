(* ::Package:: *)
(* Produce strict reducer outputs for DDI1--DDI10 at symbolic V(z)=z1^n1 z2^n2 z3^n3. No reference formulas are loaded here. *)
If[$InputFileName =!= "", SetDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
Get[FileNameJoin[{Directory[], "package", "AdSDDIReducer.wl"}]];

opts = {"MaxCycles" -> 40, "MaxSubsteps" -> 10000, "Lambda" -> \[Lambda]AdS, "Dim" -> 3, "Mass" -> (m[#] &)};
fields = If[TrueQ[ToExpression[Environment["RUN_CYCLIC_SYMBOLIC"] /. $Failed -> "False"]], {1,2,3}, {1}];
keys = Switch[Environment["KINDS"] /. $Failed -> "ALL", "B1B2", {"DDI1","DDI2"}, "B3B5", {"DDI3","DDI4","DDI5"}, "B6B10", {"DDI6","DDI7","DDI8","DDI9","DDI10"}, _, Table["DDI"<>ToString[r], {r,1,10}]];
n = SymbolicZPowers[];
outDir = FileNameJoin[{Directory[], "ddi_vz_outputs"}];
If[DirectoryQ[outDir], DeleteDirectory[outDir, DeleteContents -> True]]; CreateDirectory[outDir];
write[name_, expr_] := Export[FileNameJoin[{outDir, name}], ToString[expr, InputForm], "Text"];
runOne[k_, i_] := Module[{raw, red, poly, report, base},
  base = k <> "_i" <> ToString[i] <> "_symbolic_n1_n2_n3";
  Print["Running ", base];
  raw = DDIRawByKey[k, i, n];
  red = ReduceAdS[raw, Sequence @@ opts];
  poly = ToYZPolynomial[red];
  report = NonCanonicalReport[red];
  write[base <> "_polynomial.txt", poly];
  write[base <> "_report.txt", report];
  write[base <> "_canonical.txt", CanonicalQ[red]];
  write[base <> "_raw_term_count.txt", TermCount[raw]];
  write[base <> "_reduced_term_count.txt", TermCount[red]];
  <|"Kind"->k,"i"->i,"Powers"->n,"RawTermCount"->TermCount[raw],"ReducedTermCount"->TermCount[red],"CanonicalQ"->CanonicalQ[red],"Polynomial"->poly,"Report"->report|>
];
results = Association@Flatten[Table[(k<>"_i"<>ToString[i])->runOne[k,i], {i,fields},{k,keys}],1];
write["summary.wl", results];
write["README.txt", "Strict reducer outputs for symbolic V(z)=z1^n1 z2^n2 z3^n3. No reference formulas are used by this script.\n"];
Print["Outputs written to folder: ", outDir];
