(* ::Package:: *)

(* run.wl -- the Kuramoto-Sivashinsky benchmark, headless.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KuramotoSivashinsky/run.wl
   Exit code 0 iff every check passes (a failed library check inside the pipeline exits 2 by itself). About 30 s. *)

ksDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[ksDir, 2], "FRGscalarMSR.wl"}]];
Get[FileNameJoin[{ksDir, #}]] & /@ {"Pipeline.wl", "Reference.wl", "Checks.wl"};

ksResults = ksChecks[];
Print[""];
Scan[Print[If[#Pass, "PASS  ", "FAIL  "], #Name, "  [", #Paper, "]  ", #Detail] &, ksResults];
Print["==== KS benchmark: ", Count[ksResults, KeyValuePattern["Pass" -> True]], " PASS, ", Count[ksResults, KeyValuePattern["Pass" -> False]], " FAIL ===="];
Exit[If[AllTrue[ksResults, #Pass &], 0, 1]];
