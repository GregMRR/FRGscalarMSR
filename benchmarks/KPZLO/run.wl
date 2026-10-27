(* ::Package:: *)

(* run.wl -- the KPZ LO benchmark, headless.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZLO/run.wl
   Exit code 0 iff every check passes (a failed library check inside the pipeline exits 2 by itself). A few seconds. *)

kpzloDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[kpzloDir, 2], "FRGscalarMSR.wl"}]];
Get[FileNameJoin[{kpzloDir, #}]] & /@ {"Pipeline.wl", "Reference.wl", "Checks.wl"};

kpzloResults = kpzloChecks[];
Print[""];
Scan[Print[If[#Pass, "PASS  ", "FAIL  "], #Name, "  [", #Paper, "]  ", #Detail] &, kpzloResults];
Print["==== KPZ LO benchmark: ", Count[kpzloResults, KeyValuePattern["Pass" -> True]], " PASS, ", Count[kpzloResults, KeyValuePattern["Pass" -> False]], " FAIL ===="];
Exit[If[AllTrue[kpzloResults, #Pass &], 0, 1]];
