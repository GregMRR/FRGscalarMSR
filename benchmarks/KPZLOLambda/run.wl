(* ::Package:: *)

(* run.wl -- the KPZ LO + f_lambda benchmark, headless.
   Run:  math -script ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/run.wl
   Exit code 0 iff every check passes (a failed library check inside the pipeline exits 2 by itself). Under a minute. *)

kpzlamDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[kpzlamDir, 2], "FRGscalarMSR.wl"}]];
Get[FileNameJoin[{ParentDirectory[kpzlamDir], "KPZLO", "Pipeline.wl"}]];
Get[FileNameJoin[{kpzlamDir, #}]] & /@ {"Pipeline.wl", "Reference.wl", "Checks.wl"};

kpzlamResults = kpzlamChecks[];
Print[""];
Scan[Print[If[#Pass, "PASS  ", "FAIL  "], #Name, "  [", #Paper, "]  ", #Detail] &, kpzlamResults];
Print["==== KPZ LO + f_lambda benchmark: ", Count[kpzlamResults, KeyValuePattern["Pass" -> True]], " PASS, ",
  Count[kpzlamResults, KeyValuePattern["Pass" -> False]], " FAIL ===="];
Exit[If[AllTrue[kpzlamResults, #Pass &], 0, 1]];
