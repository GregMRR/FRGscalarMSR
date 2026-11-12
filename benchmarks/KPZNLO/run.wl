(* ::Package:: *)

(* run.wl -- the KPZ NLO benchmark, headless.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZNLO/run.wl
   Exit code 0 iff every check passes (a failed library check inside the pipeline exits 2 by itself). About a minute. *)

kpznloDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[kpznloDir, 2], "FRGscalarMSR.wl"}]];
Get[FileNameJoin[{ParentDirectory[kpznloDir], #, "Pipeline.wl"}]] & /@ {"KPZLO", "KPZLOLambda"};
Get[FileNameJoin[{kpznloDir, #}]] & /@ {"Pipeline.wl", "Reference.wl", "Checks.wl"};

kpznloResults = kpznloChecks[];
Print[""];
Scan[Print[If[#Pass, "PASS  ", "FAIL  "], #Name, "  [", #Paper, "]  ", #Detail] &, kpznloResults];
Print["==== KPZ NLO benchmark: ", Count[kpznloResults, KeyValuePattern["Pass" -> True]], " PASS, ",
  Count[kpznloResults, KeyValuePattern["Pass" -> False]], " FAIL ===="];
Exit[If[AllTrue[kpznloResults, #Pass &], 0, 1]];
