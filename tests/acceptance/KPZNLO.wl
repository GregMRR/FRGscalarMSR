(* ::Package:: *)

(* Suite hook of benchmarks/KPZNLO (design/NLOPIPELINE.md): KPZ at NLO, f_lambda = 1 and flowing. The benchmark's logic lives in
   benchmarks/KPZNLO; this file wraps its checks in test[...], adds a discrimination test and the Python round trips (v4, w = 0 included). *)

Get[FileNameJoin[{testRoot, "reference", "KPZNLOSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "PythonRoundTrip.wl"}]];
Get[FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZNLO", #}]] & /@ {"Reference.wl", "Checks.wl"};

Scan[test["KPZ NLO benchmark: " <> #Name <> " [" <> #Paper <> "]", #Pass, True] &, kpznloChecks[]];

(* the reference comparison discriminates: the w^2 in the A_a kernel replaced by 2 w^2 is detected *)
test["KPZ NLO benchmark: a perturbed reference is rejected",
  ExactlyZeroQ[(kpzlamnloAdim["Loop"][[1]] /. kpznloToRef) - (kpzlamnloRefLoop["D"] /. w -> 2 w)], False];

pythonRoundTrip["KPZ NLO", kpznloAdim, "KPZ", alamk^2];
pythonRoundTrip["KPZ NLO + f_lambda", kpzlamnloAdim, "KPZLambda", alamk^2];
