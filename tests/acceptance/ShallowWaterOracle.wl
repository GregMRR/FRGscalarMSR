(* ::Package:: *)

(* Milestone 5 acceptance: the numerical oracle agrees with the 9 assembled equations in d = 2 and d = 3 *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

swaOracle = quiet[RunOracle[swaEquations]];
test["oracle: 18 comparisons (9 equations, d = 2 and 3)", Length[swaOracle["Results"]], 18];
test["oracle: relative error below 1e-8 (finite differences, h = 1e-5)", Max[swaOracle["Results"][[All, "Error"]]] < 10^-8, True];
