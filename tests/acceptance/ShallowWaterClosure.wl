(* ::Package:: *)

(* Milestone 5 acceptance: closure of the shallow-water DE4plus ansatz *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

(* normalisations = the script's jobs LHS factors / (2 Dim)^n *)
swNorms = {{nk Ak, 1}, {gk, 1}, {g2k, 1}, {Dhk aDhk, -2}, {Hk lk, -1}, {Kk, 1}, {nk ank, 1}, {Dphik aDphik, -2}, {Dphik D2phik, -2}};
test["9 equations with the script's flowing combinations and normalisations",
  Sort[{#Flowing, #Normalisation} & /@ swaClosure["Equations"]], Sort[swNorms]];
neglectedValue[pair_, m_, n_] := First[SelectFirst[swaClosure["Report"], #Pair === pair && #Frequency == m && #Momentum == n &]["Samples"]];
test["script section 6b: i w and w^2 in Gamma_{H Psibar}, w^2 in Gamma_{Psibar Psibar} reported as generated (nonzero at a valid point)",
  FreeQ[{neglectedValue[{H, Psibar}, 1, 0], neglectedValue[{H, Psibar}, 2, 0], neglectedValue[{Psibar, Psibar}, 2, 0]}, 0, {1}], True];
test["section 6b terms are WARN lines of the card", Count[swaClosure["Warnings"], KeyValuePattern[{"Status" -> "WARN"}]] >= 3, True];
test["Protected, DeclaredNumeric and Causality checks passed",
  SubsetQ[swaClosure["Checks"][[All, "Check"]], {"Closure.protected", "Closure.declaredNumeric", "Closure.causality"}], True];
test["4 monomials vanish by structure", Count[swaClosure["Report"], KeyValuePattern["Class" -> "Structure"]], 4];
