(* ::Package:: *)

(* Declared d = 2 regression (design/LEVICIVITA.md): shallow water with "Dimension" -> FixedDimension[2] gives the same equations
   as the symbolic-Dim pipeline at Dim -> 2. Compared with ExactlyZeroQ, not ===: Assemble picks its output form by LeafCount. *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

Module[{th, an, sy, ve, re, pr, lo, eq, old = swaEquations["Equations"], new, key},
  th = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]];
  an = quiet[DefineAnsatz[th, "Density" -> swaAnsatz["Density"], "Couplings" -> swaAnsatz["Couplings"]]];
  ve = quiet[BuildVertices[an]];
  re = quiet[DefineRegulator[an,
     {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2],
     {Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2],
     "Shapes" -> {arnh, arnphi, armh, armphi}]];
  pr = quiet[BuildPropagator[ve, re]];
  sy = quiet[DeclareSymmetries[an,
     Protected[{Hbar, _}, "Momentum" -> 0, "Reason" -> "Hbar -> Hbar + c(t) (mass conservation)"],
     Protected[{Psi, _}, "Momentum" -> 0, "Reason" -> "Psi -> Psi + c(t)"],
     Neglected[{H, Psibar}, "Frequency" -> {1, 2}, "Momentum" -> 0, "Reason" -> "truncation, z = 1 power counting"],
     Neglected[{Psibar, Psibar}, "Frequency" -> 2, "Momentum" -> 0, "Reason" -> "coloured forcing, dropped"],
     Neglected[_, "Frequency" -> {1, 2}, "Momentum" -> 1, "Reason" -> "mixed i w q^2 / w^2 q^2 terms, dropped"],
     NonRenormalised[lk, "Reason" -> "Galilean invariance"],
     "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
  lo = quiet[BuildLoopIntegrands[ve, pr, re]];
  eq = quiet[AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[lo, sy]], sy]]];
  new = eq["Equations"];
  key = {#Pair, #Frequency, #Momentum, #Flowing} &;
  test["declared d = 2: same equations (pair, orders, flowing combination)", key /@ new, key /@ old];
  test["declared d = 2: classification unchanged", sy["Classification"][[All, {"Pair", "Frequency", "Momentum", "Class"}]],
   swaSymmetries["Classification"][[All, {"Pair", "Frequency", "Momentum", "Class"}]]];
  test["declared d = 2: every equation = symbolic-Dim equation at Dim -> 2 (exact)",
   AllTrue[Transpose[{new, old}], ExactlyZeroQ[#[[1]]["Integrand"] - (#[[2]]["Integrand"] /. Dim -> 2)] &], True];
  test["declared d = 2: no Epsilon key", KeyExistsQ[an, "Epsilon"], False]];
