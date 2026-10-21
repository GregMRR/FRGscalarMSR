(* ::Package:: *)

(* Oracle: the dimensions it runs in (design/LEVICIVITA.md) *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

Module[{th, an, sy, ve, re, pr, lo, eq},
  th = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> FixedDimension[2]]];
  an = quiet[DefineAnsatz[th, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
  sy = quiet[DeclareSymmetries[an, Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
     "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>]];
  ve = quiet[BuildVertices[an]];
  re = quiet[DefineRegulator[an, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]];
  pr = quiet[BuildPropagator[ve, re]];
  lo = quiet[BuildLoopIntegrands[ve, pr, re]];
  eq = quiet[AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[lo, sy]], sy]]];
  test["declared d = 2: equations carry Dimension 2", eq["Dimension"], 2];
  test["declared d = 2: Automatic runs the oracle in d = 2 only", Union[quiet[RunOracle[eq]]["Results"][[All, "Dimension"]]], {2}];
  expectFail["declared d = 2: {2, 3} aborts", RunOracle[eq, "Dimensions" -> {2, 3}], "Oracle.declaredDimension"];
  expectFail["declared d = 2: {3} aborts", RunOracle[eq, "Dimensions" -> {3}], "Oracle.declaredDimension"]];

expectFail["symbolic d: {4} aborts", RunOracle[kpzEquations, "Dimensions" -> {4}], "Oracle.declaredDimension"];
test["symbolic d: Automatic = {2, 3}", Union[quiet[RunOracle[kpzEquations]]["Results"][[All, "Dimension"]]], {2, 3}];
