(* ::Package:: *)

(* 1-scalar KPZ pipeline (dimensionful), shared by unit and acceptance tests (built once per run) *)

If[! ValueQ[kpzEquations],
  kpzTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
  kpzAnsatz = quiet[DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
  kpzSymmetries = quiet[DeclareSymmetries[kpzAnsatz,
     Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
     Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
     Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
     "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
  kpzVertices = quiet[BuildVertices[kpzAnsatz]];
  kpzRegulator = quiet[DefineRegulator[kpzAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]];
  kpzPropagator = quiet[BuildPropagator[kpzVertices, kpzRegulator]];
  kpzLoops = quiet[BuildLoopIntegrands[kpzVertices, kpzPropagator, kpzRegulator]];
  kpzIntegrated = quiet[IntegrateFrequencies[ProjectOnMonomials[kpzLoops, kpzSymmetries]]];
  kpzClosure = quiet[CheckClosure[kpzIntegrated, kpzSymmetries]];
  kpzEquations = quiet[AssembleEquations[kpzClosure]]];
