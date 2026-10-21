(* ::Package:: *)

(* 1-scalar KPZ pipeline (dimensionful), shared by unit and acceptance tests (built once per run) *)

(* the Litim shape (1/y - 1) theta(1 - y) is only C^0 at y = 1: its r''' and r'''' (delta, delta') are lost by the termwise y < 1
   integrals of the KPZ tests. p^2 flows stop at r''; p^(2n >= 4) flows need r^(2n) (design/MOMENTUMPROJECTION.md step 4) *)
litimSmoothGuard[e_] := FRGAssert["Litim.smoothShape", FreeQ[e, Derivative[kk_ /; kk >= 3][rn | rd]], Union[Cases[e, Derivative[kk_ /; kk >= 3][rn | rd], {0, Infinity}]],
   "an explicit non-smooth (Litim) shape cannot be used with p^4 or higher flows"];

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
