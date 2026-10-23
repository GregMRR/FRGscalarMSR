(* ::Package:: *)

(* KPZ with momentum functions nu_k Fnu(p^2), D_k FD(p^2) (design/MOMENTUMFUNCTIONS.md; the LO scheme: f^lambda = 1), through the
   Equations stage, with the regulator of KPZSetup.wl (shapes rn, rd), shared by unit and acceptance tests (built once per run) *)

If[! ValueQ[kpzloEquations],
  kpzloTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
  kpzloAnsatz = quiet[DefineAnsatz[kpzloTheory, "Density" -> hbar (dt[h] - nuk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2,
     "Couplings" -> {nuk, lamk, Dk}, "MomentumFunctions" -> {Fnu, FD}]];
  kpzloSymmetries = quiet[DeclareSymmetries[kpzloAnsatz,
     Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
     Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
     Neglected[_, "Reason" -> "beyond the LO ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
     "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
  kpzloVertices = quiet[BuildVertices[kpzloAnsatz]];
  kpzloRegulator = quiet[DefineRegulator[kpzloAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]];
  kpzloPropagator = quiet[BuildPropagator[kpzloVertices, kpzloRegulator]];
  kpzloLoops = quiet[BuildLoopIntegrands[kpzloVertices, kpzloPropagator, kpzloRegulator]];
  kpzloProjection = quiet[ProjectOnMonomials[kpzloLoops, kpzloSymmetries]];
  kpzloIntegrated = quiet[IntegrateFrequencies[kpzloProjection]];
  kpzloClosure = quiet[CheckClosure[kpzloIntegrated, kpzloSymmetries]];
  kpzloEquations = quiet[AssembleEquations[kpzloClosure]]];
