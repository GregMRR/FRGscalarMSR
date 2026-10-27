(* ::Package:: *)

(* KPZ at LO (functions of p^2), as in BMWflow's ModelKPZ(approximation = "LO") and notes/LO_flow_equations.tex:
     Gamma = Int hbar {dt h - nu fnu(-Laplacian) lap h - lam/2 (grad h)^2} - D hbar fD(-Laplacian) hbar
   so Gamma^(1,1) = i w + nu p^2 fnu(p^2) and Gamma^(0,2) = -2 D fD(p^2) (design/LOPIPELINE.md section 4).
   Regulators as in kpz_notation.tex, with two independent shapes: R_nu = nu q^2 rn(q^2/k^2), R_D = D rd(q^2/k^2).
   Dim stays symbolic. Built once per kernel (kpzloAdim). Needs FRGscalarMSR loaded. *)

If[! ValueQ[kpzloAdim],
  kpzloTheory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim];
  kpzloAnsatz = DefineAnsatz[kpzloTheory,
    "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
    "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}];
  kpzloSymmetries = DeclareSymmetries[kpzloAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity (the 1 in front of dt h)"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    NonRenormalised[lamk, "Reason" -> "Galilean invariance"]];
  kpzloVertices = BuildVertices[kpzloAnsatz];
  kpzloRegulator = DefineRegulator[kpzloAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
  kpzloPropagator = BuildPropagator[kpzloVertices, kpzloRegulator];
  kpzloLoops = BuildLoopIntegrands[kpzloVertices, kpzloPropagator, kpzloRegulator];
  kpzloEquations = BuildLOEquations[kpzloLoops, kpzloSymmetries];
  kpzloDims = DeriveDimensions[kpzloAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}];
  kpzloAdim = AdimensionalEquations[kpzloEquations, kpzloDims]];
