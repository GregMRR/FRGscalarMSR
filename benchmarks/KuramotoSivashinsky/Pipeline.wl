(* ::Package:: *)

(* Kuramoto-Sivashinsky, simplest DE ansatz of Gosteva, Wschebor, Canet, arXiv:2607.15784, Eq. (16):
     Gamma = Int hbar {dt h - (nu - tau dx^2) dx^2 h - lam/2 (dx h)^2 - (D - zeta dx^2) hbar}
   so Gamma^(1,1) = i w + p^2 (nu + tau p^2) and Gamma^(0,2) = -2 (D + zeta p^2) (Eq. 17).
   Regulators of Eq. (31) with two independent shapes: R^nu = D rn(q^2/k^2) (added as p^2 R^nu), R^D = D rd(q^2/k^2).
   Dim stays symbolic: the paper is d = 1, which the checks substitute; d = 2 is the control of the time-reversal check.
   Built once per kernel (ksEquations). Needs FRGscalarMSR loaded. *)

If[! ValueQ[ksEquations],
  ksTheory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim];
  ksAnsatz = DefineAnsatz[ksTheory,
    "Density" -> hbar (dt[h] - nuk lap[h] + tauk lap[lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar^2 - zetak grad[hbar].grad[hbar],
    "Couplings" -> {nuk, tauk, lamk, Dk, zetak}];
  ksSymmetries = DeclareSymmetries[ksAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity (the 1 in front of dt h)"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    Neglected[_, "Reason" -> "beyond the simplest DE ansatz"],
    NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
    "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 2|>];
  ksVertices = BuildVertices[ksAnsatz];
  ksRegulator = DefineRegulator[ksAnsatz, {hbar, h} -> Dk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
  ksPropagator = BuildPropagator[ksVertices, ksRegulator];
  ksLoops = BuildLoopIntegrands[ksVertices, ksPropagator, ksRegulator];
  ksProjected = ProjectOnMonomials[ksLoops, ksSymmetries];
  ksIntegrated = IntegrateFrequencies[ksProjected];
  ksClosure = CheckClosure[ksIntegrated, ksSymmetries];
  ksEquations = AssembleEquations[ksClosure]];
