(* ::Package:: *)

(* KPZ at LO with a flowing f_lambda (design/LOFLAMBDA.md), the Kloss-Canet-Wschebor LO ansatz:
     Gamma = Int hbar f_lambda(-Laplacian) {dt h - lam/2 (grad h)^2} - nu hbar fnu(-Laplacian) lap h - D hbar fD(-Laplacian) hbar
   so Gamma^(1,1) = i w f_lambda(p^2) + nu p^2 fnu(p^2), Gamma^(0,2) = -2 D fD(p^2), and the cubic vertex carries f_lambda at the
   momentum of its hbar leg (Galilean invariance ties it to the f_lambda of dt h). f_lambda has reference 1: f_lambda(0) = 1 is the
   Galilean Ward identity, checked at p = 0; it has no eta. Its flow is d/d(i nu) of d_s Gamma^(1,1)(nu, p) at nu = 0, taken before
   the frequency integral. Regulators as in benchmarks/KPZLO. Dim stays symbolic. Built once per kernel (kpzlamAdim).
   Needs FRGscalarMSR loaded. *)

If[! ValueQ[kpzlamAdim],
  kpzlamTheory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim];
  kpzlamAnsatz = DefineAnsatz[kpzlamTheory,
    "Density" -> fn[flam, hbar] (dt[h] - lamk/2 grad[h].grad[h]) - nuk hbar fn[fnu, lap[h]] - Dk hbar fn[fD, hbar],
    "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk, flam -> 1}];
  kpzlamSymmetries = DeclareSymmetries[kpzlamAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity: f_lambda(0) = 1"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    NonRenormalised[lamk, "Reason" -> "Galilean invariance"]];
  kpzlamVertices = BuildVertices[kpzlamAnsatz];
  kpzlamRegulator = DefineRegulator[kpzlamAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
  kpzlamPropagator = BuildPropagator[kpzlamVertices, kpzlamRegulator];
  kpzlamLoops = BuildLoopIntegrands[kpzlamVertices, kpzlamPropagator, kpzlamRegulator];
  kpzlamEquations = BuildLOEquations[kpzlamLoops, kpzlamSymmetries];
  kpzlamDims = DeriveDimensions[kpzlamAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}];
  kpzlamAdim = AdimensionalEquations[kpzlamEquations, kpzlamDims]];

(* The same model in the derivative expansion with f_lambda(s) = 1 + xk s, fnu = fD = 1: (hbar - xk lap hbar)(dt h - lam/2 (grad h)^2).
   Its flows of xk (the i w p^2 coefficient), nuk and Dk are the p^2 / p^0 Taylor coefficients of the LO rows at those functions:
   an independent route (DE projection, single-family residues, Assemble's normalisation) for the sign and normalisation of the
   d/d(i nu) row (Checks.wl, KPZLOLambda.de.*; design/reviews/LOFLambda-r1.md E1). *)
If[! ValueQ[kpzlamDEEquations],
  kpzlamDEAnsatz = DefineAnsatz[kpzlamTheory,
    "Density" -> (hbar - xk lap[hbar]) (dt[h] - lamk/2 grad[h].grad[h]) - nuk hbar lap[h] - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk, xk}];
  kpzlamDESymmetries = DeclareSymmetries[kpzlamDEAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
    "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>];
  kpzlamDEVertices = BuildVertices[kpzlamDEAnsatz];
  kpzlamDERegulator = DefineRegulator[kpzlamDEAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
  kpzlamDEPropagator = BuildPropagator[kpzlamDEVertices, kpzlamDERegulator];
  kpzlamDELoops = BuildLoopIntegrands[kpzlamDEVertices, kpzlamDEPropagator, kpzlamDERegulator];
  kpzlamDEEquations = AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[kpzlamDELoops, kpzlamDESymmetries]], kpzlamDESymmetries]]];
