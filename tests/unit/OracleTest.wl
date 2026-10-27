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

(* (Laplacian_p)^n / N_n by nested stencils (design/MOMENTUMPROJECTION.md step 5) on F((p + q)^2), F = 1/(1 + x): exact coefficient
   F''/2 + 2 Q F'''/d + 2 Q^2 F''''/(d (d + 2)) at n = 2; values at 60 digits, h = 10^-5, Richardson (h, h/2) *)
Module[{F = Function[x, 1/(1 + x)], q = {3/10, -1/2, 2/5}, exact, num, f},
  Do[With[{qd = Take[q, d]}, Module[{Q = qd . qd},
     f = Function[pv, N[F[(pv + qd) . (pv + qd)], 60]];
     exact = F''[Q]/2 + 2 Q F'''[Q]/d + 2 Q^2 F''''[Q]/(d (d + 2));
     num = FRGscalarMSR`Private`oracleLaplacianN[f, d, 2, 10^-5]/FRGscalarMSR`Private`oracleLapNorm[2, d];
     test["oracle stencil n = 2, d = " <> ToString[d] <> ": rel. error < 10^-9", Abs[(num - exact)/exact] < 10^-9, True];
     test["oracle stencil n = 1, d = " <> ToString[d] <> ": N_1 = 2 d", Abs[FRGscalarMSR`Private`oracleLaplacianN[f, d, 1, 10^-5]/(2 d) - (F'[Q] + 2 Q F''[Q]/d)] < 10^-9, True]]],
   {d, 2, 3}]];

(* LO branch (design/LOPIPELINE.md step 6): residue kernels at random p, q, c = numerical w integral of the loop, d = 2 and 3 *)
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
oLO = quiet[RunOracle[kpzloEquations, "ZeroFlows" -> True]];
test["LO oracle: both functions and the zero entry, d = 2 and 3", {oLO["Pipeline"], Length[oLO["Results"]], Union[oLO["Results"][[All, "Dimension"]]]}, {"LO", 6, {2, 3}}];
test["LO oracle: every result within 10^-8", AllTrue[oLO["Results"], #Error < 10^-8 &], True];
