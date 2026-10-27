(* ::Package:: *)

(* Propagator: inversion in atoms, Delta, validity conditions, degree assertions *)

pTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
pAnsatz = quiet[DefineAnsatz[pTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
pReg = quiet[DefineRegulator[pAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]];
kp = quiet[BuildPropagator[quiet[BuildVertices[pAnsatz]], pReg]];
fA = GammaAtom[h, hbar, 0][s0]; fDhh = GammaAtom[hbar, hbar, 0][s0];

test["KPZ atoms", Sort[kp["Atoms"]], Sort[{GammaAtom[h, hbar, 0], GammaAtom[hbar, hbar, 0]}]];
test["KPZ atom value", Expand[GammaAtom[h, hbar, 0][s0] /. AtomRules[kp]], Expand[nuk s0 + nuk s0 rn[s0/k^2]]];
test["KPZ Delta(w) = I w + f", kp["Delta"][ω, s0], I ω + fA];
test["KPZ G (response function and correlation)", PropagatorMatrix[kp, ω, s0],
  {{-fDhh/(\[CapitalDelta][-ω, s0] \[CapitalDelta][ω, s0]), 1/\[CapitalDelta][-ω, s0]}, {1/\[CapitalDelta][ω, s0], 0}}];
test["KPZ validity c0/c1 > 0", kp["ValidityConditions"][s0], {DeltaCoefficient[0][s0]/DeltaCoefficient[1][s0] > 0}];
test["KPZ Delta coefficients", kp["DeltaCoefficients"][s0], {fA, 1, 0}];

(* second-order time derivative, N = 1: Delta of degree 2 is allowed *)
p2 = quiet[BuildPropagator[quiet[BuildVertices[#]], quiet[DefineRegulator[#, {hbar, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}]]]] &[
   quiet[DefineAnsatz[pTheory, "Density" -> hbar (dt[h] + tau dt[dt[h]] - nuk lap[h]) - Dk hbar^2, "Couplings" -> {nuk, tau, Dk}]]];
test["dt dt: c2 = atom of (I w)^2, value tau", {p2["DeltaCoefficients"][s0][[3]], GammaAtom[h, hbar, 2][s0] /. AtomRules[p2]}, {GammaAtom[h, hbar, 2][s0], tau}];
test["dt dt: atom for the (I w)^2 coefficient", MemberQ[p2["Atoms"], GammaAtom[h, hbar, 2]], True];

(* negative: degree 3, degree 0, mismatched stages *)
twoTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
deg3 = quiet[DefineAnsatz[twoTheory, "Density" -> Hbar (dt[H] + tau dt[dt[H]] - nk lap[H]) + Psibar (dt[Psi] - nk lap[Psi]) - Dk Hbar^2,
    "Couplings" -> {nk, tau, Dk}]];
expectFail["degree-3 Delta aborts", BuildPropagator[quiet[BuildVertices[deg3]], quiet[DefineRegulator[deg3, {Hbar, H} -> nk q^2 rn[q^2/k^2], "Shapes" -> {rn}]]],
  "Propagator.deltaDegree"];
deg0 = quiet[DefineAnsatz[pTheory, "Density" -> hbar (-nuk lap[h]) - Dk hbar^2, "Couplings" -> {nuk, Dk}]];
expectFail["Delta without w aborts", BuildPropagator[quiet[BuildVertices[deg0]], quiet[DefineRegulator[deg0, {hbar, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}]]],
  "Propagator.deltaDegree"];
expectFail["regulator from another ansatz", BuildPropagator[quiet[BuildVertices[deg0]], pReg], "Propagator.sameAnsatz"];
expectFail["arguments swapped", BuildPropagator[pReg, quiet[BuildVertices[pAnsatz]]], "Stage.input"];
