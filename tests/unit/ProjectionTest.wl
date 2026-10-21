(* ::Package:: *)

(* Projection: the p^(2n) directional average (design/MOMENTUMPROJECTION.md M1-M3), its normalisation, its runtime checks *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

prj[e_, n_, dim_ : Dim] := FRGscalarMSR`Private`project[e, 0, n, {tnu, tp}, {tw, tq}, dim, 0];
fpq = ff[scal[tp + tq, tp + tq]];

(* closed forms of the (p^2)^n coefficient of <f((p + q)^2)> (q^2 = Q) *)
test["p^2: f'(Q) + 2 Q f''(Q)/d", Together[prj[fpq, 1] - (ff'[QQ] + 2 QQ ff''[QQ]/Dim) /. QQ -> LoopMomentumSquared], 0];
test["p^4: f''/2 + 2 Q f'''/d + 2 Q^2 f''''/(d(d+2))",
  Together[prj[fpq, 2] - (ff''[QQ]/2 + 2 QQ ff'''[QQ]/Dim + 2 QQ^2 ff''''[QQ]/(Dim (Dim + 2))) /. QQ -> LoopMomentumSquared], 0];
test["p^6: f'''/6 + Q f''''/d + 2 Q^2 f^(5)/(d(d+2)) + 4 Q^3 f^(6)/(3 d(d+2)(d+4))",
  Together[prj[fpq, 3] - (ff'''[QQ]/6 + QQ ff''''[QQ]/Dim + 2 QQ^2 Derivative[5][ff][QQ]/(Dim (Dim + 2)) +
        4 QQ^3 Derivative[6][ff][QQ]/(3 Dim (Dim + 2) (Dim + 4))) /. QQ -> LoopMomentumSquared], 0];
test["p^4 closed form = explicit (Laplacian)^2/N_2 in d = 1..4 (sound ExactlyZeroQ on derivative atoms)",
  AllTrue[Range[4], Function[d, With[{x = FRGscalarMSR`Private`explicitProjection[fpq, 0, 2, {tnu, tp}, {tw, tq}, d, 0]},
     ExactlyZeroQ[x[[1]] - (prj[fpq, 2] /. LoopMomentumSquared -> x[[2]] /. Dim -> d)]]]], True];

(* M2: N_n = (Laplacian)^n (p^2)^n = 4^n n! (d/2)_n on explicit vectors *)
test["N_n = explicit (Laplacian)^n (p^2)^n, n = 1..3, d = 1..4",
  Table[First[FRGscalarMSR`Private`explicitProjection[scal[tp, tp]^n, 0, n, {tnu, tp}, {tw, tq}, d, 0]], {n, 3}, {d, 4}], ConstantArray[1, {3, 4}]];
test["N_1 = 2 d, N_2 = 8 d (d + 2)", {FRGscalarMSR`Private`lapNorm[1, Dim], Expand[FRGscalarMSR`Private`lapNorm[2, Dim]]}, {2 Dim, Expand[8 Dim (Dim + 2)]}];

(* M3: the part odd in wedge[p, q] averages to zero in d = 2, for every n *)
test["wedge-odd term: explicit p^4 projection in d = 2 vanishes",
  ExactlyZeroQ[First[FRGscalarMSR`Private`explicitProjection[wedge[tp, tq] scal[tp, tq] fpq, 0, 2, {tnu, tp}, {tw, tq}, 2, 0]]], True];
test["wedge-odd term dropped at p^4 (checked)", FRGscalarMSR`Private`project[wedge[tp, tq] scal[tp, tq] fpq + fpq, 0, 2, {tnu, tp}, {tw, tq}, 2, 0, True],
  prj[fpq, 2, 2]];
expectFail["wedgeOdd at n = 2: a wedge with a third vector", FRGscalarMSR`Private`project[wedge[tp, tr] scal[tp, tp]^2, 0, 2, {tnu, tp}, {tw, tq}, 2, 0, True],
  "Projection.wedgeOdd"];

(* the runtime checks pass on the correct table *)
test["normalisation check passes, n <= 3, symbolic d", FRGscalarMSR`Private`checkNormalisation[3, Dim], "Projection.normalisation"];
test["normalisation check passes, n <= 3, d = 2", FRGscalarMSR`Private`checkNormalisation[3, 2], "Projection.normalisation"];

(* a p^4 KPZ pipeline: the ansatz is p^2, every p^4 monomial is Neglected (no end-to-end p^4 acceptance test, by design) *)
kpz4Symmetries = quiet[DeclareSymmetries[kpzAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
    "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 2|>]];
kpz4 = quiet[ProjectOnMonomials[kpzLoops, kpz4Symmetries]];
test["KPZ p^4: projection builds", Count[kpz4["Integrands"], r_ /; r["Momentum"] == 2] > 0, True];
test["KPZ p^4: explicitVectors and normalisation on the card",
  SubsetQ[kpz4["Checks"][[All, "Check"]], {"Projection.explicitVectors", "Projection.normalisation", "Projection.LHS"}], True];
test["KPZ p^4: shape smoothness ASSUMPTION", MemberQ[kpz4["Warnings"][[All, "Check"]], "Projection.shapeSmoothness"], True];
test["KPZ p^4: the n <= 1 integrands are those of the p^2 pipeline",
  Module[{key = {#Pair, #Frequency, #Momentum} &, old = kpzIntegrated["Projection"]["Integrands"]},
   AllTrue[old, Function[r, Expand[SelectFirst[kpz4["Integrands"], key[#] === key[r] &]["Integrand"] - r["Integrand"]] === 0]]], True];

(* mutation (reviews E1, E4, E17): a wrong <c^4> = 3 Q^2/(7 d - 6) (= 3 Q^2/(d (d + 2)) at d = 2, 3) injected into projCAverage only *)
mutAvg[k_, dim_, qq_] := Which[OddQ[k], 0, k == 4, 3 qq^2/(7 dim - 6), True, qq^(k/2) (k - 1)!!/Product[dim + 2 j, {j, 0, k/2 - 1}]];
kpzLoop = LoopIntegrand[kpzLoops, {hbar, h}, tnu, tp];
kpzDeg = kpzLoops["Propagator"]["DeltaDegree"]; kpzWq = LoopVariables[kpzLoops];
kpzP4 = FRGscalarMSR`Private`project[kpzLoop, 0, 2, {tnu, tp}, kpzWq, Dim, kpzDeg];
test["explicitVectors passes on the correct KPZ p^4 integrand",
  FRGscalarMSR`Private`checkExplicitVectors[kpzLoop, kpzP4, 0, 2, {tnu, tp}, kpzWq, Dim, kpzDeg], "Projection.explicitVectors"];
Block[{FRGscalarMSR`Private`projCAverage = mutAvg},
  expectFail["mutated <c^4>: normalisation fails", FRGscalarMSR`Private`checkNormalisation[2, Dim], "Projection.normalisation"];
  expectFail["mutated <c^4>: explicitVectors fails",
   FRGscalarMSR`Private`checkExplicitVectors[kpzLoop, FRGscalarMSR`Private`project[kpzLoop, 0, 2, {tnu, tp}, kpzWq, Dim, kpzDeg], 0, 2, {tnu, tp}, kpzWq, Dim, kpzDeg],
   "Projection.explicitVectors"];
  expectFail["mutated <c^4>: ProjectOnMonomials aborts", ProjectOnMonomials[kpzLoops, kpz4Symmetries], "Projection.normalisation"]];
expectFail["explicitVectors: a wrong integrand", FRGscalarMSR`Private`checkExplicitVectors[kpzLoop, 2 kpzP4, 0, 2, {tnu, tp}, kpzWq, Dim, kpzDeg],
  "Projection.explicitVectors"];
expectFail["explicitVectors at a declared d = 2: a wrong integrand",
  FRGscalarMSR`Private`checkExplicitVectors[kpzLoop, (kpzP4 /. Dim -> 2) + LoopMomentumSquared^3, 0, 2, {tnu, tp}, kpzWq, 2, kpzDeg], "Projection.explicitVectors"];
