(* ::Package:: *)

(* FrequencyIntegration: residue integrals against the script's I(0,1,1) and against NIntegrate *)

cRule[c0_, c1_, c2_] := {FrequencyIntegralC[0] -> c0, FrequencyIntegralC[1] -> c1, FrequencyIntegralC[2] -> c2};
test["I(0,1,1) = Pi/(2 A C) as in the primer (c2 = 1, c1 = 2A, c0 = C)", Together[FrequencyIntegral[0, 1, 1, 2] /. cRule[CC, 2 A, 1]], Pi/(2 A CC)];
test["a = 0: no pole in the upper half plane", FrequencyIntegral[0, 0, 2, 2], 0];
test["results rational (no radicals)", FreeQ[Table[FrequencyIntegral[n, a, b, 2], {a, 1, 3}, {b, 0, 2}, {n, 0, 2 (a + b) - 2}], Power[_, _Rational]], True];
numI[n_, a_, b_, {c0_, c1_, c2_}] := NIntegrate[w^n/((c0 + c1 I w + c2 (I w)^2)^a (c0 - c1 I w + c2 (I w)^2)^b), {w, -Infinity, Infinity},
   WorkingPrecision -> 30, PrecisionGoal -> 15, AccuracyGoal -> 20] // Quiet;
Do[With[{c = {7/5, 3/2, 2/3}},   (* c1/c2 > 0, c0/c2 > 0: both roots in the upper half plane; overdamped *)
   test["deg 2 I" <> ToString[nab] <> " = NIntegrate", Abs[N[FrequencyIntegral[Sequence @@ nab, 2] /. cRule @@ c, 30] - numI[Sequence @@ nab, c]] < 10^-12, True]],
  {nab, {{0, 1, 1}, {2, 2, 1}, {1, 1, 2}, {3, 2, 2}, {0, 3, 0}}}];
With[{c = {3, 1/5, 2}},   (* underdamped: complex roots, 4 c0 c2 > c1^2 *)
  test["deg 2 underdamped I(2,2,1) = NIntegrate", Abs[N[FrequencyIntegral[2, 2, 1, 2] /. cRule @@ c, 30] - numI[2, 2, 1, c]] < 10^-12, True]];
With[{c = {5/3, 2/7, 0}},
  Do[test["deg 1 I" <> ToString[nab] <> " = NIntegrate", Abs[N[FrequencyIntegral[Sequence @@ nab, 1] /. cRule @@ c, 30] - numI[Sequence @@ nab, c]] < 10^-12, True],
   {nab, {{0, 1, 1}, {1, 2, 1}, {0, 2, 2}}}]];
expectFail["divergent integral aborts", FrequencyIntegral[4, 1, 1, 2], "FrequencyIntegration.convergence"];

(* two pole families, Delta at q and at Q = p + q, external frequency 0 (LO, design/LOPIPELINE.md step 4) *)
famC = FrequencyIntegralC;
test["two families: 1/(Delta_q(w) Delta_Q(-w)) = 2 Pi/(l_q + l_Q) (deg 1, c1 = 1)",
  Together[(FrequencyIntegralFamilies[0, {{1, 0, 1}, {0, 1, 1}}] /. {famC[1, 1] -> 1, famC[2, 1] -> 1}) - 2 Pi/(famC[1, 0] + famC[2, 0])], 0];
numI2[n_, fams_, vals_] := Module[{dd},
   dd[j_, x_] := Sum[vals[[j, m + 1]] (I x)^m, {m, 0, 2}];
   Quiet[NIntegrate[w^n/Product[dd[j, w]^fams[[j, 1]] dd[j, -w]^fams[[j, 2]], {j, 2}], {w, -Infinity, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 15, AccuracyGoal -> 20]]];
fiFam[n_, fams_, vals_] := N[FrequencyIntegralFamilies[n, fams] /. famC[j_, m_] :> vals[[j, m + 1]], 30];
Do[test["two families deg 1 " <> ToString[nf] <> " = NIntegrate", Abs[fiFam[Sequence @@ nf, {{3/2, 1, 0}, {7/5, 2/3, 0}}] - numI2[Sequence @@ nf, {{3/2, 1, 0}, {7/5, 2/3, 0}}]] < 10^-12, True],
  {nf, {{0, {{2, 2, 1}, {1, 1, 1}}}, {2, {{1, 0, 1}, {1, 2, 1}}}, {2, {{2, 2, 1}, {1, 1, 1}}}}}];
Do[test["two families deg 2 " <> ToString[nf] <> " = NIntegrate", Abs[fiFam[Sequence @@ nf, {{3/2, 1, 1/2}, {7/5, 2/3, 1/3}}] - numI2[Sequence @@ nf, {{3/2, 1, 1/2}, {7/5, 2/3, 1/3}}]] < 10^-12, True],
  {nf, {{0, {{2, 1, 2}, {1, 1, 2}}}, {2, {{1, 1, 2}, {1, 1, 2}}}}}];
(* at coincidence (Q data = q data) the two-family integral is the one-family integral of the merged powers *)
Do[test["two families at coincidence = FrequencyIntegral " <> ToString[nf], Together[(FrequencyIntegralFamilies @@ nf /. famC[_, m_] :> famC[m]) -
     FrequencyIntegral[nf[[1]], nf[[2, 1, 1]] + nf[[2, 2, 1]], nf[[2, 1, 2]] + nf[[2, 2, 2]], nf[[2, 1, 3]]]], 0],
  {nf, {{0, {{2, 2, 1}, {1, 1, 1}}}, {1, {{1, 0, 1}, {1, 2, 1}}}, {0, {{2, 1, 2}, {1, 1, 2}}}}}];
expectFail["a mixed-momentum Delta in the single-momentum dispatch aborts (it would give pi/l_q instead of 2 pi/(l_q + l_Q))",
  FRGscalarMSR`Private`integrate[1/(\[CapitalDelta][w, qq2] \[CapitalDelta][-w, QQ2]), w, 1, {}], "FrequencyIntegration.singleMomentum"];
expectFail["a two-family sum singular at coincidence aborts", FRGscalarMSR`Private`loCheckRegularAtCoincidence[1/(famC[1, 0] - famC[2, 0])],
  "FrequencyIntegration.regularAtCoincidence"];
expectFail["two families: divergent integral aborts", FrequencyIntegralFamilies[3, {{1, 0, 1}, {0, 1, 1}}], "FrequencyIntegration.convergence"];

(* shifted pole family (NLO, design/NLOPIPELINE.md N1): Delta_2(w + s), Delta_2(-w - s), real s *)
numShift[n_, fams_, vals_, s_] := Module[{dd},
   dd[j_, x_] := Sum[vals[[j, m + 1]] (I x)^m, {m, 0, 2}];
   Quiet[NIntegrate[w^n/(dd[1, w]^fams[[1, 1]] dd[1, -w]^fams[[1, 2]] dd[2, w + s]^fams[[2, 1]] dd[2, -w - s]^fams[[2, 2]]), {w, -Infinity, Infinity},
     WorkingPrecision -> 30, PrecisionGoal -> 15, AccuracyGoal -> 20]]];
fiShift[n_, fams_, vals_, s_] := N[FrequencyIntegralFamilies[n, {fams[[1]], Append[fams[[2]], "Shifted"]}] /. famC[j_, m_] :> vals[[j, m + 1]] /.
     FrequencyIntegralShift[2] -> s, 30];
Do[With[{n = nf[[1]], fams = nf[[2]], vals = nf[[3]]}, Do[
    test["shifted family " <> ToString[{n, fams}] <> " at s = " <> ToString[s, InputForm] <> " = NIntegrate",
     Abs[fiShift[n, fams, vals, s] - numShift[n, fams, vals, s]] < 10^-10 Abs[numShift[n, fams, vals, s]], True], {s, {3/5, -11/7}}]],
  {nf, {{0, {{2, 2, 1}, {1, 1, 1}}, {{3/2, 1, 0}, {7/5, 2/3, 0}}}, {2, {{1, 0, 1}, {1, 2, 1}}, {{3/2, 1, 0}, {7/5, 2/3, 0}}},
    {1, {{2, 1, 2}, {1, 1, 2}}, {{3/2, 1, 1/2}, {7/5, 2/3, 1/3}}}, {2, {{1, 1, 2}, {1, 1, 2}}, {{3, 1/5, 2}, {7/5, 2/3, 1/3}}}}}];
test["shifted family at s = 0 = unshifted family (FrequencyIntegration.shiftToZero runs on every shifted pattern)",
  ExactlyZeroQ[(FrequencyIntegralFamilies[2, {{2, 1, 1}, {1, 1, 1, "Shifted"}}] /. FrequencyIntegralShift[2] -> 0) - FrequencyIntegralFamilies[2, {{2, 1, 1}, {1, 1, 1}}]], True];
test["shifted family: regular at coincidence for s != 0 (Q data = q data)",
  FreeQ[Together[FrequencyIntegralFamilies[0, {{1, 0, 1}, {0, 1, 1, "Shifted"}}] /. famC[2, m_] :> famC[1, m]], Indeterminate | ComplexInfinity], True];
expectFail["a shifted result that differs from the unshifted one at s = 0 aborts (a factor 2: independent of s)",
  FRGscalarMSR`Private`loCheckShiftToZero[2 FrequencyIntegralFamilies[0, {{1, 0, 1}, {0, 1, 1, "Shifted"}}], 0, {{1, 0, 1}, {0, 1, 1, "Shifted"}}],
  "FrequencyIntegration.shiftToZero"];

(* NLO dispatch: even part, odd part / (I nu), exact; the odd part is finite at nu = 0 *)
fiNLO = (3 + nu w + nu^2) w/(\[CapitalDelta][w, s1]^2 \[CapitalDelta][-w, s1] \[CapitalDelta][w + nu, s2] \[CapitalDelta][-w - nu, s2]) +
   nu^3/(\[CapitalDelta][w + nu, s2]^2 \[CapitalDelta][-w - nu, s2]);
fiPart[part_] := FRGscalarMSR`Private`loIntegrate[fiNLO, w, {s1, s2}, 1, nu, part];
test["NLO dispatch: Even = (F(nu) + F(-nu))/2, I nu Odd = (F(nu) - F(-nu))/2, Odd finite at nu = 0",
  {ExactlyZeroQ[fiPart["Even"] - (fiPart["Full"] + (fiPart["Full"] /. nu -> -nu))/2], ExactlyZeroQ[I nu fiPart["Odd"] - (fiPart["Full"] - (fiPart["Full"] /. nu -> -nu))/2],
   ! ExactlyZeroQ[Denominator[Together[fiPart["Odd"]]] /. nu -> 0]}, {True, True, True}];
test["NLO dispatch: without nu the LO dispatch is unchanged", FRGscalarMSR`Private`loIntegrate[fiNLO /. nu -> 0, w, {s1, s2}, 1] ===
   FRGscalarMSR`Private`loIntegrate[fiNLO /. nu -> 0, w, {s1, s2}, 1, None, "Full"], True];
expectFail["NLO dispatch: Delta(w + 2 nu, Q^2) aborts", FRGscalarMSR`Private`loIntegrate[w/(\[CapitalDelta][w, s1]^2 \[CapitalDelta][w + 2 nu, s2] \[CapitalDelta][-w - 2 nu, s2]),
   w, {s1, s2}, 1, nu, "Full"], "FrequencyIntegration.ready"];
