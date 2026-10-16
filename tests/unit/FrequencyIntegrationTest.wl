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
