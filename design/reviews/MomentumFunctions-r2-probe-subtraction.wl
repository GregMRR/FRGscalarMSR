(* r2-probe-subtraction: is [flow(p) - sum_{n<n0} flow_n p^2n] / p^(2 n0) O(1) pointwise after q <-> -p-q symmetrisation?
   Generic rotation-invariant integrand I(p,q) = f(q^2, |p+q|^2) + g(q^2) (p.q)^2 (not symmetric), d symbolic.
   Lower orders subtracted as the plan says: today's operator (n=0: value at p=0; n=1: direction-averaged s^2 coefficient),
   the subtracted integrand L(p,q) symmetrised too (best case). *)
f[a_, b_] := 1/((1 + a) (2 + b)^2); g[a_] := 1/(3 + a);
qq = 7/10;                                  (* |q| *)
P2[s_, x_] := s^2 + qq^2 + 2 s qq x;        (* |p+q|^2, x = cos theta *)
I1[s_, x_] := f[qq^2, P2[s, x]] + g[qq^2] (s qq x)^2;
(* q -> -p-q : q^2 -> |p+q|^2, |p+q|^2 -> q^2, p.q -> -p^2 - p.q *)
I2[s_, x_] := f[P2[s, x], qq^2] + g[P2[s, x]] (-s^2 - s qq x)^2;
sym[s_, x_] := (I1[s, x] + I2[s, x])/2;
(* lower orders with today's operator, as functions of the magnitude Q of the loop momentum *)
c0F[Q2_] := f[Q2, Q2];
c2F[Q2_] := Module[{e = f[Q2, sg^2 + Q2 + 2 sg Sqrt[Q2] xx] + g[Q2] (sg Sqrt[Q2] xx)^2},
   Expand[SeriesCoefficient[e, {sg, 0, 2}]] /. xx^2 -> 1/d /. xx -> 0];
Lfun[Q_, s_] := c0F[Q^2] + s^2 c2F[Q^2];
Lsym[s_, x_] := (Lfun[qq, s] + Lfun[Sqrt[P2[s, x]], s])/2;
(* case n0 = 1 (subtract p^0 only, symmetrised) *)
rem1 = Series[sym[s, x] - (Lfun[qq, 0] + Lfun[Sqrt[P2[s, x]], 0])/2, {s, 0, 2}];
Print["n0=1  coefficient of s^1 (must be 0): ", Simplify[SeriesCoefficient[rem1, 1]]];
(* case n0 = 2 (subtract p^0 and p^2, symmetrised) *)
rem2 = Series[sym[s, x] - Lsym[s, x], {s, 0, 3}];
Print["n0=2  s^1: ", Simplify[SeriesCoefficient[rem2, 1]]];
Print["n0=2  s^2 (must be 0 for O(p^4)): ", Simplify[SeriesCoefficient[rem2, 2]]];
Print["n0=2  s^2 at d=3, x=1 and x=0: ", N[SeriesCoefficient[rem2, 2] /. d -> 3 /. x -> 1], "  ", N[SeriesCoefficient[rem2, 2] /. d -> 3 /. x -> 0]];
Print["n0=2  s^3: ", N[SeriesCoefficient[rem2, 3] /. d -> 3 /. x -> 1/2]];
Print["n0=2  angular average at fixed aq of the s^2 remainder (nonzero: cancels only over the whole q-integral): ",
  Simplify[Expand[Expand[Normal[SeriesCoefficient[rem2, 2]] ] /. x^4 -> 3/(d (d + 2)) /. x^3 -> 0 /. x^2 -> 1/d /. x -> 0]]];
