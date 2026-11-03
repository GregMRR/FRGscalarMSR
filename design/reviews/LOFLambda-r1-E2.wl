(* the-enemy, LOFLAMBDA r1, E2: in d = 1 the KPZ time-reversal symmetry (hbar -> hbar + (nu/D) d_x^2 h, t -> -t, h -> -h) rests on
   Int dx h'' (h')^2 = 0, which makes the lam term invariant. With the f_lambda ansatz the lam term becomes
   -lam/2 hbar flam(-d_x^2)[(h')^2]; its variation under hbar -> hbar + c h'' is -c lam/2 Int h'' flam(-d_x^2)[(h')^2].
   Test: periodic h on [0, 2 Pi], flam(k^2) = 1 + a k^2 (first Taylor order of any flam).
   Usage: math -script design/reviews/LOFLambda-r1-E2.wl *)
h[x_] := Sin[x] + Cos[2 x]/3 + Sin[3 x]/5;
flamOp[g_, a_] := g + a (-D[g, {x, 2}]);
var[a_] := Integrate[D[h[x], {x, 2}] flamOp[D[h[x], x]^2, a], {x, 0, 2 Pi}];
Print["flam = 1:          Int h'' (h')^2           = ", Simplify[var[0]]];
Print["flam = 1 + a k^2:  Int h'' flam(-d^2)(h')^2 = ", Simplify[var[a]]];
Print["closed form -a Int (h'')^3 = ", Simplify[-a Integrate[D[h[x], {x, 2}]^3, {x, 0, 2 Pi}]]];
