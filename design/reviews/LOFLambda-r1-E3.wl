(* the-enemy, LOFLAMBDA r1, E3: BMWflow's ModelKPZ LO kernel keeps Kloss-2012 f_lambda placeholders (model_kpz.py:78-101,
   eta_calc :193-204). Spot check that they are a usable reference at generic f_lambda: the elementary frequency integral behind
   A3a = fl/denom_a = 1/(2 lq lQ (fq lQ + fQ lq)) equals Int dw/(2 Pi) 1/((w^2 fq^2 + lq^2)(w^2 fQ^2 + lQ^2)) for the KPZ LO
   propagator with Gamma^(1,1) = i w f_lambda + l.
   Usage: math -script design/reviews/LOFLambda-r1-E3.wl *)
lhs = Integrate[1/((w^2 fq^2 + lq^2) (w^2 fQ^2 + lQ^2)), {w, -Infinity, Infinity},
    Assumptions -> fq > 0 && fQ > 0 && lq > 0 && lQ > 0]/(2 Pi);
fl = fq lQ + fQ lq; a3a = fl/(2 lq lQ fl^2);
Print["Int dw/2Pi C C / (4 kq kQ) - A3a = ", Simplify[lhs - a3a]];
