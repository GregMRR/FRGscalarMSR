(* ::Package:: *)

(* Reference: the NLO loop integrands of I_D and I_nu with the external frequency w and Kloss-Canet-Wschebor's f_lambda slots, as
   kept in BMWflow flow/models/model_kpz.py Integral_pfixed_dD_NLO (f_lambdaq, f_lambdaQ, f_lambda_p; there set to 1; kernels A3a, A3c,
   A3d from kloss2012_omega_integration.nb; notes/NLO_flow_equations.tex), transcribed with f_lambda at p, q and Q (zero frequency).
   Variables of benchmarks/KPZLOLambda/Reference.wl plus w = the dimensionless external frequency (BMWflow's w, the library's aw).
   At w = 0 this is kpzlamRefLoop. No reference exists for the f_lambda row: it is checked by its w -> 0 limit (the KPZLOLambda row,
   itself checked against the derivative expansion), its p = 0 value and the oracle. *)

kpzlamnloRefLoop = Module[{q2 = qh^2, Q2 = ph^2 + qh^2 + 2 ph qh ch, lq, lQ, kq, kQ, flq, flQ, flp, fl, fl2, wff2, da, dc, dd, A3a, A3c, A3d,
    dsRD, dsRnu, pqc = ph qh ch},
   lq = q2 (fnu[q2] + r[q2]); lQ = Q2 (fnu[Q2] + r[Q2]); kq = fD[q2] + r[q2]; kQ = fD[Q2] + r[Q2];
   flq = flam[q2]; flQ = flam[Q2]; flp = flam[ph^2];
   fl = flq lQ + flQ lq; wff2 = (w flq flQ)^2;
   da = 2 lq lQ (fl^2 + wff2); A3a = fl/da;
   dc = da^2 lq/lQ; fl2 = fl + flQ lq; A3c = (fl^2 fl2 + wff2 flq lQ)/dc;
   dd = dc flq/lq^2; A3d = (fl^2 lQ + (w flQ)^2 fl2 flq)/dd;
   dsRD = -etaD r[q2] - 2 q2 r'[q2]; dsRnu = -etanu r[q2] - 2 q2 r'[q2];
   <|"D" -> 2 g flp^2 (q2 + pqc)^2 kQ (A3a dsRD - A3c dsRnu 2 q2 lq kq),
    "nu" -> -2 g flp/ph^2 (q2 + pqc) (-pqc flQ lQ A3a dsRD + (2 pqc flQ lQ lq kq A3c + (ph^2 + pqc) flq kQ (flq^2 A3d - lq^2 A3c)) q2 dsRnu)|>];
