(* ::Package:: *)

(* Reference at generic f_lambda: the LO loop integrands of I_D and I_nu with Kloss-Canet-Wschebor's f_lambda slots, as kept in
   BMWflow flow/models/model_kpz.py Integral_pfixed_dD_LO (f_lambdaq, f_lambdaQ, f_lambda_p; there set to 1; kernels A3a, A3c, A3d
   from kloss2012_omega_integration.nb), transcribed here with f_lambda at p, q and Q. Same variables as benchmarks/KPZLO/Reference.wl:
   ph, qh, ch = cos(theta), one shape r (r' = dr/d(q^2)), g = lambda^2 D/nu^3 kappa^(d-2). No reference exists for the f_lambda row
   itself: it is checked against the derivative expansion (Checks.wl, the KPZLOLambda.de checks). *)

kpzlamRefLoop = Module[{q2 = qh^2, Q2 = ph^2 + qh^2 + 2 ph qh ch, lq, lQ, kq, kQ, flq, flQ, flp, fl, fl2, A3a, A3c, A3d, dsRD, dsRnu, pqc = ph qh ch},
   lq = q2 (fnu[q2] + r[q2]); lQ = Q2 (fnu[Q2] + r[Q2]); kq = fD[q2] + r[q2]; kQ = fD[Q2] + r[Q2];
   flq = flam[q2]; flQ = flam[Q2]; flp = flam[ph^2];
   fl = flq lQ + flQ lq;
   A3a = fl/(2 lq lQ fl^2);
   fl2 = fl + flQ lq;
   A3c = fl^2 fl2/((2 lq lQ fl^2)^2 lq/lQ);
   A3d = fl^2 lQ/((2 lq lQ fl^2)^2 lq/lQ flq/lq^2);
   dsRD = -etaD r[q2] - 2 q2 r'[q2]; dsRnu = -etanu r[q2] - 2 q2 r'[q2];
   <|"D" -> 2 g flp^2 (q2 + pqc)^2 kQ (A3a dsRD - A3c dsRnu 2 q2 lq kq),
    "nu" -> -2 g flp/ph^2 (q2 + pqc) (-pqc flQ lQ A3a dsRD + (2 pqc flQ lQ lq kq A3c + (ph^2 + pqc) flq kQ (flq^2 A3d - lq^2 A3c)) q2 dsRnu)|>];
