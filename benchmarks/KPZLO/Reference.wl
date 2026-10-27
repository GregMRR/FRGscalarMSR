(* ::Package:: *)

(* Reference: KPZ LO flow as written in BMWflow notes/LO_flow_equations.tex + notes/kpz_notation.tex (Kloss, Canet, Wschebor,
   PRE 86, 051124 (2012), f_lambda = 1, omega = 0) and coded by hand in BMWflow flow/models/model_kpz.py (Integral_pfixed_dD_LO,
   eta_calc, g_rhs_calc), transcribed here. Dimensionless: p, q in units of kappa, c = cos(theta), one regulator shape r with
   r' = dr/d(q^2), primes on f = d/dq, g = lambda^2 D/nu^3 kappa^(d-2).

   kpzloRefLoop[X]: the integrand F_X of I_X(p) = v_{d-1}/(2 Pi) Int dq q^(d-1) Int dtheta sin^(d-2)(theta) F_X (= Int d^dq/(2 Pi)^d F_X),
     I_D (eq. for I_D of LO_flow_equations.tex without the prefactor's measure), I_nu (with its 1/p^2);
   kpzloRefEta: the integrands of the eta system (1 + I) eta = -I0 with <.> = v_d Int dq q^(d-1) (kpz_notation.tex),
     {{I_DD, I_Dnu}, {I_nuD, I_nunu}} and {I_D0, I_nu0} without v_d and q^(d-1);
   kpzloRefLinear[X] = eta_X f_X + p f_X'(p);  d_s g = g (d - 2 - eta_D + 3 eta_nu). *)

kpzloRefLoop = Module[{q2 = qh^2, Q2 = ph^2 + qh^2 + 2 ph qh ch, lq, lQ, kq, kQ, L, Aa, Ac, Ad, dsRD, dsRnu, pqc = ph qh ch},
   lq = q2 (fnu[q2] + r[q2]); lQ = Q2 (fnu[Q2] + r[Q2]); kq = fD[q2] + r[q2]; kQ = fD[Q2] + r[Q2]; L = lq + lQ;
   Aa = 1/(2 lq lQ L); Ac = (lQ + 2 lq)/(4 lq^3 lQ L^2); Ad = 1/(4 lq L^2);
   dsRD = -etaD r[q2] - 2 q2 r'[q2]; dsRnu = -etanu r[q2] - 2 q2 r'[q2];
   <|"D" -> 2 g (q2 + pqc)^2 kQ (Aa dsRD - 2 q2 lq kq Ac dsRnu),
    "nu" -> -2 g/ph^2 (q2 + pqc) (-pqc lQ Aa dsRD + (2 pqc lQ lq kq Ac + (ph^2 + pqc) kQ (Ad - lq^2 Ac)) q2 dsRnu)|>];

(* f at q as a function of q^2 inside the library: f'(q) = d f/dq = 2 q F'(q^2) *)
kpzloRefEta = Module[{q2 = qh^2, kk, l, qdl, qdk, rr = r[qh^2], rr1 = r'[qh^2], d = Dim, fDp = 2 qh fD'[qh^2], fnup = 2 qh fnu'[qh^2]},
   kk = fD[q2] + rr; l = q2 (fnu[q2] + rr);
   qdl = 2 l + q2 (qh fnup + 2 q2 rr1); qdk = qh fDp + 2 q2 rr1;
   <|"Matrix" -> {{-g/2 qh^4 rr kk/l^3, 3 g/4 qh^6 rr kk^2/l^4}, {g/(4 d) qh^2 rr (qdl - 2 l)/l^3, -g/(4 d) qh^4 rr (qdk - (2 - d) kk)/l^3}},
    "Vector" -> {g/2 qh^6 rr1 kk (3 q2 kk - 2 l)/l^4, -g/(2 d) qh^4 rr1 (q2 qdk - qdl + (d - 2) q2 kk + 2 l)/l^3}|>];

kpzloRefLinear = <|"D" -> etaD fD[ph^2] + ph fDprime, "nu" -> etanu fnu[ph^2] + ph fnuprime|>;
kpzloRefGFlow = g (Dim - 2 - etaD + 3 etanu);
