(* ::Package:: *)

(* Reference flows: L. Gosteva, N. Wschebor, L. Canet, "Scaling regimes of the Kuramoto-Sivashinsky equation from the functional
   renormalization group", arXiv:2607.15784 (2026), Sec. III B, Eqs. (34)-(37), transcribed from the paper; the transcription was
   checked to equal the authors' notebook DE/supmat.nb (gitlab.com/goste/kuramoto-sivashinsky-lpa, "before integration by parts")
   exactly.

   Each ksRef[c] is the integrand J(qh) of Int_{-Inf}^{Inf} dqh/(2 Pi) J(qh), qh = q/k, for
     ksRef["D"]    = d_s D / D          (34)        ksRef["zeta"] = k^2 d_s zeta / D   (36)
     ksRef["nu"]   = d_s nu / D         (35)        ksRef["tau"]  = k^2 d_s tau / D    (37)
   in the dimensionless variables of Eq. (25): Gh = lam^2/(k D^2) (the paper's D_Lambda is a unit, set to 1), Fh = nu/D,
   Ch = k^2 zeta/D, Th = k^2 tau/D, and eta = -d_s D/D (30). Notation below Eq. (37), with x = qh^2 and ' = d/dx:
     hn = Fh + Th x + rnu(x),  hd = 1 + Ch x + rD(x),  sn = d_s R^nu/D,  sd = d_s R^D/D.
   (The paper writes sn = d_s R^nu/nu; with R^nu = D rnu, Eq. (31), its own sn = alpha sd = -alpha (eta r + 2 x r') of the
   alpha-regulator requires d_s R^nu/D, which is what the comparison confirms.)
   For the regulators of Eq. (31), R = D r(x): sn = -(eta rnu + 2 x rnu'), sd = -(eta rD + 2 x rD') (ksSRules). *)

ksRefTerms[x_] := Module[{hn = Fh + Th x + rnu[x], hd = 1 + Ch x + rD[x], hn1 = Th + rnu'[x], hd1 = Ch + rD'[x],
    rn2 = rnu''[x], rd2 = rD''[x], rn3 = rnu'''[x], rd3 = rD'''[x], sn = snu[x], sd = sD[x]},
   <|
    "D" -> -Gh hd (3 sn hd - 2 sd hn)/(4 x hn^4),
    "nu" -> -Gh (2 x (sd hn1 - sn hd1) + hd sn)/(4 x hn^3),
    "zeta" -> Gh/(4 x^2 hn^6) (
       hd hn (hn x (5 sn hd1 - 6 x (rd2 sn + rn2 sd) + hn1 sd) + 2 x^2 hn1 (10 sn hd1 + 7 sd hn1) + hn^2 sd)
       - 2 hn^2 x sd (hd1 (hn + 6 x hn1) - 2 x hn rd2)
       - hd^2 sn (3 hn^2 + 5 hn x (hn1 - 2 x rn2) + 25 x^2 hn1^2)),
    "tau" -> Gh/(12 hn^6 x^2) (
       hn^3 (x (6 hd1 sn + 3 hn1 sd + 4 x^2 (rd3 sn - rn3 sd) - 6 x (rd2 sn + 2 rn2 sd)) - 9 hd sn)
       + 6 hn^2 x (2 hn1 (sn (hd1 x - hd) + x^2 (rn2 sd - 2 rd2 sn)) + 4 x rn2 sn (hd - hd1 x) + 3 hn1^2 x sd)
       + 3 hn hn1 x^2 (hn1 sn (22 hd1 x - 13 hd) + 12 hd x rn2 sn - 2 hn1^2 x sd)
       - 60 hd hn1^3 x^3 sn + 3 hn^4 sd)
   |>];

ksRef[c_String] := ksRefTerms[qh^2][c];

(* s_nu, s_D for R^nu = D rnu, R^D = D rD *)
ksSRules = {snu -> Function[x, -(eta rnu[x] + 2 x rnu'[x])], sD -> Function[x, -(eta rD[x] + 2 x rD'[x])]};
