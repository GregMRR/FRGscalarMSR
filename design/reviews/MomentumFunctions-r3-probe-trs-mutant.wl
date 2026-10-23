(* R2: is the d = 1 time-reversal identity of the Gamma2 flows pointwise in the loop momentum at generic p (constant-coupling KPZ,
   rn = rd, eta_nu = eta_D)? I11/p^2 - (nu/D) (-I02/2) at fixed q, raw and q <-> -p-q symmetrised, and after the q-integral *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
prep[pair_] := Module[{li = LoopIntegrand[kpzLoops, pair, 0, pp]},
   li = li /. DsRRules[kpzRegulator, scal[qv, qv]] //. DeltaRules[kpzPropagator];
   li /. AtomRules[kpzPropagator]];
nuv = 2; Dv = 1; eta = 3/10;
pars = {rn -> (Exp[-#] &), rd -> (Exp[-#] &), nuk -> nuv, Dk -> Dv, lamk -> 1, k -> 1, dsCoupling[nuk] -> -eta nuv, dsCoupling[Dk] -> -(eta+1/5) Dv, dsCoupling[lamk] -> 0};
l11 = prep[{hbar, h}] /. pars; l02 = prep[{hbar, hbar}] /. pars;
num[e_, q_, p_] := e /. scal[a_, b_] :> ((a /. {qv -> q, pp -> p}) (b /. {qv -> q, pp -> p}));
wint[e_] := With[{w = wv}, NIntegrate[e, {w, -Infinity, Infinity}, WorkingPrecision -> 25, PrecisionGoal -> 12]];
f[q_?NumericQ, p_?NumericQ] := wint[num[l11, q, p]]/p^2 - (nuv/Dv) (-wint[num[l02, q, p]]/2);
fs[q_?NumericQ, p_?NumericQ] := (f[q, p] + f[-p - q, p])/2;
p0 = 3/5;
Print["pointwise raw  at q = 0.7, 1.3: ", N[{f[7/10, p0], f[13/10, p0]}, 6]];
Print["pointwise sym  at q = 0.7, 1.3: ", N[{fs[7/10, p0], fs[13/10, p0]}, 6]];
