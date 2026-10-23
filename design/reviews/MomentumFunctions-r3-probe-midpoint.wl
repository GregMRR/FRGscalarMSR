(* R1: is the library's KPZ loop integrand even in p at fixed midpoint r = q + p/2 (d = 1, c = +1, nu = 0)? *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
li = LoopIntegrand[kpzLoops, {hbar, h}, 0, pp];
li = li /. DsRRules[kpzRegulator, scal[qv, qv]];
li = li //. DeltaRules[kpzPropagator];
atoms = {GammaAtom[h, hbar, 0] -> (# /(1 + #) + # Exp[-#] &), GammaAtom[hbar, h, 0] -> (# /(1 + #) + # Exp[-#] &),
   GammaAtom[hbar, hbar, 0] -> (-2/(1 + 2 #) - 2 Exp[-#] &)};
li = li /. atoms /. AtomRules[kpzPropagator];
num[r_, p_] := li /. scal[a_, b_] :> ((a /. {qv -> r - p/2, pp -> p}) (b /. {qv -> r - p/2, pp -> p})) /.
   {rn -> (Exp[-#] &), rd -> (Exp[-#] &), nuk -> 1, Dk -> 1, lamk -> 1, k -> 1, dsCoupling[_] -> 0};
Print["free symbols left: ", Complement[Cases[num[7/10, 1/10], _Symbol, Infinity], {wv, Pi, E}]];
wint[r_, p_] := With[{w = wv, e = num[r, p]}, NIntegrate[e, {w, -Infinity, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 12]];
Do[With[{a = wint[7/10, p], b = wint[7/10, -p]},
   Print["p = ", N[p], "  midpoint raw I/p^2 = ", N[a/p^2, 6], "   (I(p)+I(-p))/2/p^2 = ", N[(a + b)/2/p^2, 6], "   odd part (I(p)-I(-p))/2/p = ", N[(a - b)/2/p, 6]]],
 {p, {1/10, 1/100, 1/1000, 1/10000}}];
