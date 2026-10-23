(* R1 + R4: which shape derivatives does the p^2 Taylor coefficient of the KPZ (hbar,h) loop integrand contain,
   at fixed q (today's routing) vs at fixed midpoint r (symmetrised)? d = 1, c = +1, nu = 0, before the w-integral *)
Get["/home/greg/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
quiet[x_] := Block[{Print}, x];
Get["/home/greg/Physics/Sources/FRGscalarMSR/tests/reference/KPZSetup.wl"] // quiet;
{wv, qv} = LoopVariables[kpzLoops];
li = LoopIntegrand[kpzLoops, {hbar, h}, 0, pp] /. DsRRules[kpzRegulator, scal[qv, qv]] //. DeltaRules[kpzPropagator];
li = li /. AtomRules[kpzPropagator] /. k -> 1;
at[q_, p_] := li /. scal[a_, b_] :> ((a /. {qv -> q, pp -> p}) (b /. {qv -> q, pp -> p}));
shapeDerivs[e_] := Union[Cases[e, Derivative[n_][rn | rd] :> n, {0, Infinity}, Heads -> True]];
today = D[at[x, p], {p, 2}] /. p -> 0;
mid = D[(at[x - p/2, p] + at[-x - p/2, p])/2, {p, 2}] /. p -> 0;
Print["orders of rn, rd derivatives in the p^2 coefficient, today's routing: ", shapeDerivs[today]];
Print["orders of rn, rd derivatives in the p^2 coefficient, symmetrised midpoint: ", shapeDerivs[mid]];
Print[Union[Cases[at[x, p], (rn | rd | Derivative[_][rn | rd])[_], {0, Infinity}, Heads -> True]]];
Print[Union[Cases[li, _GammaAtom | _DeltaCoefficient | _\[CapitalDelta], {0, Infinity}, Heads -> True]]];
