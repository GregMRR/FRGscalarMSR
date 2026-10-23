(* ::Package:: *)

(* FrequencyIntegration.wl -- IntegrateFrequencies: every projected integrand is a sum of monomials w^n Delta(w)^-a Delta(-w)^-b;
   each is integrated over the real w-axis exactly by residues (closing in the upper half plane, where the roots of Delta(w) lie),
   with a result rational in the coefficients c0, c1, c2 of Delta (no square roots). Generalises Iint/omegaIntegrate of DE4plus_fast.wl.
   Objective, assumptions and checks: docs/MODULES.md section "FrequencyIntegration". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

IntegrateFrequencies::usage = "IntegrateFrequencies[projected] integrates every projected integrand over w. Declared monomials are \
integrated symbolically (result in GammaAtom, DeltaCoefficient and dsR atoms); the others (Protected, DeclaredNumeric, Causality, \
Neglected) are evaluated exactly at random rational points of all atoms. Returns the \"Integrated\" stage.";
FrequencyIntegral::usage = "FrequencyIntegral[n, a, b, deg] = Integrate[w^n / (Delta(w)^a Delta(-w)^b), {w, -Infinity, Infinity}] with \
Delta(w) = C0 + C1 (I w) + C2 (I w)^2 (deg = 2) or C0 + C1 (I w) (deg = 1), all roots of Delta(w) in the upper half plane; \
rational in the coefficients FrequencyIntegralC[0..2].";
FrequencyIntegralC::usage = "FrequencyIntegralC[m] is the coefficient c_m of Delta in FrequencyIntegral results.";
FrequencyIntegralFamilies::usage = "FrequencyIntegralFamilies[n, {{a1, b1}, {a2, b2}, ...}, deg] = Integrate[w^n / Prod_f (Delta_f(w)^a_f \
Delta_f(-w)^b_f), {w, -Infinity, Infinity}] for several propagator denominators Delta_f (same degree, every root of each Delta_f(w) in the \
upper half plane), rational in their coefficients FrequencyIntegralC[f, m]; the residues of all families are combined and cancelled, so \
the result is finite where two families coincide (design/MOMENTUMFUNCTIONS.md C3).";
ShowFlow::usage = "ShowFlow[integrated, {a, b}, {m, n}] prints the integrated flow of a monomial (symbolic, or its sampled values).";

Begin["`Private`"];

cc[m_] := FrequencyIntegralC[m];

(* residues at the roots z1, z2 of Delta(w) (upper half plane); z1,2 = I alpha +- beta, the result is even in beta *)
FrequencyIntegral[n_Integer, a_Integer, b_Integer, deg_Integer] := FrequencyIntegral[n, a, b, deg] = Block[{$FRGModule = "FrequencyIntegration"},
   Module[{w, z1, z2, z, f1, f2, r, al, be, x, num, den},
    FRGAssert["FrequencyIntegration.convergence", n - deg (a + b) <= -2, {n, a, b, deg},
     "divergent frequency integral: the integrand does not decay faster than 1/w^2"];
    Which[
     a == 0, 0,   (* no pole in the upper half plane *)
     deg == 1,
      z = I cc[0]/cc[1];
      f1 = w^n/((I cc[1])^a (-I cc[1])^b (w + z)^b);
      Together[2 Pi I (D[f1, {w, a - 1}] /. w -> z)/(a - 1)!],
     deg == 2,
      z1 = I al + be; z2 = I al - be;
      f1 = w^n/((-cc[2])^(a + b) (w - z2)^a ((w + z1) (w + z2))^b);
      f2 = w^n/((-cc[2])^(a + b) (w - z1)^a ((w + z1) (w + z2))^b);
      r = Together[2 Pi I ((D[f1, {w, a - 1}] /. w -> z1) + (D[f2, {w, a - 1}] /. w -> z2))/(a - 1)!];
      (* make the denominator even in beta, then beta^2 -> x *)
      num = Expand[Numerator[r] (Denominator[r] /. be -> -be)]; den = Expand[Denominator[r] (Denominator[r] /. be -> -be)];
      FRGAssert["FrequencyIntegration.evenInBeta", Expand[(num /. be -> -be) - num] === 0 && Expand[(den /. be -> -be) - den] === 0, {n, a, b},
       "internal: the residue sum is not symmetric in the two roots"];
      r = Together[((num /. be -> Sqrt[x]) /(den /. be -> Sqrt[x])) /. {x -> (4 cc[0] cc[2] - cc[1]^2)/(4 cc[2]^2), al -> cc[1]/(2 cc[2])}];
      FRGAssert["FrequencyIntegration.rational", FreeQ[r, be | x | al | Power[_, _Rational]], r, "internal: square roots survived"];
      r,
     True, FRGAssert["FrequencyIntegration.degree", False, deg, "deg Delta must be 1 or 2"]]]];

(* several families (generic p): residues at the upper-half-plane roots of every Delta_f(w); each family's roots enter symmetrically,
   so the sum is even in every beta_f and rational after beta_f^2 -> its value; Together cancels the factors that vanish where two
   families coincide (E3), checked by FrequencyIntegration.coincidentFinite *)
FrequencyIntegralFamilies[n_Integer, exps_List, deg_Integer] := FrequencyIntegralFamilies[n, exps, deg] = Block[{$FRGModule = "FrequencyIntegration"},
   Module[{nf = Length[exps], w, al, be, xx, roots, fac, f, res, num, den, cf},
    cf[i_, m_] := FrequencyIntegralC[i, m];
    FRGAssert["FrequencyIntegration.convergence", n - deg Total[Flatten[exps]] <= -2, {n, exps, deg},
     "divergent frequency integral: the integrand does not decay faster than 1/w^2"];
    (* Delta_f(w) = -c2 (w - z1)(w - z2) (deg 2), I c1 (w - z) (deg 1); Delta_f(-w) likewise with w -> -w *)
    roots[i_] := If[deg == 2, {I al[i] + be[i], I al[i] - be[i]}, {I cf[i, 0]/cf[i, 1]}];
    fac[i_, x_] := If[deg == 2, -cf[i, 2] (x - roots[i][[1]]) (x - roots[i][[2]]), I cf[i, 1] (x - roots[i][[1]])];
    f = w^n/Product[fac[i, w]^exps[[i, 1]] fac[i, -w]^exps[[i, 2]], {i, nf}];
    res = Sum[If[exps[[i, 1]] == 0, 0,
       Sum[With[{z = roots[i][[j]]}, Module[{g = Together[f (w - z)^exps[[i, 1]]]},
          (D[g, {w, exps[[i, 1]] - 1}] /. w -> z)/(exps[[i, 1]] - 1)!]], {j, Length[roots[i]]}]], {i, nf}];
    res = Together[2 Pi I res];
    If[deg == 2, Do[
       num = Expand[Numerator[res] (Denominator[res] /. be[i] -> -be[i])]; den = Expand[Denominator[res] (Denominator[res] /. be[i] -> -be[i])];
       FRGAssert["FrequencyIntegration.evenInBeta", Expand[(num /. be[i] -> -be[i]) - num] === 0 && Expand[(den /. be[i] -> -be[i]) - den] === 0, {n, exps},
        "internal: the residue sum is not symmetric in the two roots of one family"];
       res = Together[(num /. be[i] -> Sqrt[xx[i]])/(den /. be[i] -> Sqrt[xx[i]])], {i, nf}];
     res = Together[res /. Flatten[Table[{xx[i] -> (4 cf[i, 0] cf[i, 2] - cf[i, 1]^2)/(4 cf[i, 2]^2), al[i] -> cf[i, 1]/(2 cf[i, 2])}, {i, nf}]]]];
    FRGAssert["FrequencyIntegration.rational", FreeQ[res, be | xx | al | Power[_, _Rational]], res, "internal: square roots survived"];
    (* E3: finite where any two families coincide *)
    Do[FRGAssert["FrequencyIntegration.coincidentFinite", ! ExactlyZeroQ[Denominator[res] /. Table[cf[j, m] -> cf[i, m], {m, 0, 2}], "Seeds" -> {1}],
       {n, exps, {i, j}}, "internal: the combined residues are singular where two pole families coincide"], {i, nf}, {j, i + 1, nf}];
    res]];

(* dispatch on monomials w^n Prod_f X_f^a_f Y_f^b_f, X_f = 1/Delta(w, s_f), Y_f = 1/Delta(-w, s_f): one family per distinct momentum
   argument s_f; a monomial with one family uses FrequencyIntegral, several use FrequencyIntegralFamilies. coefOf[s] = {c0, c1, c2} at s.
   No table, nothing skipped *)
integrate[e_, w_, deg_, coefOf_] := Module[{args, exy, xs, ys, res, nf},
   args = DeleteDuplicates[Expand /@ Cases[e, \[CapitalDelta][w | -w, s_] :> s, {0, Infinity}]];
   nf = Length[args];
   xs = Array[famX, nf]; ys = Array[famY, nf];
   exy = e /. \[CapitalDelta][x_, s_] :> With[{i = FirstPosition[args, Expand[s], {0}, {1}][[1]]}, Which[i == 0, \[CapitalDelta][x, s], x === w, 1/xs[[i]], x === -w, 1/ys[[i]], True, \[CapitalDelta][x, s]]];
   FRGAssert["FrequencyIntegration.ready", FreeQ[exy, \[CapitalDelta] | Derivative[__][\[CapitalDelta]]], Cases[exy, _\[CapitalDelta], Infinity],
    "integrand not ready: Delta with another argument"];
   res = Total[Function[cr, Module[{ex = Transpose[{cr[[1, 2 ;; nf + 1]], cr[[1, nf + 2 ;;]]}], used},
        used = Flatten[Position[ex, {a_, b_} /; a + b > 0, {1}, Heads -> False]];
        cr[[2]] Which[
          used === {}, FRGAssert["FrequencyIntegration.convergence", False, cr, "a term without propagator: the frequency integral diverges"],
          Length[used] == 1, FrequencyIntegral[cr[[1, 1]], ex[[used[[1]], 1]], ex[[used[[1]], 2]], deg] /. Thread[cc /@ Range[0, 2] -> coefOf[args[[used[[1]]]]]],
          True, FrequencyIntegralFamilies[cr[[1, 1]], ex[[used]], deg] /.
            Flatten[Table[Thread[Table[FrequencyIntegralC[j, m], {m, 0, 2}] -> coefOf[args[[used[[j]]]]]], {j, Length[used]}]]]]] /@
      CoefficientRules[exy, Join[{w}, xs, ys]]];
   FRGAssert["FrequencyIntegration.complete", FreeQ[res, w | famX | famY | \[CapitalDelta] | Power[_, _Rational]], res, "leftover w, Delta or radicals"];
   res];

(* random rational point of every atom, Delta coefficients consistent with the atoms, w and Delta kept. The momentum invariants of the
   generic-p integrands are a physical configuration: q^2, p^2 > 0 and p.q = |p| |q| c with -1 < c < 1 *)
momentumInvariants = {LoopMomentumSquared, ExternalMomentumSquared, MomentumProduct};
momentumPoint[seed_] := Module[{qm, pm, cs},
   SeedRandom[seed + 31337];
   qm = RandomInteger[{1, 40}]/RandomInteger[{10, 40}]; pm = RandomInteger[{1, 40}]/RandomInteger[{10, 40}]; cs = RandomInteger[{-19, 19}]/20;
   {LoopMomentumSquared -> qm^2, ExternalMomentumSquared -> pm^2, MomentumProduct -> pm qm cs}];
sampleAtoms = x : (Derivative[__][GammaAtom[__]][_] | GammaAtom[__][_] | Derivative[__][dsR[__]][_] | dsR[__][_] | dsR[__] | dsCoupling[_]);

samplePoint[e_, prop_, w_, seed_] := Module[{e1, syms, pt, memo, mom},
   e1 = e /. {DeltaCoefficient[m_][s_] :> prop["DeltaCoefficients"][s][[m + 1]],
       Derivative[k_][DeltaCoefficient[m_]][s_] :> (D[prop["DeltaCoefficients"][xs][[m + 1]], {xs, k}] /. xs -> s)} /. \[CapitalDelta] -> dd;
   syms = DeleteCases[Complement[Union[Cases[e1, x_Symbol /; Context[x] =!= "System`", {0, Infinity}, Heads -> False]], momentumInvariants], w];
   pt = RandomRationalPoint[syms, seed];
   mom = momentumPoint[seed];
   SeedRandom[seed + 104729];
   memo[x_] := memo[x] = RandomInteger[{2, 60}]/RandomInteger[{61, 97}];
   e1 = e1 /. mom /. pt /. sampleAtoms :> memo[x];
   {e1 /. dd -> \[CapitalDelta], Join[pt, mom], memo}];

(* {value, valid}: valid = every momentum argument of a Delta is > 0 and satisfies the Propagator validity conditions there *)
sampleValue[e_, prop_, w_, deg_, seed_] := Module[{e1, pt, memo, coefOf, args},
   {e1, pt, memo} = samplePoint[e, prop, w, seed];
   coefOf = Function[s, prop["DeltaCoefficients"][s] /. x : (GammaAtom[__][_]) :> memo[x]];
   args = DeleteDuplicates[Cases[e1, \[CapitalDelta][_, s_] :> s, {0, Infinity}]];
   {integrate[e1, w, deg, coefOf], AllTrue[args, validPointQ[prop, #, coefOf[#]] &]}];

(* the point is physical: s > 0 and every root of Delta(w) in the upper half plane (Propagator validity conditions) *)
validPointQ[prop_, s_, cv_] := TrueQ[s > 0] && TrueQ[And @@ (prop["ValidityConditions"][s] /. DeltaCoefficient[m_][_] :> cv[[m + 1]])];

(* a sample at a valid point, for values that are reported as sizes (Neglected); zero tests do not need valid points *)
validSample[e_, prop_, w_, deg_, seed_] := Module[{r},
   Do[r = sampleValue[e, prop, w, deg, seed + 1009 j]; If[r[[2]], Return[r[[1]], Module]], {j, 0, 199}];
   FRGAssert["FrequencyIntegration.validPoint", False, e, "no random point satisfies the validity conditions in 200 tries: check the signs of the ansatz"]];

symbolicCoefficients = Function[s, DeltaCoefficient[#][s] & /@ Range[0, 2]];

IntegrateFrequencies[projected_] := Block[{$FRGModule = "FrequencyIntegration"},
  Module[{t0 = AbsoluteTime[], prop, w, deg, res, checks = {}, warns = {}, stage, nI, nF, zeroSamples},
   validateStage[projected, "Projection"];
   prop = projected["Loops"]["Propagator"]; w = projected["LoopFrequency"]; deg = prop["DeltaDegree"];
   zeroSamples[ig_] := Table[First[sampleValue[ig, prop, w, deg, sd]], {sd, {1, 2, 3}}];
   res = Function[r, Module[{t1 = AbsoluteTime[], out},
        out = Switch[r["Class"],
          "Declared", <|"Value" -> integrate[r["Integrand"], w, deg, symbolicCoefficients]|>,
          "Functional", Join[<|"Value" -> integrate[r["Integrand"], w, deg, symbolicCoefficients],
             "AtZeroValue" -> integrate[r["AtZero"], w, deg, symbolicCoefficients]|>,
            If[r["ProtectedAtZero"] === None, <||>, <|"AtZeroSamples" -> zeroSamples[r["AtZero"]]|>]],
          "Neglected", <|"Samples" -> {validSample[r["Integrand"], prop, w, deg, 1]}, "SampleValid" -> True|>,
          _, <|"Samples" -> zeroSamples[r["Integrand"]]|>];
        Join[KeyDrop[r, {"Integrand", "AtZero"}], out, <|"Seconds" -> AbsoluteTime[] - t1|>]]] /@ projected["Integrands"];
   nI = Length[Select[DownValues[FrequencyIntegral], FreeQ[#, Pattern] &]];
   nF = Length[Select[DownValues[FrequencyIntegralFamilies], FreeQ[#, Pattern] &]];
   AppendTo[checks, assertCheck["FrequencyIntegration.allIntegrated", True, Null, "",
     ToString[Count[res, KeyValuePattern["Class" -> "Declared" | "Functional"]]] <> " symbolic, " <> ToString[Count[res, _?(KeyExistsQ[#, "Samples"] &)]] <>
      " sampled; no leftover w, Delta or radicals; " <> ToString[nI] <> " distinct residue integrals" <>
      If[nF > 0, ", " <> ToString[nF] <> " with several pole families (combined residues, finite where families coincide)", ""]]];
   AppendTo[warns, infoCheck["INFO", "FrequencyIntegration.samples", "zero tests use 3 arbitrary random points; Neglected values use one point where q^2 > 0 and the validity conditions hold"]];
   AppendTo[warns, infoCheck["ASSUMPTION", "FrequencyIntegration.poles", "all roots of Delta(w) in the upper half plane (Propagator validity conditions)" <>
      If[nF > 0, ", at every momentum argument (|q -+ p/2|^2 and q^2 for the generic-p flows)", ""]]];
   stage = stampStage[<|"Stage" -> "Integrated", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {projected["Fingerprint"]},
      "Projection" -> projected, "Flows" -> res, "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

ShowFlow[it_Association, {a_, b_}, {m_, n_}] := Module[{r = SelectFirst[it["Flows"], #Pair === {a, b} && #Frequency == m && #Momentum == n &]},
   If[MissingQ[r], Print["no such monomial"],
    Print["d_s of (I nu)^", m, " (p^2)^", n, " in Gamma2[", a, ", ", b, "] (", r["Class"], "): ",
     If[KeyExistsQ[r, "Value"], displayString[r["Value"]], "samples " <> inputForm[r["Samples"]]]]]];

End[];
EndPackage[];
