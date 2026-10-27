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
FrequencyIntegralC::usage = "FrequencyIntegralC[m] is the coefficient c_m of Delta in FrequencyIntegral results; FrequencyIntegralC[j, m] \
the coefficient c_m of the Delta of pole family j in FrequencyIntegralFamilies results.";
FrequencyIntegralFamilies::usage = "FrequencyIntegralFamilies[n, {{a1, b1, deg}, {a2, b2, deg}}] = Integrate[w^n / Prod_j (Delta_j(w)^aj \
Delta_j(-w)^bj), {w, -Infinity, Infinity}] for two pole families Delta_j(w) = sum_m FrequencyIntegralC[j, m] (I w)^m (the propagator \
Delta at q and at Q = p + q, external frequency 0), all roots of each Delta_j(w) in the upper half plane; rational in the coefficients \
(LO pipeline, design/LOPIPELINE.md step 4).";
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

(* ---------- two pole families (LO): Delta(w, q^2) and Delta(w, Q^2), external frequency 0 ---------- *)

(* the roots of Delta_j(w) (upper half plane) and the leading coefficient: Delta_j(w) = lead (w - z_1)...; Delta_j(-w) = lead (-1)^deg (w + z_1)... *)
loFamilyRoots[j_, 1] := {{I FrequencyIntegralC[j, 0]/FrequencyIntegralC[j, 1]}, I FrequencyIntegralC[j, 1]};
loFamilyRoots[j_, 2] := {{I loAl[j] + loBe[j], I loAl[j] - loBe[j]}, -FrequencyIntegralC[j, 2]};

(* r is even in loBe[j] separately for each family: make the denominator even, then beta^2 and alpha in the coefficients *)
loEliminateBeta[r0_, j_] := Module[{be = loBe[j], num, den, x},
   num = Expand[Numerator[r0] (Denominator[r0] /. be -> -be)]; den = Expand[Denominator[r0] (Denominator[r0] /. be -> -be)];
   FRGAssert["FrequencyIntegration.evenInBeta", Expand[(num /. be -> -be) - num] === 0 && Expand[(den /. be -> -be) - den] === 0, j,
    "internal: the residue sum is not symmetric in the two roots of one pole family"];
   Together[((num /. be -> Sqrt[x]) /(den /. be -> Sqrt[x])) /. {x -> (4 FrequencyIntegralC[j, 0] FrequencyIntegralC[j, 2] - FrequencyIntegralC[j, 1]^2)/(4 FrequencyIntegralC[j, 2]^2),
      loAl[j] -> FrequencyIntegralC[j, 1]/(2 FrequencyIntegralC[j, 2])}]];

(* the residue sum is regular when the two families coincide (Q data -> q data): single residues carry 1/(z_q - z_Q) factors that
   must cancel in the Together-simplified sum (design/LOPIPELINE.md step 4) *)
loCheckRegularAtCoincidence[r_] := Module[{den = Denominator[Together[r]], co},
   co = den /. FrequencyIntegralC[2, m_] :> FrequencyIntegralC[1, m];
   FRGAssert["FrequencyIntegration.regularAtCoincidence", ! ExactlyZeroQ[co], den,
    "internal: the two-family residue sum is singular when the Delta at Q equals the Delta at q (uncancelled 1/(z_q - z_Q))"];
   r];

FrequencyIntegralFamilies[n_Integer, fams : {{_Integer, _Integer, _Integer} ..}] := FrequencyIntegralFamilies[n, fams] = Block[{$FRGModule = "FrequencyIntegration"},
   Module[{w, roots, leads, f, r = 0, nf = Length[fams]},
    FRGAssert["FrequencyIntegration.convergence", n - Total[(#[[1]] + #[[2]]) #[[3]] & /@ fams] <= -2, {n, fams},
     "divergent frequency integral: the integrand does not decay faster than 1/w^2"];
    FRGAssert["FrequencyIntegration.degree", AllTrue[fams[[All, 3]], MemberQ[{1, 2}, #] &], fams, "deg Delta must be 1 or 2 per pole family"];
    {roots, leads} = Transpose[Table[loFamilyRoots[j, fams[[j, 3]]], {j, nf}]];
    f = w^n/Product[(leads[[j]]^fams[[j, 1]] Product[(w - z)^fams[[j, 1]], {z, roots[[j]]}]) *
        ((leads[[j]] (-1)^fams[[j, 3]])^fams[[j, 2]] Product[(w + z)^fams[[j, 2]], {z, roots[[j]]}]), {j, nf}];
    (* residues at the upper-half-plane roots of every Delta_j(w), pole order a_j *)
    Do[With[{a = fams[[j, 1]]}, If[a > 0, Do[
        r += (D[Together[f (w - z)^a], {w, a - 1}] /. w -> z)/(a - 1)!, {z, roots[[j]]}]]], {j, nf}];
    r = Together[2 Pi I r];
    Do[If[fams[[j, 3]] == 2, r = loEliminateBeta[r, j]], {j, nf}];
    FRGAssert["FrequencyIntegration.rational", FreeQ[r, _loBe | _loAl | Power[_, _Rational]], r, "internal: square roots survived"];
    loCheckRegularAtCoincidence[r]]];

(* LO dispatch: monomials w^n X1^a1 Y1^b1 X2^a2 Y2^b2 with X_j = 1/Delta(w, s_j), Y_j = 1/Delta(-w, s_j); a single family goes to
   FrequencyIntegral. cvals[j] are the rules FrequencyIntegralC[m] -> c_m of family j *)
loIntegrate[e_, w_, {s1_, s2_}, deg_] := Module[{exy, res, cv},
   exy = e /. {\[CapitalDelta][w, s1] -> 1/loX1, \[CapitalDelta][-w, s1] -> 1/loY1, \[CapitalDelta][w, s2] -> 1/loX2, \[CapitalDelta][-w, s2] -> 1/loY2};
   FRGAssert["FrequencyIntegration.ready", FreeQ[exy, \[CapitalDelta] | Derivative[__][\[CapitalDelta]]], Cases[exy, _\[CapitalDelta], Infinity],
    "integrand not ready: Delta with a momentum other than q and Q"];
   cv[j_, s_] := Table[FrequencyIntegralC[j, m] -> DeltaCoefficient[m][s], {m, 0, 2}];
   res = Total[Function[cr, Module[{n, a1, b1, a2, b2},
        {n, a1, b1, a2, b2} = cr[[1]];
        cr[[2]] Which[
          a2 == 0 && b2 == 0, FrequencyIntegral[n, a1, b1, deg] /. FrequencyIntegralC[m_] :> DeltaCoefficient[m][s1],
          a1 == 0 && b1 == 0, FrequencyIntegral[n, a2, b2, deg] /. FrequencyIntegralC[m_] :> DeltaCoefficient[m][s2],
          True, FrequencyIntegralFamilies[n, {{a1, b1, deg}, {a2, b2, deg}}] /. Join[cv[1, s1], cv[2, s2]]]]] /@
      CoefficientRules[exy, {w, loX1, loY1, loX2, loY2}]];
   FRGAssert["FrequencyIntegration.complete", FreeQ[res, w | loX1 | loY1 | loX2 | loY2 | \[CapitalDelta] | Power[_, _Rational]], res, "leftover w, Delta or radicals"];
   res];

(* dispatch on monomials w^n X^a Y^b, X = 1/Delta(w), Y = 1/Delta(-w): no table, nothing skipped *)
integrate[e_, w_, deg_, cvals_] := Module[{exy, res},
   (* one momentum only: a loop with Delta(w, q^2) and Delta(w, Q^2) needs the two-family residues (review E4 of design/LOPIPELINE.md) *)
   With[{moms = Union[Cases[e, \[CapitalDelta][w | -w, s_] :> s, {0, Infinity}]]},
    FRGAssert["FrequencyIntegration.singleMomentum", Length[moms] <= 1, moms,
     "Delta appears at several momenta: every pole factor would get the same residue (wrong); this needs FrequencyIntegralFamilies (LO pipeline)"]];
   exy = e /. {\[CapitalDelta][w, _] -> 1/xX, \[CapitalDelta][-w, _] -> 1/yY};
   FRGAssert["FrequencyIntegration.ready", FreeQ[exy, \[CapitalDelta] | Derivative[__][\[CapitalDelta]]], Cases[exy, _\[CapitalDelta], Infinity],
    "integrand not ready: Delta with another argument"];
   res = Total[#[[2]] (FrequencyIntegral[#[[1, 1]], #[[1, 2]], #[[1, 3]], deg] /. cvals) & /@ CoefficientRules[exy, {w, xX, yY}]];
   FRGAssert["FrequencyIntegration.complete", FreeQ[res, w | xX | yY | \[CapitalDelta] | Power[_, _Rational]], res, "leftover w, Delta or radicals"];
   res];

(* random rational point of every atom, Delta coefficients consistent with the atoms, w and Delta kept *)
samplePoint[e_, prop_, w_, seed_] := Module[{e1, syms, atoms, pt, memo},
   e1 = e /. Join[{}, {DeltaCoefficient[m_][s_] :> prop["DeltaCoefficients"][s][[m + 1]],
       Derivative[k_][DeltaCoefficient[m_]][s_] :> (D[prop["DeltaCoefficients"][xs][[m + 1]], {xs, k}] /. xs -> s)}] /. \[CapitalDelta] -> dd;
   syms = DeleteCases[Union[Cases[e1, x_Symbol /; Context[x] =!= "System`", {0, Infinity}, Heads -> False], {LoopMomentumSquared}], w];
   pt = RandomRationalPoint[syms, seed];
   SeedRandom[seed + 104729];
   memo[x_] := memo[x] = RandomInteger[{2, 60}]/RandomInteger[{61, 97}];
   e1 = e1 /. pt /. x : (Derivative[__][GammaAtom[__]][_] | GammaAtom[__][_] | dsR[__] | dsCoupling[_]) :> memo[x];
   {e1 /. dd -> \[CapitalDelta], pt, memo}];

sampleValue[e_, prop_, w_, deg_, seed_] := Module[{e1, pt, memo, s, cv},
   {e1, pt, memo} = samplePoint[e, prop, w, seed];
   s = LoopMomentumSquared /. pt;
   cv = Thread[cc /@ Range[0, 2] -> (prop["DeltaCoefficients"][s] /. x : (GammaAtom[__][_]) :> memo[x])];
   {integrate[e1, w, deg, cv], validPointQ[prop, s, cv]}];

(* the point is physical: q^2 > 0 and every root of Delta(w) in the upper half plane (Propagator validity conditions) *)
validPointQ[prop_, s_, cv_] := TrueQ[s > 0] && TrueQ[And @@ (prop["ValidityConditions"][s] /. DeltaCoefficient[m_][_] :> (cc[m] /. cv))];

(* a sample at a valid point, for values that are reported as sizes (Neglected); zero tests do not need valid points *)
validSample[e_, prop_, w_, deg_, seed_] := Module[{r},
   Do[r = sampleValue[e, prop, w, deg, seed + 1009 j]; If[r[[2]], Return[r[[1]], Module]], {j, 0, 199}];
   FRGAssert["FrequencyIntegration.validPoint", False, e, "no random point satisfies the validity conditions in 200 tries: check the signs of the ansatz"]];

IntegrateFrequencies[projected_] := Block[{$FRGModule = "FrequencyIntegration"},
  Module[{t0 = AbsoluteTime[], prop, w, deg, res, checks = {}, warns = {}, stage, symbolicC, nI},
   validateStage[projected, "Projection"];
   prop = projected["Loops"]["Propagator"]; w = projected["LoopFrequency"]; deg = prop["DeltaDegree"];
   symbolicC = Thread[cc /@ Range[0, 2] -> (DeltaCoefficient[#][LoopMomentumSquared] & /@ Range[0, 2])];
   res = Function[r, Module[{t1 = AbsoluteTime[], out},
        out = If[r["Class"] === "Declared",
          <|"Value" -> integrate[r["Integrand"], w, deg, symbolicC]|>,
          If[r["Class"] === "Neglected",
           <|"Samples" -> {validSample[r["Integrand"], prop, w, deg, 1]}, "SampleValid" -> True|>,
           <|"Samples" -> Table[First[sampleValue[r["Integrand"], prop, w, deg, sd]], {sd, {1, 2, 3}}]|>]];
        Join[KeyDrop[r, "Integrand"], out, <|"Seconds" -> AbsoluteTime[] - t1|>]]] /@ projected["Integrands"];
   nI = Length[Select[DownValues[FrequencyIntegral], FreeQ[#, Pattern] &]];
   AppendTo[checks, assertCheck["FrequencyIntegration.allIntegrated", True, Null, "",
     ToString[Count[res, KeyValuePattern["Class" -> "Declared"]]] <> " symbolic, " <> ToString[Count[res, _?(KeyExistsQ[#, "Samples"] &)]] <>
      " sampled; no leftover w, Delta or radicals; " <> ToString[nI] <> " distinct residue integrals"]];
   AppendTo[warns, infoCheck["INFO", "FrequencyIntegration.samples", "zero tests use 3 arbitrary random points; Neglected values use one point where q^2 > 0 and the validity conditions hold"]];
   AppendTo[warns, infoCheck["ASSUMPTION", "FrequencyIntegration.poles", "all roots of Delta(w) in the upper half plane (Propagator validity conditions)"]];
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
