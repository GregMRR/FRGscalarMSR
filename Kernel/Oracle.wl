(* ::Package:: *)

(* Oracle.wl -- RunOracle: a numerical cross-check of the assembled equations (DE4plus_fast.wl section 6c, generalised).
   At one random point it builds Gamma2 + R numerically, inverts it as a matrix, integrates the loop trace over w with NIntegrate
   and takes the frequency and momentum projections by finite differences in nu and in an explicit d-vector p, in d = 2 and 3.
   It is independent of Propagator, Projection, FrequencyIntegration and Assemble; it shares the vertices and the loop formula,
   so it does not test the loop normalisation (E3/E4).
   Objective, assumptions and checks: docs/MODULES.md section "Oracle". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

RunOracle::usage = "RunOracle[equations, opts] compares every equation (and, with \"ZeroFlows\" -> True, every flow that must vanish) \
with a numerical evaluation of the loop at one random point: numeric matrix inverse, NIntegrate over w, finite differences in nu \
and p. Options: \"Dimensions\" -> {2, 3}, \"Tolerance\" -> 10^-8 (relative; absolute for zero flows), \"Seed\" -> 1, \
\"ZeroFlows\" -> False. Returns the \"Oracle\" stage.";

Begin["`Private`"];

(* concrete shapes for the numerical point, as in the script's randomPoint *)
oracleShapes = {Function[y, 1/(1 + y^2)], Function[y, Exp[-y]], Function[y, 2/(2 + y)], Function[y, 1/(1 + y)], Function[y, 3/(3 + y^2)]};

oraclePoint[eqs_, seed_, d_] := Module[{cps = eqs["Couplings"], dsc, pt},
   dsc = dsCoupling /@ cps;
   SeedRandom[seed];
   pt = Join[Thread[cps -> RandomInteger[{1, 30}, Length[cps]]/17], Thread[dsc -> RandomInteger[{-30, 30}, Length[dsc]]/17],
     {eqs["Scale"] -> RandomInteger[{10, 30}]/17, eqs["Dimension"] -> d},
     MapIndexed[#1 -> oracleShapes[[Mod[#2[[1]] - 1, Length[oracleShapes]] + 1]] &, eqs["Shapes"]]];
   pt];

sci[x_, n_] := If[x == 0, "0", ToString[CForm[N[Chop[x, 10^-30], n]]]];

(* finite-difference stencils {shift, weight} of d^m/dx^m at 0 *)
stencil[0, h_] := {{0, 1}};
stencil[1, h_] := {{h, 1/(2 h)}, {-h, -1/(2 h)}};
stencil[2, h_] := {{h, 1/h^2}, {0, -2/h^2}, {-h, 1/h^2}};

Options[RunOracle] = {"Dimensions" -> {2, 3}, "Tolerance" -> 10^-8, "Seed" -> 1, "ZeroFlows" -> False};

RunOracle[eqs_, OptionsPattern[]] := Block[{$FRGModule = "Oracle"},
  Module[{t0 = AbsoluteTime[], closure, loops, v, reg, fields, hasG4, recs, checks = {}, warns = {}, results = {}, stage, tol = OptionValue["Tolerance"],
    h = 10^-5, prec = 40},
   validateStage[eqs, "Equations"];
   closure = eqs["Closure"]; loops = closure["Integrated"]["Projection"]["Loops"];
   v = loops["Vertices"]; reg = loops["Regulator"]; fields = v["Fields"]; hasG4 = loops["HasGamma4"];
   recs = Join[Function[e, <|"Pair" -> e["Pair"], "Frequency" -> e["Frequency"], "Momentum" -> e["Momentum"], "Label" -> inputForm[e["Flowing"]],
         "Symbolic" -> e["Integrand"], "Normalisation" -> e["Normalisation"]|>] /@ eqs["Equations"],
     If[TrueQ[OptionValue["ZeroFlows"]],
      Function[r, <|"Pair" -> r["Pair"], "Frequency" -> r["Frequency"], "Momentum" -> r["Momentum"],
          "Label" -> r["Class"] <> " " <> ToString[r["Pair"]] <> " (I w)^" <> ToString[r["Frequency"]] <> " (q^2)^" <> ToString[r["Momentum"]], "Symbolic" -> 0, "Normalisation" -> 1|>] /@
       Select[closure["Integrated"]["Flows"], MemberQ[{"Protected", "DeclaredNumeric", "Causality"}, #Class] &], {}]];

   Do[Module[{pt, sd = OptionValue["Seed"], qVec, q2, valid, mat, gNum, g3, g4, dsr, fLoop, memo, vecs, dot},
      qVec = Take[{3/10, -1/2, 2/5, 1/7}, d]; q2 = qVec . qVec;
      (* a point where the poles are causal: the residue formula of the symbolic side assumes it *)
      Do[pt = oraclePoint[eqs, sd + 7919 j, d];
       valid = TrueQ[And @@ (eqs["ValidityConditions"]["Explicit"] /. pt /. eqs["LoopMomentum"] -> Sqrt[q2])];
       If[valid, Break[]], {j, 0, 99}];
      FRGAssert["Oracle.validPoint", valid, eqs["ValidityConditions"]["Explicit"], "no random point satisfies the validity conditions"];

      (* numeric Gamma2 + R as a polynomial in w, at momentum squared x2 *)
      mat = N[(Gamma2Matrix[v, wO, xO] /. scal[xO, xO] -> x2O) + RegulatorMatrix[reg, x2O] /. pt, prec + 10];
      gNum[w_, x2_] := Inverse[mat /. {wO -> w, x2O -> x2}];
      vecs = {qO -> qVec, pO -> pVecO};
      dot[e_] := e /. scal[a_, b_] :> ((a /. vecs) . (b /. vecs));
      g3[i_, {nu_, p_}, {w_, q_}] := g3[i, {nu, p}, {w, q}] = N[dot[Gamma3Matrix[v, i, {nu, p}, {w, q}]] /. pt, prec + 10];
      g4[a_, b_] := g4[a, b] = N[dot[Table[VertexValue[v, {a, b, c, e}, {{nuO, pO}, {-nuO, -pO}, {wO, qO}}], {c, fields}, {e, fields}]] /. pt, prec + 10];
      dsr = N[DsRegulatorMatrix[reg, q2] /. pt, prec + 10];
      (* the loop at external (nu, pVec), integrated over w: the same trace formula as Loops *)
      fLoop[{a_, b_}, nu_, pVec_] := memo[{a, b}, nu, pVec] = Module[{g0, gm, tr},
         g0 = gNum[wO, q2]; gm = gNum[wO + nu, (pVec + qVec) . (pVec + qVec)];
         tr = Tr[dsr . g0 . (g3[a, {nuO, pO}, {wO, qO}] /. {nuO -> nu, pVecO -> pVec}) . gm .
             (g3[b, {-nuO, -pO}, {wO + nuO, pO + qO}] /. {nuO -> nu, pVecO -> pVec}) . g0];
         If[hasG4, tr = tr - Tr[dsr . g0 . (g4[a, b] /. {nuO -> nu, pVecO -> pVec}) . g0]/2];
         (* a rational function of w with numeric coefficients; AccuracyGoal stops the search on integrands that vanish *)
         tr = Together[tr/(2 Pi)];
         If[Chop[tr, 10^-(prec - 5)] === 0, 0,
         Quiet[NIntegrate[tr, {wO, -Infinity, -10, 0, 10, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 22, AccuracyGoal -> prec - 8, MaxRecursion -> 30],
          {NIntegrate::izero, NIntegrate::precw}]]];
      Do[Module[{m = r["Frequency"], n = r["Momentum"], units = IdentityMatrix[d], num, sym, err, ok, t1 = AbsoluteTime[]},
         (* (d/d(I nu))^m / m!  and  Laplacian_p / (2 d), by central differences; equations are d_s(comb) = coefficient / norm *)
         num = Total[Flatten[Table[sn[[2]] If[n == 0, fLoop[r["Pair"], sn[[1]], ConstantArray[0, d]],
               Total[Table[sp[[2]] fLoop[r["Pair"], sn[[1]], sp[[1]] units[[c]]], {c, d}, {sp, stencil[2, h]}], 2]],
             {sn, stencil[m, h]}]]]/(I^m m! (2 d)^n r["Normalisation"]);
         sym = N[r["Symbolic"] /. pt /. eqs["LoopMomentum"] -> Sqrt[q2], prec];
         err = If[sym == 0, Abs[num], Abs[(num - sym)/sym]];
         ok = err < tol;
         AppendTo[results, <|"Dimension" -> d, "Label" -> r["Label"], "Numeric" -> num, "Symbolic" -> sym, "Error" -> err, "Seconds" -> AbsoluteTime[] - t1|>];
         AppendTo[checks, assertCheck["Oracle." <> r["Label"] <> " (d = " <> ToString[d] <> ")", ok, {num, sym},
           "the symbolic equation and the numerical loop disagree at the oracle point",
           "numeric " <> sci[num, 10] <> ", symbolic " <> sci[sym, 10] <> If[sym == 0, ", |numeric| = ", ", rel. diff "] <> sci[err, 2]]]],
       {r, recs}]],
    {d, OptionValue["Dimensions"]}];
   AppendTo[warns, infoCheck["ASSUMPTION", "Oracle.scope",
     "the oracle shares the vertices and the loop formula (prefactor, tadpole -1/2, 1/(2 Pi)) with the symbolic chain: it does not test the normalisation (E3/E4)"]];
   AppendTo[warns, infoCheck["INFO", "Oracle.method", "numeric inverse of Gamma2 + R, NIntegrate over w (" <> ToString[prec] <>
      " digits), central differences h = 10^-5 in nu and in an explicit p vector, d in " <> ToString[OptionValue["Dimensions"]]]];
   stage = stampStage[<|"Stage" -> "Oracle", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {eqs["Fingerprint"]}, "Results" -> results,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

End[];
EndPackage[];
