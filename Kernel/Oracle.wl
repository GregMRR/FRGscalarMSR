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
and p. Options: \"Dimensions\" -> Automatic ({2, 3} for a symbolic d, {2} for a declared d = 2), \"Tolerance\" -> 10^-8 (relative; absolute for zero flows), \"Seed\" -> 1, \
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

(* (Laplacian_p)^n at p = 0 by nested central second differences over direction tuples (c1, ..., cn), and its normalisation
   N_n = (Laplacian)^n (p^2)^n = 4^n n! (d/2)_n (written here independently of Projection). For n >= 1, Richardson extrapolation
   (h, h/2) removes the O(h^2) error (review E19; at n = 1 it keeps the relative error below the tolerance when the coefficient is small). f is a function of the explicit d-vector p. *)
oracleLapNorm[n_, d_] := 4^n n! Pochhammer[d/2, n];
oracleLapStencil[f_, d_, n_, h_] := With[{units = IdentityMatrix[d], st = stencil[2, h]},
   Total[Flatten[Table[Times @@ sh[[All, 2]] f[ConstantArray[0, d] + Sum[sh[[i, 1]] units[[cs[[i]]]], {i, n}]], {cs, Tuples[Range[d], n]}, {sh, Tuples[st, n]}]]]];
oracleLaplacianN[f_, d_, n_, h_] := If[n == 0, oracleLapStencil[f, d, n, h], (4 oracleLapStencil[f, d, n, h/2] - oracleLapStencil[f, d, n, h])/3];

(* eps_ij in the numeric vertices: wedge[a, b] -> a1 b2 - a2 b1, evaluated once both vectors are explicit (p is numeric only in fLoop) *)
wedgeN[a_?VectorQ, b_?VectorQ] /; Length[a] == Length[b] == 2 := a[[1]] b[[2]] - a[[2]] b[[1]];

(* ---------- LO branch (design/LOPIPELINE.md step 6) ----------
   At a random point (couplings, d_s couplings, k, concrete shapes and concrete positive test functions F) and explicit d-vectors p, q:
   numeric inverse of Gamma2 + R, the loop at nu = 0 integrated over w with NIntegrate, against the symbolic residue result
   (the equation's Integrand times norm p^(2n)). Independent of FrequencyIntegralFamilies and of the atom bookkeeping. *)
oracleFunctions = {Function[x, (4 + 3 x)/(4 + x)], Function[x, (5 + x^2)/(5 + 2 x)], Function[x, 2/(2 + x) + x/(3 + x)]};

oracleLO[eqs_, dims_, tol_, seed_, zeroQ_] := Module[{t0 = AbsoluteTime[], loops, v, reg, fields, hasG4, funcs, qU, pU, recs, results = {}, checks = {}, warns = {},
    stage, prec = 40, pairs, negPats},
   loops = eqs["Loops"]; v = loops["Vertices"]; reg = loops["Regulator"]; fields = v["Fields"]; hasG4 = loops["HasGamma4"];
   funcs = eqs["Functions"][[All, 1]]; qU = eqs["LoopMomentum"]; pU = eqs["ExternalMomentum"];
   pairs = Select[Tuples[fields, 2], OrderedQ[{First[FirstPosition[fields, #[[1]]]], First[FirstPosition[fields, #[[2]]]]}] &];
   negPats = Select[eqs["Symmetries"]["Tags"], #Kind === "Neglected" &][[All, "Pattern"]];
   recs = Join[Function[e, <|"Pair" -> e["Pair"], "Label" -> inputForm[e["Flowing"]], "Symbolic" -> e["Integrand"] e["Normalisation"] pU^(2 e["P2Divided"])|>] /@ eqs["Equations"],
     If[TrueQ[zeroQ], <|"Pair" -> #, "Label" -> "zero entry " <> ToString[#], "Symbolic" -> 0|> & /@ Select[pairs,
        ! MemberQ[eqs["Equations"][[All, "Pair"]], #] && NoneTrue[negPats, Function[pat, pairMatchQ[pat, #]]] &], {}]];
   Do[Module[{pt, qVec, pVec, q2, Q2, c, valid, mat, gNum, g3, g4, dsr, num, sym, err, ok, t1},
      qVec = Take[{3/10, -1/2, 2/5}, d]; pVec = Take[{2/7, 1/3, -1/5}, d];
      q2 = qVec . qVec; Q2 = (pVec + qVec) . (pVec + qVec); c = pVec . qVec/Sqrt[(pVec . pVec) q2];
      Do[pt = Join[oraclePoint[eqs, seed + 7919 j, d], MapIndexed[#1 -> oracleFunctions[[Mod[#2[[1]] - 1, Length[oracleFunctions]] + 1]] &, funcs]];
       valid = TrueQ[And @@ (eqs["ValidityConditions"]["Explicit"] /. pt /. {qU -> Sqrt[q2], pU -> Sqrt[pVec . pVec], CosTheta -> c})];
       If[valid, Break[]], {j, 0, 99}];
      FRGAssert["Oracle.validPoint", valid, eqs["ValidityConditions"]["Explicit"], "no random point satisfies the validity conditions at q and Q"];
      mat = N[(Gamma2Matrix[v, wO, xO] /. scal[xO, xO] -> x2O) + RegulatorMatrix[reg, x2O] /. pt, prec + 10];
      gNum[w_, x2_] := Inverse[mat /. {wO -> w, x2O -> x2}];
      g3[i_, {nu_, p_}, {w_, q_}] := N[Gamma3Matrix[v, i, {nu, p}, {w, q}] /. {scal[a_, b_] :> ((a /. {qO -> qVec, pO -> pVec}) . (b /. {qO -> qVec, pO -> pVec}))} /. pt, prec + 10];
      g4[a_, b_] := N[Table[VertexValue[v, {a, b, cc, e}, {{0, pO}, {0, -pO}, {wO, qO}}], {cc, fields}, {e, fields}] /.
          {scal[x_, y_] :> ((x /. {qO -> qVec, pO -> pVec}) . (y /. {qO -> qVec, pO -> pVec}))} /. pt, prec + 10];
      dsr = N[DsRegulatorMatrix[reg, q2] /. pt, prec + 10];
      Do[t1 = AbsoluteTime[];
       num = Module[{g0 = gNum[wO, q2], gm = gNum[wO, Q2], tr},
         tr = Tr[dsr . g0 . g3[r["Pair"][[1]], {0, pO}, {wO, qO}] . gm . g3[r["Pair"][[2]], {0, -pO}, {wO, pO + qO}] . g0];
         If[hasG4, tr = tr - Tr[dsr . g0 . g4 @@ r["Pair"] . g0]/2];
         tr = Together[tr/(2 Pi)];
         If[Chop[tr, 10^-(prec - 5)] === 0, 0,
          Quiet[NIntegrate[tr, {wO, -Infinity, -10, 0, 10, Infinity}, WorkingPrecision -> prec, PrecisionGoal -> 22, AccuracyGoal -> prec - 8, MaxRecursion -> 30],
           {NIntegrate::izero, NIntegrate::precw}]]];
       sym = N[r["Symbolic"] /. pt /. {qU -> Sqrt[q2], pU -> Sqrt[pVec . pVec], CosTheta -> c}, prec];
       err = If[sym == 0, Abs[num], Abs[(num - sym)/sym]];
       ok = err < tol;
       AppendTo[results, <|"Dimension" -> d, "Label" -> r["Label"], "Numeric" -> num, "Symbolic" -> sym, "Error" -> err, "Seconds" -> AbsoluteTime[] - t1|>];
       AppendTo[checks, assertCheck["Oracle." <> r["Label"] <> " (d = " <> ToString[d] <> ")", ok, {num, sym},
         "the LO residue result and the numerical w integral of the loop disagree at the oracle point",
         "numeric " <> sci[num, 10] <> ", symbolic " <> sci[sym, 10] <> If[sym == 0, ", |numeric| = ", ", rel. diff "] <> sci[err, 2]]], {r, recs}]],
    {d, dims}];
   AppendTo[warns, infoCheck["ASSUMPTION", "Oracle.scope",
     "the oracle shares the vertices and the loop formula with the symbolic chain: it tests the two-family residues, the atoms and the explicit integrand"]];
   AppendTo[warns, infoCheck["INFO", "Oracle.methodLO", "nu = 0, explicit p and q vectors, numeric inverse of Gamma2 + R with concrete shapes and test functions F, " <>
      "NIntegrate over w (" <> ToString[prec] <> " digits), d in " <> ToString[dims]]];
   stage = stampStage[<|"Stage" -> "Oracle", "Pipeline" -> "LO", "Upstream" -> {eqs["Fingerprint"]}, "Results" -> results,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage];

Options[RunOracle] = {"Dimensions" -> Automatic, "Tolerance" -> 10^-8, "Seed" -> 1, "ZeroFlows" -> False};

RunOracle[eqs_, OptionsPattern[]] := Block[{$FRGModule = "Oracle"},
  Module[{t0 = AbsoluteTime[], closure, loops, v, reg, fields, hasG4, recs, checks = {}, warns = {}, results = {}, stage, tol = OptionValue["Tolerance"],
    h = 10^-5, prec = 40, dims},
   validateStage[eqs, "Equations"];
   dims = Replace[OptionValue["Dimensions"], Automatic -> If[eqs["Dimension"] === 2, {2}, {2, 3}]];
   FRGAssert["Oracle.declaredDimension", If[eqs["Dimension"] === 2, dims === {2}, ListQ[dims] && dims =!= {} && SubsetQ[{2, 3}, dims]], dims,
    If[eqs["Dimension"] === 2, "the dimension is declared: d = 2 only (\"Dimensions\" -> {2} or Automatic)",
     "\"Dimensions\" must be a non-empty subset of {2, 3}"]];
   If[eqs["Pipeline"] === "LO", Return[oracleLO[eqs, dims, tol, OptionValue["Seed"], OptionValue["ZeroFlows"]], Module]];
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
      (* numerics per working precision wp: 40 digits (PrecisionGoal 22) for n = 0 as before; more for n >= 1, where 2n-th
         differences at h = 10^-5 eat 10 n digits (reviews E3, E15: the precision is part of every memo key) *)
      mat[wp_] := mat[wp] = N[(Gamma2Matrix[v, wO, xO] /. scal[xO, xO] -> x2O) + RegulatorMatrix[reg, x2O] /. pt, wp + 10];
      gNum[w_, x2_, wp_] := Inverse[mat[wp] /. {wO -> w, x2O -> x2}];
      vecs = {qO -> qVec, pO -> pVecO};
      dot[e_] := e /. {scal[a_, b_] :> ((a /. vecs) . (b /. vecs)), wedge[a_, b_] :> wedgeN[a /. vecs, b /. vecs]};
      g3[i_, {nu_, p_}, {w_, q_}, wp_] := g3[i, {nu, p}, {w, q}, wp] = N[dot[Gamma3Matrix[v, i, {nu, p}, {w, q}]] /. pt, wp + 10];
      g4[a_, b_, wp_] := g4[a, b, wp] = N[dot[Table[VertexValue[v, {a, b, c, e}, {{nuO, pO}, {-nuO, -pO}, {wO, qO}}], {c, fields}, {e, fields}]] /. pt, wp + 10];
      dsr[wp_] := dsr[wp] = N[DsRegulatorMatrix[reg, q2] /. pt, wp + 10];
      (* the loop at external (nu, pVec), integrated over w: the same trace formula as Loops *)
      fLoop[{a_, b_}, nu_, pVec_, {wp_, pg_}] := memo[{a, b}, nu, pVec, {wp, pg}] = Module[{g0, gm, tr},
         g0 = gNum[wO, q2, wp]; gm = gNum[wO + nu, (pVec + qVec) . (pVec + qVec), wp];
         tr = Tr[dsr[wp] . g0 . (g3[a, {nuO, pO}, {wO, qO}, wp] /. {nuO -> nu, pVecO -> pVec}) . gm .
             (g3[b, {-nuO, -pO}, {wO + nuO, pO + qO}, wp] /. {nuO -> nu, pVecO -> pVec}) . g0];
         If[hasG4, tr = tr - Tr[dsr[wp] . g0 . (g4[a, b, wp] /. {nuO -> nu, pVecO -> pVec}) . g0]/2];
         (* a rational function of w with numeric coefficients; AccuracyGoal stops the search on integrands that vanish *)
         tr = Together[tr/(2 Pi)];
         If[Chop[tr, 10^-(wp - 5)] === 0, 0,
         Quiet[NIntegrate[tr, {wO, -Infinity, -10, 0, 10, Infinity}, WorkingPrecision -> wp, PrecisionGoal -> pg, AccuracyGoal -> wp - 8, MaxRecursion -> 30],
          {NIntegrate::izero, NIntegrate::precw}]]];
      Do[Module[{m = r["Frequency"], n = r["Momentum"], units = IdentityMatrix[d], num, sym, err, ok, t1 = AbsoluteTime[], pspec},
         (* (d/d(I nu))^m / m!  and  (Laplacian_p)^n / N_n, by central differences; equations are d_s(comb) = coefficient / norm *)
         If[n == 0,
          pspec = {prec, 22};
          num = Total[Flatten[Table[sn[[2]] If[n == 0, fLoop[r["Pair"], sn[[1]], ConstantArray[0, d], pspec],
                Total[Table[sp[[2]] fLoop[r["Pair"], sn[[1]], sp[[1]] units[[c]], pspec], {c, d}, {sp, stencil[2, h]}], 2]],
              {sn, stencil[m, h]}]]]/(I^m m! (2 d)^n r["Normalisation"]),
          pspec = {8 + 10 n + 6 + 18, 8 + 10 n + 6};
          num = Total[Table[sn[[2]] oracleLaplacianN[fLoop[r["Pair"], sn[[1]], #, pspec] &, d, n, h], {sn, stencil[m, h]}]]/
            (I^m m! oracleLapNorm[n, d] r["Normalisation"])];
         sym = N[r["Symbolic"] /. pt /. eqs["LoopMomentum"] -> Sqrt[q2], First[pspec]];
         err = If[sym == 0, Abs[num], Abs[(num - sym)/sym]];
         ok = err < tol;
         AppendTo[results, <|"Dimension" -> d, "Label" -> r["Label"], "Numeric" -> num, "Symbolic" -> sym, "Error" -> err, "Seconds" -> AbsoluteTime[] - t1|>];
         AppendTo[checks, assertCheck["Oracle." <> r["Label"] <> " (d = " <> ToString[d] <> ")", ok, {num, sym},
           "the symbolic equation and the numerical loop disagree at the oracle point",
           "numeric " <> sci[num, 10] <> ", symbolic " <> sci[sym, 10] <> If[sym == 0, ", |numeric| = ", ", rel. diff "] <> sci[err, 2]]]],
       {r, recs}]],
    {d, dims}];
   AppendTo[warns, infoCheck["ASSUMPTION", "Oracle.scope",
     "the oracle shares the vertices and the loop formula (prefactor, tadpole -1/2, 1/(2 Pi)) with the symbolic chain: it does not test the normalisation (E3/E4)"]];
   AppendTo[warns, infoCheck["INFO", "Oracle.method", "numeric inverse of Gamma2 + R, NIntegrate over w (" <> ToString[prec] <>
      " digits), central differences h = 10^-5 in nu and in an explicit p vector, d in " <> ToString[dims]]];
   If[Max[0, recs[[All, "Momentum"]]] >= 1, AppendTo[warns, infoCheck["INFO", "Oracle.highMomentum",
      "(q^2)^(n >= 1): nested stencils / N_n with Richardson (h, h/2), " <> ToString[8 + 10 # + 6 + 18] <> " digits at n = " <> ToString[#] <>
       "; the oracle runs d = 2, 3 only: the d-dependence at n >= 2 is certified by Projection.explicitVectors"] &[Max[recs[[All, "Momentum"]]]]]];
   stage = stampStage[<|"Stage" -> "Oracle", "Pipeline" -> eqs["Pipeline"], "Upstream" -> {eqs["Fingerprint"]}, "Results" -> results,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

End[];
EndPackage[];
