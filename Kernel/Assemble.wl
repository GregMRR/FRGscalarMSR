(* ::Package:: *)

(* Assemble.wl -- AssembleEquations: the declared equations in compact form, written in scale-stripped atoms generated from the
   propagator atoms (replaces the hand-written toAtoms/atomVals of DE4plus_fast.wl), then in explicit form with the derived d_s R.
   Objective, assumptions and checks: docs/MODULES.md section "Assemble". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

AssembleEquations::usage = "AssembleEquations[closure] writes every declared equation d_s(comb) = RHS in compact form (scale-stripped \
atoms ScaleAtom[...], abstract dsR) and in explicit form (\"Integrand\": the loop q-integrand in the couplings, dsCoupling[...], \
the shapes r(q^2/k^2) and their derivatives, Dim). Returns the \"Equations\" stage.";
ScaleAtom::usage = "ScaleAtom[a, b, m, j] is a scale-stripped atom of the compact equations: GammaAtom[a, b, m](s) = scale(s) ScaleAtom[a, b, m, 0], \
and its j-th s-derivative = (scale without s) ScaleAtom[a, b, m, j]; the scale is the common factor (number, couplings, power of s) \
of the terms of the atom. Explicit values: the \"Atoms\" key of the Equations stage.";
ShowEquation::usage = "ShowEquation[equations, comb] prints d_s(comb) = RHS in compact form, the atoms it uses and the validity conditions.";
EquationOf::usage = "EquationOf[equations, comb] is the equation record whose flowing combination is comb.";

Begin["`Private`"];

SetAttributes[ScaleAtom, Protected];

(* ---- scale-stripped atoms ---- *)

(* for an atom value v(s): {cscale, n, X} with v = cscale s^n X, cscale = number x couplings common to all terms *)
stripScale[v_, s_, cps_, shapes_] := Module[{val = Expand[v], terms, tm, num, ex},
   terms = If[Head[val] === Plus, List @@ val, {val}];
   tm = terms /. (Alternatives @@ shapes)[__] -> shapeMark /. Derivative[__][Alternatives @@ shapes][__] -> shapeMark;
   num = First[FactorTermsList[val]];
   ex = Table[Min[Exponent[#, x] & /@ tm], {x, Append[cps, s]}];
   {num Times @@ (cps^Most[ex]), Last[ex], Expand[val/(num Times @@ (cps^Most[ex]) s^Last[ex])]}];

(* atom form of GammaAtom[a, b, m] and of its s-derivatives up to jmax: rules at s = LoopMomentumSquared, and atom values *)
atomSystem[prop_, cps_, shapes_, jmax_] := Module[{rules = {}, vals = <||>, forms = <||>},
   Do[Module[{cs, n, x, a, b, m, xj, f},
      {a, b, m} = List @@ g;
      {cs, n, x} = stripScale[prop["AtomValues"][g][xs], xs, cps, shapes];
      Do[
       xj = If[j == 0, x, Expand[D[xs^n x, {xs, j}]]];
       f = If[Head[xj] === Plus, (vals[ScaleAtom[a, b, m, j]] = Function[Evaluate[xs], Evaluate[xj]]; ScaleAtom[a, b, m, j]), xj /. xs -> LoopMomentumSquared];
       f = If[j == 0, cs LoopMomentumSquared^n f, cs f];
       forms[{g, j}] = f;
       AppendTo[rules, If[j == 0, g[LoopMomentumSquared] -> f, Derivative[j][g][LoopMomentumSquared] -> f]],
       {j, 0, jmax}]],
    {g, prop["Atoms"]}];
   {rules, vals, forms}];

(* LeafCount with every atom ScaleAtom[..] and dsR[..] counted as one leaf, as the script's xE, dsNh *)
atomLeaves[e_] := LeafCount[e /. (ScaleAtom | dsR)[__] -> atomLeaf];

(* ---- compact form (DE4plus_fast.wl section 5) ---- *)

collectNum[o_][num_] := Module[{ftl = FactorTermsList[num, o]}, (Times @@ Most[ftl]) Collect[Last[ftl], o, Factor]];
compactFrac[o_][r_] := collectNum[o][Numerator[r]]/Factor[Denominator[r]];

(* candidates: per dsR coefficient or as one fraction, numerator content pulled out and nested-Collected over the atoms; small results
   also get a time-limited Simplify per dsR coefficient; the smallest explicit LeafCount wins and is checked exactly *)
compactify[a_, dsList_, orders_, explicit_] := Module[{tg, cands, best},
   If[a === 0, Return[{0, "zero"}]];
   tg = Together[a];
   cands = Join[
     Table[{"perDs" <> ToString[i], Collect[a, dsList, Collect[Numerator[#], orders[[i]], Factor]/Factor[Denominator[#]] &]}, {i, Length[orders]}],
     Table[{"whole" <> ToString[i], compactFrac[orders[[i]]][tg]}, {i, Length[orders]}]];
   best = First[SortBy[cands, LeafCount[explicit[#[[2]]]] &]];
   If[LeafCount[explicit[best[[2]]]] < 15000,
    AppendTo[cands, {"perDsSimplify", Collect[a, dsList, TimeConstrained[Simplify[#], 10, #] &]}];
    best = First[SortBy[cands, LeafCount[explicit[#[[2]]]] &]]];
   FRGAssert["Assemble.compactExact", ExactlyZeroQ[best[[2]] - a], best[[1]], "internal: compactify changed the expression"];
   {best[[2]], best[[1]]}];

AssembleEquations[closure_] := Block[{$FRGModule = "Assemble"},
  Module[{t0 = AbsoluteTime[], prop, reg, k, qU, cps, shapes, raw, deltaRules, jmax, rules, vals, forms, atomVals, toQ, d0, d1, orders,
    eqs, checks = {}, warns = {}, stage, explicitAtoms, validity, dsrRules, bad, private},
   validateStage[closure, "Closure"];
   prop = closure["Integrated"]["Projection"]["Loops"]["Propagator"];
   reg = closure["Integrated"]["Projection"]["Loops"]["Regulator"];
   k = reg["Scale"]; cps = reg["Couplings"]; shapes = reg["Shapes"];
   qU = Replace[reg["LoopMomentum"], None -> Symbol["Global`q"]];
   toQ = LoopMomentumSquared -> qU^2;

   (* Delta coefficients and their s-derivatives in GammaAtom atoms *)
   deltaRules = {DeltaCoefficient[m_][x_] :> prop["DeltaCoefficients"][x][[m + 1]],
     Derivative[j_][DeltaCoefficient[m_]][x_] :> (D[prop["DeltaCoefficients"][xs][[m + 1]], {xs, j}] /. xs -> x)};
   raw = (#RHS /. deltaRules) & /@ closure["Equations"];
   jmax = Max[0, Cases[raw, Derivative[j_][_GammaAtom][_] :> j, Infinity]];

   (* scale-stripped atoms and the sync check: atom form = explicit value, for every atom and derivative *)
   {rules, vals, forms} = atomSystem[prop, cps, shapes, jmax];
   atomVals = KeyValueMap[#1 -> #2[qU^2] &, vals];
   bad = Select[Keys[forms], ! ExactlyZeroQ[(forms[#] /. toQ /. atomVals) - (D[prop["AtomValues"][#[[1]]][xs], {xs, #[[2]]}] /. xs -> qU^2)] &];
   AppendTo[checks, assertCheck["Assemble.atomSync", bad === {}, bad, "internal: the scale-stripped atoms do not reproduce the propagator atoms",
     ToString[Length[vals]] <> " scale-stripped atoms; atom form = explicit value for " <> ToString[Length[forms]] <> " atoms and derivatives"]];

   (* collection orders: derivative atoms first (they pair with a dsR), then the others in two orders *)
   d1 = Select[Keys[vals], Last[#] > 0 &]; d0 = Select[Keys[vals], Last[#] == 0 &];
   orders = {Join[d1, Reverse[d0]], Join[d1, d0]};
   dsrRules = DsRRules[reg, qU^2];
   AppendTo[checks, assertCheck["Assemble.dsRLinear", AllTrue[dsrRules[[All, 2]], PolynomialQ[#, Cases[#, _dsCoupling, {0, Infinity}]] &&
        Max[0, Exponent[#, Cases[#, _dsCoupling, {0, Infinity}]]] <= 1 &], dsrRules, "internal: d_s R is not linear in the dsCoupling symbols",
     "d_s R linear in dsCoupling[c]"]];

   eqs = MapThread[Function[{r, rw}, Module[{t1 = AbsoluteTime[], a, dsList, cpt, how, integrand, lin, ref},
        a = rw /. rules /. toQ;
        FRGAssert["Assemble.atomsOnly", FreeQ[a, GammaAtom | DeltaCoefficient | \[CapitalDelta] | LoopMomentumSquared], a, "internal: leftover propagator atoms"];
        dsList = Union[Cases[a, _dsR, {0, Infinity}]];
        a = Collect[a, dsList, Together];
        If[a =!= 0 && ExactlyZeroQ[a], a = 0];
        lin = (a /. d_dsR :> 2 d) - 2 a;
        FRGAssert["Assemble.linearInDsR", ExactlyZeroQ[{lin, a /. _dsR -> 0}, "Seeds" -> {1}], r["Flowing"], "the RHS must be linear and homogeneous in the d_s R entries"];
        {cpt, how} = compactify[a, dsList, orders, (# /. atomVals) &];
        integrand = cpt /. atomVals /. dsrRules;
        private = Union[Cases[integrand, x_Symbol /; StringContainsQ[Context[x], "Private`"], {0, Infinity}, Heads -> True]];
        FRGAssert["Assemble.noIntermediates", private === {} && FreeQ[integrand, GammaAtom | ScaleAtom | DeltaCoefficient | \[CapitalDelta] | dsR | LoopMomentumSquared],
         {r["Flowing"], private}, "internal: intermediate symbols left in the explicit equation"];
        Join[KeyTake[r, {"Pair", "Frequency", "Momentum", "GammaK", "Normalisation", "Flowing"}],
         <|"Compact" -> cpt, "Integrand" -> integrand, "DsR" -> dsList, "DsCouplings" -> Union[Cases[integrand, _dsCoupling, {0, Infinity}]],
          "IdenticallyZero" -> (cpt === 0), "CompactForm" -> how,
          "LeafCount" -> <|"Compact" -> atomLeaves[cpt], "Explicit" -> LeafCount[cpt /. atomVals], "Integrand" -> LeafCount[integrand]|>,
          "Seconds" -> AbsoluteTime[] - t1|>]]], {closure["Equations"], raw}];
   AppendTo[checks, assertCheck["Assemble.equations", True, Null, "",
     ToString[Length[eqs]] <> " equations, compact = explicit (exact), no intermediate symbols, linear in dsR and dsCoupling; LeafCounts compact/explicit: " <>
      StringRiffle[(inputForm[#Flowing] <> " " <> ToString[#LeafCount["Compact"]] <> "/" <> ToString[#LeafCount["Integrand"]]) & /@ eqs, ", "]]];
   If[# =!= {}, AppendTo[warns, infoCheck["INFO", "Assemble.zero", "identically zero: " <> StringRiffle[inputForm /@ #, ", "]]]] &[Select[eqs, #IdenticallyZero &][[All, "Flowing"]]];

   validity = prop["ValidityConditions"][LoopMomentumSquared] /. deltaRules /. rules /. toQ;
   explicitAtoms = validity /. atomVals;
   AppendTo[warns, infoCheck["ASSUMPTION", "Assemble.validity", displayString[validity] <> " at every q (causal poles; to be checked at run time)"]];
   AppendTo[warns, infoCheck["INFO", "Assemble.contract", "Integrand: d^dq/(2 Pi)^d integrand of d_s(comb) in the couplings, dsCoupling[c] = k dc/dk (linear), " <>
      "shapes " <> inputForm[shapes] <> " of y = " <> ToString[qU] <> "^2/" <> ToString[k] <> "^2 with derivatives, " <> ToString[closure["Symmetries"]["Dimension"]] <> " symbolic"]];

   stage = stampStage[<|"Stage" -> "Equations", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {closure["Fingerprint"]},
      "Closure" -> closure, "Equations" -> eqs, "Atoms" -> Association[atomVals], "AtomForms" -> Association[KeyValueMap[#1 -> (#2 /. toQ) &, forms]],
      "ValidityConditions" -> <|"Atoms" -> validity, "Explicit" -> explicitAtoms|>, "DsRRules" -> dsrRules,
      "LoopMomentum" -> qU, "Scale" -> k, "Shapes" -> shapes, "Couplings" -> cps, "Dimension" -> closure["Symmetries"]["Dimension"],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

EquationOf[e_Association, comb_] := SelectFirst[e["Equations"], Expand[#Flowing - comb] === 0 &];

ShowEquation[e_Association, comb_] := Module[{r = EquationOf[e, comb], used},
   If[MissingQ[r], Print["no equation for ", comb, "; equations: ", e["Equations"][[All, "Flowing"]]]; Return[Null]];
   Print["d_s(", displayString[comb], ") = Int d^dq/(2 Pi)^d of (", r["CompactForm"], ", LeafCount ", r["LeafCount"]["Compact"], ", explicit ",
    r["LeafCount"]["Integrand"], "):"];
   Print["  ", displayString[r["Compact"]]];
   used = Union[Cases[r["Compact"], _ScaleAtom, {0, Infinity}]];
   If[used =!= {}, Print["atoms:"]; Scan[Print["  ", displayString[#], " = ", displayString[e["Atoms"][#]]] &, used]];
   Print["d_s R entries:"]; Scan[Print["  ", displayString[#[[1]]], " = ", displayString[#[[2]]]] &, Select[e["DsRRules"], MemberQ[r["DsR"], #[[1]]] &]];
   Print["valid where: ", displayString[e["ValidityConditions"]["Atoms"]]]];

End[];
EndPackage[];
