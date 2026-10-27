(* ::Package:: *)

(* ExternalMomentum.wl -- BuildLOEquations: the LO pipeline (design/LOPIPELINE.md step 5). The Gamma2 couplings are a reference times
   a function F(p^2); the flow of each function is the one-loop d_s Gamma2_ab(nu = 0, p) at generic external momentum p, with the
   loop frequency integrated by residues over the two pole families Delta(w, q^2) and Delta(w, Q^2), Q = p + q. It replaces
   ProjectOnMonomials, IntegrateFrequencies, CheckClosure and AssembleEquations of the derivative expansion and returns an
   "Equations" stage with "Pipeline" -> "LO".
   Objective, assumptions and checks: docs/MODULES.md section "ExternalMomentum". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

BuildLOEquations::usage = "BuildLOEquations[loops, symmetries] builds the LO flow equations: for each function F of the ansatz \
(Gamma2_ab(0, p) = norm ref p^(2n) F(p^2)) the integrand of d_s(ref F(p^2)) = Integral d^dq/(2 Pi)^d Integrand(p, q, CosTheta), \
already divided by norm p^(2n), and the p -> 0 limit (integrand of d_s ref, which fixes eta through F(0) = 1). Checks the zero \
entries at every p, the Protected i w entries at p = 0, the low orders below p^(2n), the angles and the shape derivatives. \
Returns the \"Equations\" stage with \"Pipeline\" -> \"LO\".";
CosTheta::usage = "CosTheta is the cosine of the angle between the external momentum p and the loop momentum q in LO integrands: \
p.q = p q CosTheta, Q^2 = p^2 + q^2 + 2 p q CosTheta.";
ShowLOEquation::usage = "ShowLOEquation[equations, f] prints the LO flow of the function f: prefactor, divided power of p^2, loop \
integrand and its p -> 0 limit.";

Begin["`Private`"];

SetAttributes[CosTheta, Protected];

(* explicit values of the propagator atoms (and their s-derivatives), the Delta coefficients and d_s R at the loop momentum *)
loExplicit[e_, prop_, reg_, sQ_] := Module[{av = prop["AtomValues"], dc = prop["DeltaCoefficients"]},
   (e /. {DeltaCoefficient[m_][s_] :> dc[s][[m + 1]], Derivative[j_][DeltaCoefficient[m_]][s_] :> (D[dc[loXs][[m + 1]], {loXs, j}] /. loXs -> s)}
      /. {(h : GammaAtom[__])[s_] :> av[h][s], Derivative[j_][h : GammaAtom[__]][s_] :> (D[av[h][loXs], {loXs, j}] /. loXs -> s)})
     /. DsRRules[reg, sQ]];

(* vectors -> magnitudes: p.p -> p^2, q.q -> q^2, p.q -> p q c *)
loMagnitudes[e_, pV_, qV_, pU_, qU_] := e /. {scal[pV, pV] -> pU^2, scal[qV, qV] -> qU^2, scal[pV, qV] -> pU qU CosTheta};

(* the (I nu)^m (p^2)^n projection at nu = p = 0 (Projection.wl), integrated over w (single momentum), atoms explicit, q^2 -> qU^2 *)
loProjectedAtZero[loop_, m_, n_, {nu_, pV_}, wq_, dim_, prop_, reg_, qU_] := Module[{e, deg = prop["DeltaDegree"]},
   e = project[loop, m, n, {nu, pV}, wq, dim, deg];
   e = integrate[e, First[wq], deg, Thread[FrequencyIntegralC /@ Range[0, 2] -> (DeltaCoefficient[#][LoopMomentumSquared] & /@ Range[0, 2])]];
   loExplicit[e, prop, reg, LoopMomentumSquared] /. LoopMomentumSquared -> qU^2];

(* function atoms of the LO integrands: F[x], F'[x], ... and shapes r[x], r'[x], ... *)
loDerivativeOrders[e_, heads_] := Union[Cases[e, Derivative[k_][h : Alternatives @@ heads][_] :> {h, k}, {0, Infinity}]];

BuildLOEquations[loops_, symmetries_] := Block[{$FRGModule = "ExternalMomentum"},
  Module[{t0 = AbsoluteTime[], v, prop, reg, ansatz, funcs, fields, dim, deg, wq, wL, qL, pV, nu, qU, pU, kk, shapes, fs, pairs, flowRec, prot,
    negl, checks = {}, warns = {}, zeroBad = {}, neglected = {}, wardBad = {}, wardDone = {}, eqs, nFam, bad, stage, sq, sQ, atZero, loopOf,
    integratedOf, cps, validity},
   validateStage[loops, "Loops"]; validateStage[symmetries, "Symmetries"];
   FRGAssert["ExternalMomentum.pipeline", loops["Pipeline"] === "LO" && symmetries["Pipeline"] === "LO", {loops["Pipeline"], symmetries["Pipeline"]},
    "BuildLOEquations needs an LO ansatz (fn[...]); a derivative-expansion ansatz goes through ProjectOnMonomials"];
   v = loops["Vertices"]; prop = loops["Propagator"]; reg = loops["Regulator"]; ansatz = v["Ansatz"];
   FRGAssert["ExternalMomentum.sameAnsatz", v["Upstream"] === symmetries["Upstream"], {v["Upstream"], symmetries["Upstream"]},
    "the loops and the symmetries were built from different ansatz stages"];
   funcs = ansatz["Functions"]; fs = funcs[[All, 1]]; fields = v["Fields"]; cps = ansatz["Couplings"];
   dim = symmetries["Dimension"]; deg = prop["DeltaDegree"]; wq = LoopVariables[loops]; {wL, qL} = wq;
   kk = reg["Scale"]; shapes = reg["Shapes"];
   qU = Replace[reg["LoopMomentum"], None -> Symbol["Global`q"]]; pU = Symbol["Global`p"];
   pV = loPV; nu = loNu;
   sq = scal[qL, qL]; sQ = Expand[scal[pV + qL, pV + qL]];
   prot = Select[symmetries["Tags"], #Kind === "Protected" &]; negl = Select[symmetries["Tags"], #Kind === "Neglected" &];
   pairs = Select[Tuples[fields, 2], OrderedQ[{First[FirstPosition[fields, #[[1]]]], First[FirstPosition[fields, #[[2]]]]}] &];

   (* the flowing entries: Gamma2_ab(0, p) = norm ref (p^2)^n F(p^2) *)
   flowRec = Association[Function[pr, Module[{g = Expand[VertexValue[v, pr, {{0, pV}}]], f, pre, n, ref, norm},
        f = FirstCase[g, (h : Alternatives @@ fs)[_] :> h, None, {0, Infinity}];
        If[f === None, Nothing,
         pre = Together[g/f[scal[pV, pV]]]; ref = f /. funcs;
         n = Exponent[pre, scal[pV, pV]]; norm = Together[pre/(ref scal[pV, pV]^n)];
         FRGAssert["ExternalMomentum.prefactor", FreeQ[pre, Alternatives @@ fs] && NumericQ[norm], {pr, g},
          "internal: Gamma2_ab(0, p) is not number x reference x (p^2)^n x F(p^2)"];
         pr -> <|"Pair" -> pr, "Function" -> f, "Reference" -> ref, "Normalisation" -> norm, "P2Divided" -> n,
           "Prefactor" -> norm ref pU^(2 n), "GammaK" -> norm ref pU^(2 n) f[pU^2]|>]]] /@ pairs];
   FRGAssert["ExternalMomentum.functionsFlow", Sort[Values[flowRec][[All, "Function"]]] === Sort[fs], Values[flowRec][[All, "Function"]],
    "internal: every function must sit in exactly one Gamma2 entry"];

   (* the loop at nu = 0, generic p, integrated over w by residues (two pole families), atoms explicit, magnitudes *)
   loopOf[pr_] := loopOf[pr] = LoopIntegrand[loops, pr, nu, pV];
   integratedOf[pr_] := integratedOf[pr] = Module[{e = loopOf[pr] /. nu -> 0, x},
      x = loIntegrate[e, wL, {sq, sQ}, deg];
      FRGAssert["ExternalMomentum.isotropicAngles", FreeQ[loMagnitudes[x, pV, qL, pU, qU], scal | wedge | pV | qL], pr,
       "the loop depends on the angle other than through p.q (wedge / jac terms): not in LO (FixedDimension LO is in docs/FUTURE.md)"];
      Together[loMagnitudes[loExplicit[x, prop, reg, sq], pV, qL, pU, qU]]];

   AppendTo[checks, assertCheck["ExternalMomentum.isotropicAngles", AllTrue[pairs, (integratedOf[#]; True) &], Null, "",
     "every entry depends on the angle only through c = cos(p, q)"]];
   nFam = Length[Select[DownValues[FrequencyIntegralFamilies], FreeQ[#, Pattern] &]];
   AppendTo[checks, assertCheck["FrequencyIntegration.regularAtCoincidence", True, Null, "",
     ToString[nFam] <> " two-family residue integrals (Delta at q and at Q = p + q, nu = 0), each regular when the Q data equal the q data"]];

   (* entries outside the ansatz: zero at every p (causality, structure), or dropped by a Neglected tag *)
   Scan[Function[pr, If[! KeyExistsQ[flowRec, pr],
       Which[
        AnyTrue[negl, pairMatchQ[#Pattern, pr] &], AppendTo[neglected, pr],
        ! ExactlyZeroQ[integratedOf[pr]], AppendTo[zeroBad, pr]]]], pairs];
   AppendTo[checks, assertCheck["ExternalMomentum.zeroEntries", zeroBad === {}, zeroBad,
     "these Gamma2 entries are not in the ansatz (or vanish by causality) but their flow at nu = 0 is not zero: add a function, or a Neglected[pair] tag",
     "entries outside the ansatz vanish exactly at every p: " <> StringRiffle[ToString /@ Complement[pairs, Keys[flowRec], neglected], ", "]]];
   If[neglected =!= {}, AppendTo[warns, infoCheck["WARN", "ExternalMomentum.neglected", "Neglected entries (flow computed, dropped): " <>
       StringRiffle[(ToString[#] <> " = " <> inputForm[integratedOf[#]]) & /@ neglected, "; "]]]];

   (* Protected i w coefficients: exactly at p = 0 (d/d(I nu), then the p^(2n) average, then nu -> 0); neglected at p != 0 (L7) *)
   Scan[Function[t, Scan[Function[pr, If[t["Frequency"] === 1 && pairMatchQ[t["Pattern"], pr],
        AppendTo[wardDone, {pr, t["Momentum"]}];
        If[! ExactlyZeroQ[loProjectedAtZero[loopOf[pr], 1, t["Momentum"], {nu, pV}, wq, dim, prop, reg, qU]], AppendTo[wardBad, {pr, t["Reason"]}]]]], pairs]], prot];
   AppendTo[checks, assertCheck["ExternalMomentum.ward", wardBad === {}, wardBad,
     "a Protected i w coefficient flows at p = 0: the declared Ward identity does not hold for this ansatz",
     If[wardDone === {}, "no Protected i w coefficient", "Protected i w coefficients flow to exactly 0 at p = 0: " <> inputForm[wardDone]]]];
   If[wardDone =!= {}, AppendTo[warns, infoCheck["INFO", "ExternalMomentum.wardNeglectedLO",
      "Neglected (LO): the Protected i w coefficients at p != 0 (LO drops the frequency dependence of the functions)"]]];

   (* the flow of each function and its p -> 0 limit (eta) *)
   eqs = Function[r, Module[{pr = r["Pair"], n = r["P2Divided"], norm = r["Normalisation"], flow, eta},
       flow = Together[integratedOf[pr]/(norm pU^(2 n))];
       eta = Together[loProjectedAtZero[loopOf[pr], 0, n, {nu, pV}, wq, dim, prop, reg, qU]/norm];
       Join[r, <|"Flowing" -> r["Reference"] r["Function"][pU^2], "Integrand" -> flow, "EtaIntegrand" -> eta,
         "DsCouplings" -> Union[Cases[flow, _dsCoupling, {0, Infinity}]], "IdenticallyZero" -> ExactlyZeroQ[flow]|>]]] /@ Values[flowRec];

   (* L5: the eta rows contain only F, F', r, r' (BMWflow has no second derivatives; no integration by parts) *)
   bad = Select[eqs, Max[0, loDerivativeOrders[#EtaIntegrand, Join[fs, shapes]][[All, 2]]] > 1 &];
   AppendTo[checks, assertCheck["ExternalMomentum.firstDerivativesOnly", bad === {},
     Thread[bad[[All, "Function"]] -> (loDerivativeOrders[#EtaIntegrand, Join[fs, shapes]] & /@ bad)],
     "the p -> 0 limit (eta row) needs a second derivative of a function or of a shape: outside the LO scope (KPZ-like models; docs/FUTURE.md)",
     "eta rows contain only F, F', r, r'"]];

   (* every order below p^(2n) of a divided entry vanishes (this also checks the Protected w^0 low-p entries) *)
   bad = {};
   Scan[Function[r, Do[If[! ExactlyZeroQ[loProjectedAtZero[loopOf[r["Pair"]], 0, j, {nu, pV}, wq, dim, prop, reg, qU]], AppendTo[bad, {r["Pair"], "p^" <> ToString[2 j]}]],
      {j, 0, r["P2Divided"] - 1}]], eqs];
   Scan[Function[t, Scan[Function[r, If[t["Frequency"] === 0 && pairMatchQ[t["Pattern"], r["Pair"]] && t["Momentum"] >= r["P2Divided"],
        AppendTo[bad, {r["Pair"], "Protected p^" <> ToString[2 t["Momentum"]] <> " is not below the divided p^" <> ToString[2 r["P2Divided"]]}]]], eqs]], prot];
   AppendTo[checks, assertCheck["ExternalMomentum.dividedOrders", bad === {}, bad,
     "a row divided by p^(2n) has a nonzero lower order (it would be exported as 1/p^2, ...): the ansatz breaks a symmetry it assumes, \
or a Protected order is not below the function's prefactor",
     "orders below p^(2n) vanish exactly: " <> StringRiffle[(ToString[#Function] <> " (n = " <> ToString[#P2Divided] <> ")") & /@ eqs, ", "]]];

   validity = Join[prop["ValidityConditions"][qU^2], prop["ValidityConditions"][pU^2 + qU^2 + 2 pU qU CosTheta]] /.
     {DeltaCoefficient[m_][s_] :> prop["DeltaCoefficients"][s][[m + 1]]} /. (h : GammaAtom[__])[s_] :> prop["AtomValues"][h][s];
   AppendTo[warns, infoCheck["ASSUMPTION", "ExternalMomentum.causalPoles", "every root of Delta(w) at q and at Q in the upper half plane " <>
      "(residue assumption, not exported in LO: docs/FUTURE.md): " <> displayString[validity]]];
   AppendTo[warns, infoCheck["INFO", "ExternalMomentum.contract", "Integrand: d^dq/(2 Pi)^d integrand of d_s(ref F(p^2)) divided by norm p^(2n), in p, " <>
      ToString[qU] <> ", CosTheta, F at q^2 and Q^2, shapes r(q^2/k^2), dsCoupling[c]; EtaIntegrand: its p -> 0 limit (integrand of d_s ref)"]];

   stage = stampStage[<|"Stage" -> "Equations", "Pipeline" -> "LO", "Upstream" -> {loops["Fingerprint"], symmetries["Fingerprint"]},
      "Loops" -> loops, "Symmetries" -> symmetries, "Equations" -> eqs, "DsRRules" -> DsRRules[reg, qU^2],
      "ValidityConditions" -> <|"Explicit" -> validity|>,
      "LoopMomentum" -> qU, "ExternalMomentum" -> pU, "Cosine" -> CosTheta, "Scale" -> kk, "Shapes" -> shapes, "Functions" -> funcs,
      "Couplings" -> cps, "Dimension" -> dim, "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   Scan[ShowLOEquation[stage, #] &, fs];
   PrintCard[stage];
   stage]];

ShowLOEquation[e_Association, f_] := Module[{r = SelectFirst[e["Equations"], #Function === f &]},
   If[MissingQ[r], Print["no function ", f, "; functions: ", e["Equations"][[All, "Function"]]]; Return[Null]];
   Print["d_s(", displayString[r["Flowing"]], ") = Int d^dq/(2 Pi)^d of (Gamma2", r["Pair"], "(0, p) = ", displayString[r["GammaK"]],
    ", divided by ", displayString[r["Prefactor"]/r["Reference"]], "):"];
   Print["  ", displayString[r["Integrand"]]];
   Print["p -> 0 (d_s ", r["Reference"], "): ", displayString[r["EtaIntegrand"]]]];

End[];
EndPackage[];
