(* ::Package:: *)

(* Closure.wl -- CheckClosure: decide, monomial by monomial, whether the integrated flows close on the ansatz.
   Protected, DeclaredNumeric and Causality flows must vanish (else abort with the tag's reason); Neglected flows are reported
   with their value at one valid random point; every declared monomial becomes an equation d_s(combination) = RHS.
   Objective, assumptions and checks: docs/MODULES.md section "Closure". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

CheckClosure::usage = "CheckClosure[integrated, symmetries] asserts that every Protected, DeclaredNumeric and Causality flow vanishes, \
warns about every Neglected flow with its value at a valid random point, and turns every declared monomial with Gamma_k coefficient \
norm * comb (norm the pure number) into the equation d_s(comb) = Value / norm. Returns the \"Closure\" stage.";
ShowClosureReport::usage = "ShowClosureReport[closure] prints, for every monomial, its class and the outcome of the closure check.";

Begin["`Private`"];

monomialString[r_] := "(I w)^" <> ToString[r["Frequency"]] <> " (q^2)^" <> ToString[r["Momentum"]];
numString[x_] := ToString[CForm[N[x]]];
recordString[r_] := ToString[r["Pair"]] <> " " <> monomialString[r];

(* numeric content (sign included) and flowing combination: -2 Dhk aDhk -> {-2, Dhk aDhk} *)
splitCoefficient[c_] := Module[{f = FactorTermsList[c]}, {Times @@ Most[f], Last[f]}];

CheckClosure[integrated_, symmetries_] := Block[{$FRGModule = "Closure"},
  Module[{t0 = AbsoluteTime[], flows, checks = {}, warns = {}, zeroCl, bad, eqs, neg, report, stage},
   FRGAssert["Closure.pipeline", NoneTrue[{integrated, symmetries}, AssociationQ[#] && Lookup[#, "Pipeline", None] === "LO" &], "LO", "this is a derivative-expansion stage; an LO ansatz (fn[...]) goes through BuildLOEquations (design/LOPIPELINE.md)"];
   validateStage[integrated, "Integrated"]; validateStage[symmetries, "Symmetries"];
   FRGAssert["Closure.sameSymmetries", MemberQ[integrated["Projection"]["Upstream"], symmetries["Fingerprint"]],
    integrated["Projection"]["Upstream"], "the flows were projected with another Symmetries stage"];
   flows = integrated["Flows"];

   (* flows that must vanish: one check per class, the hint carries the reasons of the offenders *)
   Do[Module[{cl = c[[1]], name = c[[2]]},
      bad = Select[flows, #Class === cl && ! MatchQ[Lookup[#, "Samples"], {0 ..}] &];
      AppendTo[checks, assertCheck[name, bad === {}, {#Pair, #Frequency, #Momentum, N[#Samples, 4]} & /@ bad,
        "these flows must vanish but do not: " <> StringRiffle[(recordString[#] <> " (" <> #Reason <> ")") & /@ bad, "; "] <>
         ". The tag is wrong (the symmetry is not a symmetry of the ansatz and regulator), or the ansatz breaks it",
        ToString[Count[flows, KeyValuePattern["Class" -> cl]]] <> " " <> cl <> " flows vanish exactly (3 random rational points each)"]]],
    {c, {{"Protected", "Closure.protected"}, {"DeclaredNumeric", "Closure.declaredNumeric"}, {"Causality", "Closure.causality"}}}];

   (* neglected flows: reported, not asserted *)
   neg = Select[flows, #Class === "Neglected" &];
   warns = Join[warns, If[First[#Samples] === 0,
        infoCheck["INFO", "Closure.neglected " <> recordString[#], "flow vanishes at a valid random point (" <> #Reason <> ")"],
        infoCheck["WARN", "Closure.neglected " <> recordString[#],
         "generated, dropped: value " <> numString[First[#Samples]] <> " at one random valid point, not a bound (" <> #Reason <> ")"]] & /@ neg];

   (* declared monomials: d_s(comb) = Value / norm *)
   eqs = Function[r, Module[{norm, comb},
       {norm, comb} = splitCoefficient[r["GammaK"]];
       FRGAssert["Closure.normalisation", NumericQ[norm] && norm =!= 0 && Expand[norm comb - r["GammaK"]] === 0 && ! NumericQ[comb],
        {r["GammaK"], norm, comb}, "internal: the Gamma_k coefficient does not split into number x combination"];
       Join[r, <|"Normalisation" -> norm, "Flowing" -> comb, "RHS" -> r["Value"]/norm|>]]] /@ Select[flows, #Class === "Declared" &];
   AppendTo[checks, assertCheck["Closure.equations", eqs =!= {}, flows, "no declared monomial: nothing flows",
     ToString[Length[eqs]] <> " equations d_s(comb) = Value/norm: " <> StringRiffle[(inputForm[#Flowing] <> " (/" <> ToString[#Normalisation, InputForm] <> ")") & /@ eqs, ", "]]];

   AppendTo[warns, infoCheck["INFO", "Closure.structure", ToString[Count[symmetries["Classification"], KeyValuePattern["Class" -> "Structure"]]] <>
      " monomials vanish by structure (odd frequency powers of diagonal entries), not computed"]];
   warns = Join[warns, infoCheck["ASSUMPTION", "NonRenormalised[" <> ToString[#[[1]]] <> "]", #[[2, 2]] <> " (not checkable in v1: needs d_s Gamma3)"] & /@
      symmetries["NonRenormalised"]];
   AppendTo[warns, infoCheck["ASSUMPTION", "Closure.wardIdentities", "Ward identities are assumed (the Protected tags), not derived"]];

   report = Join[flows, Select[symmetries["Classification"], #Class === "Structure" &]];
   stage = stampStage[<|"Stage" -> "Closure", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {integrated["Fingerprint"], symmetries["Fingerprint"]},
      "Integrated" -> integrated, "Symmetries" -> symmetries, "Equations" -> eqs, "Report" -> report,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

outcome[r_] := Switch[r["Class"],
   "Declared", "equation for d_s(" <> inputForm[Last[splitCoefficient[r["GammaK"]]]] <> ")",
   "Neglected", If[First[r["Samples"]] === 0, "vanishes (dropped)", "generated, dropped: " <> numString[First[r["Samples"]]]],
   "Structure", "vanishes identically (not computed)",
   _, "flow = 0 (checked)"];

ShowClosureReport[c_Association] := Module[{rows},
   Print["Closure report. Declared: the coefficient is in the ansatz and gets a flow equation. Protected: a symmetry (Ward identity, \
user-declared) forbids a flow; checked to vanish. DeclaredNumeric: a pure number in the ansatz (e.g. the 1 in front of dt), its flow must vanish. \
Causality: Gamma_{phi phi} without response field stays 0. Neglected: generated by the flow but dropped by the truncation; value shown at one \
random point where the propagator poles are causal. Structure: odd frequency powers of a diagonal entry, zero by symmetry."];
   rows = {ToString[#Pair], monomialString[#], #Class, outcome[#]} & /@ c["Report"];
   Print[If[$FrontEnd === Null,
     StringRiffle[StringRiffle[MapThread[StringPadRight, {#, {18, 16, 16, 1}}], " | "] & /@ Prepend[rows, {"pair", "monomial", "class", "outcome"}], "\n"],
     Grid[Prepend[rows, {"pair", "monomial", "class", "outcome"}], Frame -> All, Alignment -> Left]]]];

End[];
EndPackage[];
