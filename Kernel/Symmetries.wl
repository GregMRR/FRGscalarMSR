(* ::Package:: *)

(* Symmetries.wl -- DeclareSymmetries: classify every monomial (I w)^m (q^2)^n of every Gamma2 entry up to the truncation
   orders as declared / Protected / Neglected / vanishing by structure, before any flow is computed.
   Objective, assumptions and checks: docs/MODULES.md section "Symmetries". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

DeclareSymmetries::usage = "DeclareSymmetries[ansatz, tags..., \"TruncationOrders\" -> <|\"Frequency\" -> mf, \"Momentum\" -> mp|>] \
classifies every monomial (I w)^m (q^2)^n, m <= mf, n <= mp, of every Gamma2 entry. Tags: Protected[pairPattern, opts], \
Neglected[pairPattern, opts], NonRenormalised[coupling, \"Reason\" -> ...]. Precedence: Protected > declared > Neglected; an untagged \
monomial aborts. Returns the \"Symmetries\" stage.";
Neglected::usage = "Neglected[pairPattern, \"Frequency\" -> m | {m..} | All, \"Momentum\" -> n | {n..} | All, \"Reason\" -> text]: \
the flow of these monomials is computed and reported as a warning, then dropped.";
NonRenormalised::usage = "NonRenormalised[coupling, \"Reason\" -> text]: a vertex coupling assumed not to flow (not checkable in v1: needs d_s Gamma3).";
ShowClassification::usage = "ShowClassification[symmetries] prints the classification table of all monomials.";
MonomialClasses::usage = "Classes of a monomial: \"Declared\" (equation), \"DeclaredNumeric\" (numeric coefficient: flow must vanish), \
\"Protected\" (flow must vanish), \"Neglected\" (flow reported), \"Causality\" (Gamma_{phi phi}: flow must vanish), \
\"Structure\" (odd frequency power of a diagonal entry: vanishes identically, not computed).";

Begin["`Private`"];

(* Protected is a System symbol (attribute); the tag is used in the unambiguous form Protected[pattern, opts] *)
SetAttributes[{Neglected, NonRenormalised}, Protected];

specMatchQ[All, _] := True;
specMatchQ[l_List, x_] := MemberQ[l, x];
specMatchQ[i_Integer, x_] := i === x;
specMatchQ[_, _] := False;

pairMatchQ[pat_, {a_, b_}] := MatchQ[{a, b}, pat] || MatchQ[{b, a}, pat];

tagRecord[tag_, kind_] := Module[{pat = First[tag], o = Association[Rest[List @@ tag]]},
   FRGAssert["Symmetries.tagOptions", SubsetQ[{"Frequency", "Momentum", "Reason"}, Keys[o]], tag,
    kind <> "[pairPattern, \"Frequency\" -> .., \"Momentum\" -> .., \"Reason\" -> \"..\"]"];
   FRGAssert["Symmetries.reason", StringQ[Lookup[o, "Reason", None]] && Lookup[o, "Reason"] =!= "", tag, "every tag needs a \"Reason\" -> \"...\""];
   <|"Kind" -> kind, "Pattern" -> pat, "Frequency" -> Lookup[o, "Frequency", All], "Momentum" -> Lookup[o, "Momentum", All],
    "Reason" -> o["Reason"], "Tag" -> tag|>];

tagMatchQ[t_, {a_, b_}, m_, n_] := pairMatchQ[t["Pattern"], {a, b}] && specMatchQ[t["Frequency"], m] && specMatchQ[t["Momentum"], n];

(* coefficient of (I nu)^m (p^2)^n in a Gamma2 entry written with scal[p, p] *)
gammaCoefficient[e_, nu_, p_, m_, n_] := Expand[Coefficient[Coefficient[Expand[e /. scal[p, p] -> pp2], pp2, n], nu, m]/I^m];

DeclareSymmetries[ansatz_, args___] := Block[{$FRGModule = "Symmetries"},
  Module[{t0 = AbsoluteTime[], theory, fields, phis, pos, opts, trunc, mf, mp, tags, prot, negl, nonren, v, g2, pairs, cls, checks = {}, warns = {},
    used = {}, stage, cps, bad, declOutside},
   validateStage[ansatz, "Ansatz"];
   theory = ansatz["Theory"]; fields = theory["Fields"]; phis = theory["Scalars"]; cps = ansatz["Couplings"];
   pos[a_] := First[FirstPosition[fields, a]];
   opts = Cases[{args}, (Rule | RuleDelayed)[_String, _]];
   tags = Complement[{args}, opts];
   FRGAssert["Symmetries.arguments", AllTrue[tags, MatchQ[#, _Protected | _Neglected | _NonRenormalised] &] && MatchQ[opts, {"TruncationOrders" -> _}],
    Select[{args}, ! MatchQ[#, _Protected | _Neglected | _NonRenormalised | ("TruncationOrders" -> _)] &],
    "arguments: Protected[...], Neglected[...], NonRenormalised[...] and \"TruncationOrders\" -> <|\"Frequency\" -> mf, \"Momentum\" -> mp|> (required)"];
   trunc = Lookup[opts, "TruncationOrders"];
   FRGAssert["Symmetries.truncationOrders", AssociationQ[trunc] && Sort[Keys[trunc]] === {"Frequency", "Momentum"} && AllTrue[Values[trunc], IntegerQ[#] && # >= 0 &],
    trunc, "\"TruncationOrders\" -> <|\"Frequency\" -> 2, \"Momentum\" -> 1|> (non-negative integers)"];
   {mf, mp} = {trunc["Frequency"], trunc["Momentum"]};
   AppendTo[checks, assertCheck["Symmetries.momentumOrder", mp <= 1, mp,
     "momentum truncation above p^2 is not supported in v1: the p-Laplacian projection is exact only for p^2 (E2, docs/FUTURE.md)",
     "truncation: (I w)^" <> ToString[mf] <> ", (q^2)^" <> ToString[mp]]];

   prot = tagRecord[#, "Protected"] & /@ Cases[tags, _Protected];
   negl = tagRecord[#, "Neglected"] & /@ Cases[tags, _Neglected];
   nonren = Cases[tags, _NonRenormalised];
   FRGAssert["Symmetries.nonRenormalised", AllTrue[nonren, MatchQ[#, NonRenormalised[c_Symbol, "Reason" -> r_String] /; MemberQ[cps, c] && r =!= ""] &], nonren,
    "NonRenormalised[coupling, \"Reason\" -> \"...\"] with a coupling of the ansatz"];

   (* Gamma_k coefficients (regulator off) from the vertices; silent: BuildVertices prints its own table *)
   v = Block[{Print}, BuildVertices[ansatz]];
   pairs = Select[Tuples[fields, 2], pos[#[[1]]] <= pos[#[[2]]] &];
   cls = Flatten[Table[Module[{a = pr[[1]], b = pr[[2]], e, coef, cl, reason, pt, nt},
       e = VertexValue[v, {a, b}, {{nu, p}}];
       coef = gammaCoefficient[e, nu, p, m, n];
       pt = Select[prot, tagMatchQ[#, {a, b}, m, n] &]; nt = Select[negl, tagMatchQ[#, {a, b}, m, n] &];
       {cl, reason} = Which[
         ! IntersectingQ[{a, b}, theory["ResponseFields"]], {"Causality", "Gamma_{phi phi} = 0 (causality); its flow is checked to vanish"},
         a === b && OddQ[m], {"Structure", "Gamma2_aa(w) = Gamma2_aa(-w): odd powers of w vanish identically"},
         coef =!= 0 && ! NumericQ[coef] && pt =!= {}, FRGAssert["Symmetries.contradiction", False, {{a, b}, m, n, coef, pt[[All, "Tag"]]},
           "this monomial carries a coupling of the ansatz and is also Protected (its flow would have to vanish): drop the tag or the coupling"],
         pt =!= {}, (used = Join[used, pt[[All, "Tag"]]]; {"Protected", First[pt]["Reason"]}),
         coef =!= 0 && NumericQ[coef], {"DeclaredNumeric", "numeric coefficient " <> ToString[coef] <> " in the ansatz: its flow must vanish"},
         coef =!= 0, {"Declared", "flow equation for " <> inputForm[coef]},
         nt =!= {}, (used = Join[used, nt[[All, "Tag"]]]; {"Neglected", First[nt]["Reason"]}),
         True, {"Untagged", ""}];
       <|"Pair" -> {a, b}, "Frequency" -> m, "Momentum" -> n, "Class" -> cl, "GammaK" -> coef, "Reason" -> reason|>],
      {pr, pairs}, {m, 0, mf}, {n, 0, mp}], 2];
   bad = Select[cls, #Class === "Untagged" &];
   AppendTo[checks, assertCheck["Symmetries.allTagged", bad === {}, {#Pair, "(I w)^" <> ToString[#Frequency], "(q^2)^" <> ToString[#Momentum]} & /@ bad,
     "every monomial up to the truncation orders must be declared in the ansatz, Protected or Neglected", ToString[Length[cls]] <> " monomials classified"]];

   (* declared monomials of the ansatz beyond the truncation orders would be silently ignored *)
   declOutside = Select[Flatten[Table[{pr, m, n}, {pr, pairs}, {m, 0, 2 mf + 4}, {n, 0, 2}], 2],
     (#[[2]] > mf || #[[3]] > mp) && gammaCoefficient[VertexValue[v, #[[1]], {{nu, p}}], nu, p, #[[2]], #[[3]]] =!= 0 &];
   AppendTo[checks, assertCheck["Symmetries.truncationCoversAnsatz", declOutside === {}, declOutside,
     "the ansatz contains Gamma2 monomials beyond \"TruncationOrders\": raise the orders", "the truncation contains every monomial of the ansatz"]];

   Scan[If[! MemberQ[used, #Tag], AppendTo[warns, infoCheck["WARN", "Symmetries.unusedTag", inputForm[#Tag] <> " matches no monomial (precedence)"]]] &, Join[prot, negl]];
   warns = Join[warns, infoCheck["ASSUMPTION", "NonRenormalised[" <> ToString[#[[1]]] <> "]", #[[2, 2]] <> " (not checkable in v1: needs d_s Gamma3)"] & /@ nonren];
   AppendTo[warns, infoCheck["ASSUMPTION", "Symmetries.wardIdentities", "Protected tags are user-declared Ward identities, not derived"]];

   stage = stampStage[<|"Stage" -> "Symmetries", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {ansatz["Fingerprint"]},
      "Fields" -> fields, "Classification" -> cls, "TruncationOrders" -> trunc,
      "Tags" -> Join[prot, negl], "NonRenormalised" -> nonren, "Dimension" -> theory["Dimension"],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowClassification[stage];
   PrintCard[stage];
   stage]];

ShowClassification[s_Association] := Module[{rows},
   rows = {ToString[#Pair], "(I w)^" <> ToString[#Frequency] <> " (q^2)^" <> ToString[#Momentum], #Class,
       If[#GammaK === 0, "", inputForm[#GammaK]], #Reason} & /@ s["Classification"];
   Print["Monomial classification (Gamma2_ab, a before b in field order; precedence Protected > declared > Neglected):"];
   Print[If[$FrontEnd === Null,
     StringRiffle[StringRiffle[MapThread[StringPadRight, {#, {18, 16, 16, 18, 1}}], " | "] & /@ Prepend[rows, {"pair", "monomial", "class", "Gamma_k coeff.", "reason"}], "\n"],
     Grid[Prepend[rows, {"pair", "monomial", "class", "Gamma_k coeff.", "reason"}], Frame -> All, Alignment -> Left]]]];

End[];
EndPackage[];
