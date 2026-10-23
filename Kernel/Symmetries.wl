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

specMatchQ[All, _Integer] := True;
specMatchQ[l_List, x_] := MemberQ[l, x];
specMatchQ[i_Integer, x_] := i === x;
specMatchQ[_, _] := False;

pairMatchQ[pat_, {a_, b_}] := MatchQ[{a, b}, pat] || MatchQ[{b, a}, pat];

tagRecord[tag_, kind_] := Module[{pat = First[tag], o = Association[Rest[List @@ tag]]},
   FRGAssert["Symmetries.tagOptions", SubsetQ[{"Frequency", "Momentum", "Reason"}, Keys[o]], tag,
    kind <> "[pairPattern, \"Frequency\" -> .., \"Momentum\" -> .., \"Reason\" -> \"..\"]"];
   FRGAssert["Symmetries.tagMomentum", MatchQ[Lookup[o, "Momentum", All], All | "AllMomenta" | _Integer?NonNegative | {__Integer?NonNegative}], tag,
    "\"Momentum\" -> n, {n, ...}, All (every order up to the truncation, at p -> 0) or \"AllMomenta\" (the whole slot at generic p)"];
   FRGAssert["Symmetries.reason", StringQ[Lookup[o, "Reason", None]] && Lookup[o, "Reason"] =!= "", tag, "every tag needs a \"Reason\" -> \"...\""];
   <|"Kind" -> kind, "Pattern" -> pat, "Frequency" -> Lookup[o, "Frequency", All], "Momentum" -> Lookup[o, "Momentum", All],
    "Reason" -> o["Reason"], "Tag" -> tag|>];

tagMatchQ[t_, {a_, b_}, m_, n_] := pairMatchQ[t["Pattern"], {a, b}] && specMatchQ[t["Frequency"], m] && specMatchQ[t["Momentum"], n];
(* a tag on the whole slot (pair, m) at generic p *)
allMomentaQ[t_, {a_, b_}, m_] := t["Momentum"] === "AllMomenta" && pairMatchQ[t["Pattern"], {a, b}] && specMatchQ[t["Frequency"], m];
(* a tag naming explicit orders n (not All): the only ones that can reach a functional slot (design/MOMENTUMFUNCTIONS.md B2) *)
explicitOrders[t_] := Replace[t["Momentum"], {i_Integer :> {i}, l_List :> l, _ -> {}}];

(* coefficient of (I nu)^m (p^2)^n in a Gamma2 entry written with scal[p, p] *)
gammaCoefficient[e_, nu_, p_, m_, n_] := Expand[Coefficient[Coefficient[Expand[e /. scal[p, p] -> pp2], pp2, n], nu, m]/I^m];

DeclareSymmetries[ansatz_, args___] := Block[{$FRGModule = "Symmetries"},
  Module[{t0 = AbsoluteTime[], theory, fields, phis, pos, opts, trunc, mf, mp, tags, prot, negl, nonren, v, g2, pairs, cls, checks = {}, warns = {},
    used = {}, stage, cps, bad, declOutside, fns, fnPat, slotRows, classify, fnSlots, siblings, nAll},
   validateStage[ansatz, "Ansatz"];
   theory = ansatz["Theory"]; fields = theory["Fields"]; phis = theory["Scalars"]; cps = ansatz["Couplings"];
   fns = Lookup[ansatz, "MomentumFunctions", {}]; fnPat = Alternatives @@ fns;
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
   AppendTo[checks, assertCheck["Symmetries.truncation", True, mp, "",
     "truncation: (I w)^" <> ToString[mf] <> ", (q^2)^" <> ToString[mp]]];

   prot = tagRecord[#, "Protected"] & /@ Cases[tags, _Protected];
   negl = tagRecord[#, "Neglected"] & /@ Cases[tags, _Neglected];
   nonren = Cases[tags, _NonRenormalised];
   FRGAssert["Symmetries.nonRenormalisedFunction", FreeQ[nonren[[All, 1]], fnPat], Select[nonren, MemberQ[fns, First[#]] &],
    "NonRenormalised takes constant couplings: a momentum function flows by its own equation"];
   FRGAssert["Symmetries.nonRenormalised", AllTrue[nonren, MatchQ[#, NonRenormalised[c_Symbol, "Reason" -> r_String] /; MemberQ[cps, c] && r =!= ""] &], nonren,
    "NonRenormalised[coupling, \"Reason\" -> \"...\"] with a coupling of the ansatz"];

   (* Gamma_k coefficients (regulator off) from the vertices; silent: BuildVertices prints its own table *)
   v = Block[{Print}, BuildVertices[ansatz]];
   pairs = Select[Tuples[fields, 2], pos[#[[1]]] <= pos[#[[2]]] &];
   classify[a_, b_, m_, n_, coef_] := Module[{cl, reason, pt, nt},
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
       <|"Pair" -> {a, b}, "Frequency" -> m, "Momentum" -> n, "Class" -> cl, "GammaK" -> coef, "Reason" -> reason|>];

   (* one slot (pair, m): ordinary rows (I w)^m (q^2)^n, n <= mp; a slot carrying a momentum function F (design/MOMENTUMFUNCTIONS.md A4,
      B1, B2): its orders below n0 are ordinary rows, the orders >= n0 are one row of class Functional; a function-free slot with an
      "AllMomenta" tag also gets one row at generic p (B3) *)
   slotRows[a_, b_, m_] := Module[{e = VertexValue[v, {a, b}, {{nu, p}}], cm, fa, f, t, n0, pref, rest, hi, ft, atZero, am},
      cm = Expand[Coefficient[Expand[e /. scal[p, p] -> pp2], nu, m]/I^m];
      fa = Union[Cases[cm, fnPat[_], {0, Infinity}][[All, 0]]];
      am = Select[Join[prot, negl], allMomentaQ[#, {a, b}, m] &];
      If[fa === {},
       Join[Table[classify[a, b, m, n, Coefficient[cm, pp2, n]], {n, 0, mp}],
        If[am === {} || ! IntersectingQ[{a, b}, theory["ResponseFields"]] || (a === b && OddQ[m]), {},
         FRGAssert["Symmetries.contradiction", First[am]["Kind"] === "Neglected" || FreeQ[cm, Alternatives @@ cps], {{a, b}, m, cm, First[am]["Tag"]},
          "this slot carries a coupling of the ansatz and is also Protected at every momentum (its flow would have to vanish): drop the tag or the coupling"];
         used = Join[used, {First[am]["Tag"]}];
         {<|"Pair" -> {a, b}, "Frequency" -> m, "Momentum" -> "AllMomenta", "Class" -> First[am]["Kind"], "GammaK" -> cm, "Reason" -> First[am]["Reason"]|>}]],
       FRGAssert["Symmetries.functionalSlot", Length[fa] == 1, {{a, b}, m, fa}, "one function per slot: two functions in one slot are one function, keep one"];
       f = First[fa];
       t = Expand[D[cm /. f[pp2] -> fz, fz]];
       n0 = Exponent[t, pp2]; pref = Expand[t/pp2^n0];
       rest = Expand[cm /. f[pp2] -> 0];
       hi = Select[Range[n0, Max[n0, Exponent[rest, pp2]]], Coefficient[rest, pp2, #] =!= 0 &];
       FRGAssert["Symmetries.functionalSlot", hi === {}, {{a, b}, m, f, "orders " -> hi},
        "a constant term at an order >= n0 of the function's slot is F's job (e.g. -nk hbar lap[h] next to -nk Fnu hbar lap[h]): drop it"];
       FRGAssert["Symmetries.contradiction", am === {}, {{a, b}, m, f, am[[All, "Tag"]]},
        "an \"AllMomenta\" tag on a slot that carries the momentum function " <> ToString[f] <> ": use \"Momentum\" -> " <> ToString[n0] <> " (Protected) for F(0), or Related"];
       (* tags naming explicit orders >= n0: only Protected at exactly n0 (the loop part of the flow at p = 0 vanishes) *)
       ft = Select[Join[prot, negl], pairMatchQ[#Pattern, {a, b}] && specMatchQ[#Frequency, m] && AnyTrue[explicitOrders[#], # >= n0 &] &];
       FRGAssert["Symmetries.functionalTag", AllTrue[ft, #Kind === "Protected" && explicitOrders[#] === {n0} &], {{a, b}, m, f, ft[[All, "Tag"]]},
        "only Protected[pair, \"Frequency\" -> " <> ToString[m] <> ", \"Momentum\" -> " <> ToString[n0] <>
         ", ...] is allowed on a functional slot (F(0) does not flow); Neglected or higher orders have no meaning there"];
       used = Join[used, ft[[All, "Tag"]]];
       atZero = If[ft === {}, None, First[ft]["Reason"]];
       Append[Table[classify[a, b, m, n, Coefficient[rest, pp2, n]], {n, 0, Min[mp, n0 - 1]}],
        <|"Pair" -> {a, b}, "Frequency" -> m, "Momentum" -> n0, "Class" -> "Functional", "GammaK" -> pref f, "Function" -> f,
         "ProtectedAtZero" -> atZero,
         "Reason" -> "flow equation for the function " <> ToString[f] <> "(p^2) at generic p" <> If[atZero === None, "", "; F(0) does not flow (" <> atZero <> ")"]|>]]];
   cls = Flatten[Table[slotRows[pr[[1]], pr[[2]], m], {pr, pairs}, {m, 0, mf}], 2];

   (* design/MOMENTUMFUNCTIONS.md A4: every declared function sits in exactly one slot with a nonzero Gamma2 *)
   fnSlots = AssociationMap[Function[f, Select[cls, Lookup[#, "Function", None] === f &][[All, {"Pair", "Frequency"}]]], fns];
   If[fns =!= {},
    bad = Select[fns, Length[fnSlots[#]] =!= 1 &];
    FRGAssert["Symmetries.functionSlot", bad === {}, <|# -> fnSlots[#] & /@ bad|>,
     "every momentum function sits in exactly one Gamma2 slot: a function whose term vanishes in Gamma2 (total derivative, odd w on a diagonal \
entry, e.g. Fnu div[hbar grad[h]]) must be dropped; a function in two slots must become two functions related by Related[...]"];
    AppendTo[checks, assertCheck["Symmetries.functionSlot", True, Null, "",
      "one slot per momentum function: " <> StringRiffle[(ToString[#] <> " in " <> ToString[fnSlots[#][[1, "Pair"]]] <> " (I w)^" <> ToString[fnSlots[#][[1, "Frequency"]]]) & /@ fns, ", "]]];
    AppendTo[warns, infoCheck["INFO", "Symmetries.functionalTruncation", "\"TruncationOrders\" \"Momentum\" applies to function-free slots and to the orders below n0 of functional slots"]];
    (* B4: siblings of a functional slot that are computed at p -> 0 only (e.g. the I w slot of LO KPZ, f^lambda = 1) *)
    siblings = Select[cls, Function[r, MemberQ[{"Declared", "DeclaredNumeric", "Protected", "Neglected"}, r["Class"]] && IntegerQ[r["Momentum"]] &&
          AnyTrue[Values[fnSlots], Function[sl, sl[[1, "Pair"]] === r["Pair"] && sl[[1, "Frequency"]] =!= r["Frequency"]]] &&
          ! MemberQ[cls, KeyValuePattern[{"Pair" -> r["Pair"], "Frequency" -> r["Frequency"], "Momentum" -> "AllMomenta"}]]]];
    siblings = DeleteDuplicates[siblings[[All, {"Pair", "Frequency"}]]];
    warns = Join[warns, infoCheck["ASSUMPTION", "Symmetries.siblingAtZeroMomentum", ToString[#Pair] <> " (I w)^" <> ToString[#Frequency] <>
         " computed at p -> 0 only (add Neglected[" <> ToString[#Pair] <> ", \"Frequency\" -> " <> ToString[#Frequency] <> ", \"Momentum\" -> \"AllMomenta\", ...] to see its size)"] & /@ siblings]];
   nAll = Count[cls, KeyValuePattern["Momentum" -> "AllMomenta"]];
   If[nAll > 0, AppendTo[checks, assertCheck["Symmetries.allMomenta", True, Null, "", ToString[nAll] <> " slot(s) tagged \"AllMomenta\": computed at generic p"]]];
   bad = Select[cls, #Class === "Untagged" &];
   AppendTo[checks, assertCheck["Symmetries.allTagged", bad === {}, {#Pair, "(I w)^" <> ToString[#Frequency], "(q^2)^" <> ToString[#Momentum]} & /@ bad,
     "every monomial up to the truncation orders must be declared in the ansatz, Protected or Neglected", ToString[Length[cls]] <> " monomials classified"]];

   (* declared monomials of the ansatz beyond the truncation orders would be silently ignored *)
   (* functional slots: only their orders below n0 (the function is a single row at any order) *)
   declOutside = Select[Flatten[Table[{pr, m, n}, {pr, pairs}, {m, 0, 2 mf + 4}, {n, 0, Max[2, Max[ansatz["Monomials"][[All, "MomentumOrder"]]]/2]}], 2],
     (#[[2]] > mf || #[[3]] > mp) && With[{fr = SelectFirst[cls, Function[r, r["Class"] === "Functional" && r["Pair"] === #[[1]] && r["Frequency"] === #[[2]]]]},
        If[MissingQ[fr], gammaCoefficient[VertexValue[v, #[[1]], {{nu, p}}] /. fnPat[_] -> 0, nu, p, #[[2]], #[[3]]] =!= 0 ||
          (#[[2]] > mf && ! FreeQ[Coefficient[Expand[VertexValue[v, #[[1]], {{nu, p}}]], nu, #[[2]]], fnPat]),
         #[[3]] < fr["Momentum"] && gammaCoefficient[VertexValue[v, #[[1]], {{nu, p}}] /. fnPat[_] -> 0, nu, p, #[[2]], #[[3]]] =!= 0]] &];
   AppendTo[checks, assertCheck["Symmetries.truncationCoversAnsatz", declOutside === {}, declOutside,
     "the ansatz contains Gamma2 monomials beyond \"TruncationOrders\": raise the orders", "the truncation contains every monomial of the ansatz"]];

   Scan[If[! MemberQ[used, #Tag], AppendTo[warns, infoCheck["WARN", "Symmetries.unusedTag", inputForm[#Tag] <> " matches no monomial (precedence)"]]] &, Join[prot, negl]];
   warns = Join[warns, infoCheck["ASSUMPTION", "NonRenormalised[" <> ToString[#[[1]]] <> "]", #[[2, 2]] <> " (not checkable in v1: needs d_s Gamma3)"] & /@ nonren];
   AppendTo[warns, infoCheck["ASSUMPTION", "Symmetries.wardIdentities", "Protected tags are user-declared Ward identities, not derived"]];

   stage = stampStage[<|"Stage" -> "Symmetries", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {ansatz["Fingerprint"]},
      "Fields" -> fields, "Classification" -> cls, "TruncationOrders" -> trunc, Sequence @@ If[fns === {}, {}, {"MomentumFunctions" -> fns}],
      "Tags" -> Join[prot, negl], "NonRenormalised" -> nonren, "Dimension" -> theory["Dimension"],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowClassification[stage];
   PrintCard[stage];
   stage]];

ShowClassification[s_Association] := Module[{rows},
   rows = {ToString[#Pair], "(I w)^" <> ToString[#Frequency] <> Switch[#Momentum, "AllMomenta", " all p", _, If[#Class === "Functional",
          " (q^2)^{>=" <> ToString[#Momentum] <> "}", " (q^2)^" <> ToString[#Momentum]]], #Class,
       If[#GammaK === 0, "", If[#Class === "Functional", With[{f = Lookup[#, "Function", None]}, displayString[#GammaK /. f -> f[Symbol["FRGscalarMSR`Symmetries`Private`p"]^2]]] , inputForm[#GammaK]]], #Reason} & /@ s["Classification"];
   Print["Monomial classification (Gamma2_ab, a before b in field order; precedence Protected > declared > Neglected):"];
   Print[If[$FrontEnd === Null,
     StringRiffle[StringRiffle[MapThread[StringPadRight, {#, {18, 16, 16, 18, 1}}], " | "] & /@ Prepend[rows, {"pair", "monomial", "class", "Gamma_k coeff.", "reason"}], "\n"],
     Grid[Prepend[rows, {"pair", "monomial", "class", "Gamma_k coeff.", "reason"}], Frame -> All, Alignment -> Left]]]];

End[];
EndPackage[];
