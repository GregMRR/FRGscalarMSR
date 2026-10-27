(* ::Package:: *)

(* Symmetries: classification, precedence, aborts; KPZ through IntegrateFrequencies *)

sTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
sAnsatz = quiet[DefineAnsatz[sTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
kpzTags = Sequence[Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
   Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"]];
sym = quiet[DeclareSymmetries[sAnsatz, kpzTags, "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
cl[pair_, m_, n_] := SelectFirst[sym["Classification"], #Pair === pair && #Frequency == m && #Momentum == n &]["Class"];

test["KPZ: 18 monomials", Length[sym["Classification"]], 18];
test["KPZ: h h is causality", cl[{h, h}, 0, 1], "Causality"];
test["KPZ: nu q^2 declared", cl[{h, hbar}, 0, 1], "Declared"];
test["KPZ: i w coefficient 1 is Protected (precedence over numeric declared)", cl[{h, hbar}, 1, 0], "Protected"];
test["KPZ: odd w of hbar hbar vanishes by structure", cl[{hbar, hbar}, 1, 1], "Structure"];
test["KPZ: D declared", cl[{hbar, hbar}, 0, 0], "Declared"];
test["KPZ: w^2 hbar hbar neglected", cl[{hbar, hbar}, 2, 0], "Neglected"];
test["NonRenormalised printed as an assumption", MemberQ[sym["Warnings"][[All, "Check"]], "NonRenormalised[lamk]"], True];

kpzFlows = quiet[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[#1, #2, #3], sym] & @@
     Module[{v = quiet[BuildVertices[sAnsatz]], r = quiet[DefineRegulator[sAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}]]},
      {v, quiet[BuildPropagator[v, r]], r}]]];
test["KPZ: Galilean i w coefficient flows to exactly 0 (Protected samples)",
  SelectFirst[kpzFlows["Flows"], #Pair === {h, hbar} && #Frequency == 1 && #Momentum == 0 &]["Samples"], {0, 0, 0}];
test["KPZ: causality flows vanish", Union[Flatten[Select[kpzFlows["Flows"], #Class === "Causality" &][[All, "Samples"]]]], {0}];
test["KPZ: nu and D flows are symbolic and nonzero", Length[Select[kpzFlows["Flows"], #Class === "Declared" && Lookup[#, "Value"] =!= 0 &]], 2];

(* aborts *)
trunc = "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>;
expectFail["untagged monomial", DeclareSymmetries[sAnsatz, Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "x"], trunc], "Symmetries.allTagged"];
expectFail["Protected pattern hitting a declared noise term", DeclareSymmetries[sAnsatz, Protected[{hbar, _}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "x"],
   Neglected[_, "Reason" -> "y"], trunc], "Symmetries.contradiction"];
expectFail["declared + Protected contradiction", DeclareSymmetries[sAnsatz, Protected[{hbar, h}, "Momentum" -> 1, "Reason" -> "x"], Neglected[_, "Reason" -> "y"], trunc],
  "Symmetries.contradiction"];
test["p^4 truncation accepted (design/MOMENTUMPROJECTION.md)",
  Max[quiet[DeclareSymmetries[sAnsatz, Neglected[_, "Reason" -> "y"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 2|>]]["Classification"][[All, "Momentum"]]], 2];
expectFail["missing reason", DeclareSymmetries[sAnsatz, Neglected[_], trunc], "Symmetries.reason"];
expectFail["truncation below the ansatz", DeclareSymmetries[sAnsatz, Neglected[_, "Reason" -> "y"], "TruncationOrders" -> <|"Frequency" -> 0, "Momentum" -> 1|>],
  "Symmetries.truncationCoversAnsatz"];
expectFail["NonRenormalised of a non-coupling", DeclareSymmetries[sAnsatz, Neglected[_, "Reason" -> "y"], NonRenormalised[zz, "Reason" -> "r"], trunc], "Symmetries.nonRenormalised"];
expectFail["missing truncation orders", DeclareSymmetries[sAnsatz, Neglected[_, "Reason" -> "y"]], "Symmetries.arguments"];

(* LO branch: declarations only (design/LOPIPELINE.md L7); BuildLOEquations checks them *)
sLOAnsatz = quiet[DefineAnsatz[sTheory, "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
    "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]];
sLOProt = Sequence[Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift symmetry"]];
sLO = quiet[DeclareSymmetries[sLOAnsatz, sLOProt, NonRenormalised[lamk, "Reason" -> "Galilean invariance"]]];
test["LO: stage without classification", {sLO["Pipeline"], KeyExistsQ[sLO, "Classification"], KeyExistsQ[sLO, "TruncationOrders"]}, {"LO", False, False}];
test["LO: tags recorded", {sLO["Tags"][[All, "Kind"]], sLO["NonRenormalised"][[All, 1]]}, {{"Protected", "Protected"}, {lamk}}];
expectFail["LO: no truncation orders", DeclareSymmetries[sLOAnsatz, sLOProt, trunc], "Symmetries.arguments"];
expectFail["LO: a Protected tag names one w order and one p^2 order", DeclareSymmetries[sLOAnsatz, Protected[{hbar, h}, "Reason" -> "x"]], "Symmetries.tagOptions"];
expectFail["LO: Neglected drops whole entries", DeclareSymmetries[sLOAnsatz, Neglected[{hbar, hbar}, "Momentum" -> 1, "Reason" -> "x"]], "Symmetries.tagOptions"];
