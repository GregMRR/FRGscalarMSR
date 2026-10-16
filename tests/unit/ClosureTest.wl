(* ::Package:: *)

(* Closure: normalisation split, Neglected values at valid points, Protected / Causality checks, aborts (KPZ) *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

kpzEq[c_] := SelectFirst[kpzClosure["Equations"], #Flowing === c &];
test["KPZ: two equations, d_s nuk and d_s Dk", Sort[kpzClosure["Equations"][[All, "Flowing"]]], Sort[{nuk, Dk}]];
test["KPZ: Gamma_k coefficient -2 Dk splits as -2 x Dk", {kpzEq[Dk]["Normalisation"], kpzEq[nuk]["Normalisation"]}, {-2, 1}];
test["KPZ: RHS = Value / normalisation", Expand[kpzEq[Dk]["RHS"] + kpzEq[Dk]["Value"]/2], 0];
test["KPZ: Neglected values sampled at valid points", AllTrue[Select[kpzIntegrated["Flows"], #Class === "Neglected" &], #SampleValid &], True];
test["KPZ: Galilean Ward identity and causality checked", SubsetQ[kpzClosure["Checks"][[All, "Check"]], {"Closure.protected", "Closure.causality"}], True];
test["KPZ: NonRenormalised printed as an assumption", MemberQ[kpzClosure["Warnings"][[All, "Check"]], "NonRenormalised[lamk]"], True];
test["KPZ: report covers every classified monomial", Length[kpzClosure["Report"]], Length[kpzSymmetries["Classification"]]];

(* a protected term with nonzero flow: w^2 in Gamma_{hbar hbar} is generated, tagging it Protected must abort (design section 10) *)
badSym = quiet[DeclareSymmetries[kpzAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift"],
    Protected[{hbar, hbar}, "Frequency" -> 2, "Momentum" -> 0, "Reason" -> "wrong on purpose"],
    Neglected[_, "Reason" -> "beyond the ansatz"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
expectFail["protected term with nonzero flow", CheckClosure[IntegrateFrequencies[ProjectOnMonomials[kpzLoops, badSym]], badSym], "Closure.protected"];
expectFail["flows projected with other symmetries", CheckClosure[kpzIntegrated, badSym], "Closure.sameSymmetries"];
expectFail["wrong input stage", CheckClosure[kpzLoops, kpzSymmetries], "Stage.input"];
