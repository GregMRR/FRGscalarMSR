(* ::Package:: *)

(* ExternalMomentum (BuildLOEquations, design/LOPIPELINE.md step 5): the KPZ LO equations, their records, every check, and the
   pipeline guards of the DE-only modules. The physics (equality with BMWflow's LO equations) is in tests/acceptance/KPZLO.wl. *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];

emRec[f_] := SelectFirst[kpzloEquations["Equations"], #Function === f &];
test["KPZ LO: stage and pipeline", {kpzloEquations["Stage"], kpzloEquations["Pipeline"]}, {"Equations", "LO"}];
test["KPZ LO: one equation per function", Sort[kpzloEquations["Equations"][[All, "Function"]]], {fD, fnu}];
test["KPZ LO: fnu sits in Gamma2(h, hbar) = nuk p^2 fnu(p^2), divided by p^2",
  emRec[fnu][[{"Pair", "Reference", "Normalisation", "P2Divided", "GammaK"}]], <|"Pair" -> {h, hbar}, "Reference" -> nuk, "Normalisation" -> 1,
   "P2Divided" -> 1, "GammaK" -> nuk p^2 fnu[p^2]|>];
test["KPZ LO: fD sits in Gamma2(hbar, hbar) = -2 Dk fD(p^2)",
  emRec[fD][[{"Pair", "Reference", "Normalisation", "P2Divided"}]], <|"Pair" -> {hbar, hbar}, "Reference" -> Dk, "Normalisation" -> -2, "P2Divided" -> 0|>];
test["KPZ LO: integrands in p, q, CosTheta only (no vector, no Delta, no atom)",
  FreeQ[kpzloEquations["Equations"][[All, "Integrand"]], scal | wedge | \[CapitalDelta] | GammaAtom | DeltaCoefficient | dsR], True];
test["KPZ LO: f at q^2 and at Q^2 in the flows, only at q^2 in the eta rows",
  {Union[Cases[emRec[fnu]["Integrand"], fnu[x_] :> x, {0, Infinity}]], Union[Cases[emRec[fnu]["EtaIntegrand"], (fnu | Derivative[1][fnu])[x_] :> x, {0, Infinity}]]},
  {Sort[{q^2, p^2 + 2 CosTheta p q + q^2}], {q^2}}];
test["KPZ LO: the eta rows are the exact p -> 0 limit of the flows (angular average at p = 0)",
  ExactlyZeroQ[(emRec[fD]["Integrand"] /. p -> 0) - emRec[fD]["EtaIntegrand"]], True];
test["KPZ LO at F = 1: the eta rows are the DE flows of Dk and nuk (tests/reference/KPZSetup.wl)",
  ExactlyZeroQ[({emRec[fD]["EtaIntegrand"], emRec[fnu]["EtaIntegrand"]} /. {fD -> (1 &), fnu -> (1 &)}) -
    {EquationOf[kpzEquations, Dk]["Integrand"], EquationOf[kpzEquations, nuk]["Integrand"]}], True];
test["KPZ LO: card", Sort[kpzloEquations["Checks"][[All, "Check"]]], Sort[{"ExternalMomentum.isotropicAngles", "FrequencyIntegration.regularAtCoincidence",
    "ExternalMomentum.zeroEntries", "ExternalMomentum.ward", "ExternalMomentum.firstDerivativesOnly", "ExternalMomentum.dividedOrders"}]];

(* every check fires on a model that breaks it *)
emTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
expectFail["a Gamma4 tadpole c4 hbar h grad[h].grad[h] gives a p^0 term in the row divided by p^2",
  loEquationsOf[emTheory, kpzloDensity + c4 hbar h grad[h].grad[h], {nuk, lamk, Dk, c4}, {fD -> Dk, fnu -> nuk}, kpzloTags, kpzloRegs], "ExternalMomentum.dividedOrders"];
expectFail["a constant vertex c3 hbar h^2: the eta row of fnu needs r'' and F''",
  loEquationsOf[emTheory, hbar (dt[h] - nuk fn[fnu, lap[h]] - c3 h^2) - Dk hbar fn[fD, hbar], {nuk, c3, Dk}, {fD -> Dk, fnu -> nuk}, {}, kpzloRegs],
  "ExternalMomentum.firstDerivativesOnly"];
expectFail["a non-Galilean vertex hbar h lap[h] breaks the declared Ward identity at p = 0",
  loEquationsOf[emTheory, hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk h lap[h]) - Dk hbar fn[fD, hbar], {nuk, lamk, Dk}, {fD -> Dk, fnu -> nuk},
   Take[kpzloTags, 1], kpzloRegs], "ExternalMomentum.ward"];
expectFail["no noise function: the generated Gamma2(hbar, hbar) is outside the ansatz",
  loEquationsOf[emTheory, hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]), {nuk, lamk, Dk}, {fnu -> nuk}, Take[kpzloTags, 2], kpzloRegs],
  "ExternalMomentum.zeroEntries"];
test["no noise function + Neglected[{hbar, hbar}]: built, the dropped entry is reported",
  MemberQ[quiet[loEquationsOf[emTheory, hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]), {nuk, lamk, Dk}, {fnu -> nuk},
       Append[Take[kpzloTags, 2], Neglected[{hbar, hbar}, "Reason" -> "noise not followed"]], kpzloRegs]]["Warnings"][[All, "Check"]], "ExternalMomentum.neglected"], True];
emTheory2D = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> FixedDimension[2]]];
expectFail["a jac vertex (d = 2): the loop depends on wedge(p, q)",
  loEquationsOf[emTheory2D, kpzloDensity + c2 hbar jac[h, lap[h]], {nuk, lamk, Dk, c2}, {fD -> Dk, fnu -> nuk}, kpzloTags, kpzloRegs], "ExternalMomentum.isotropicAngles"];
expectFail["a DE loop stage is refused", BuildLOEquations[kpzLoops, kpzSymmetries], "ExternalMomentum.pipeline"];

(* the DE-only stage functions refuse LO stages *)
expectFail["ProjectOnMonomials refuses LO stages", ProjectOnMonomials[kpzloLoops, kpzloSymmetries], "Projection.pipeline"];
expectFail["CheckClosure refuses an LO symmetries stage", CheckClosure[kpzIntegrated, kpzloSymmetries], "Closure.pipeline"];
expectFail["AssembleEquations refuses an LO stage", AssembleEquations[kpzloEquations], "Assemble.pipeline"];
expectFail["ExportJulia refuses an LO stage", ExportJulia[kpzloAdim, "Dimension" -> 2, "File" -> FileNameJoin[{$TemporaryDirectory, "lo.jl"}]], "JuliaExport.pipeline"];
