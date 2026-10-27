(* ::Package:: *)

(* Adimensional: the reference rules and the aborts, on KPZ (cheap); the physics is in tests/acceptance/*Adimensional.wl *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
adKPZDims = quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]];
adKPZ = quiet[AdimensionalEquations[kpzEquations, adKPZDims]];
test["KPZ: both references running, nothing to fix", {Values[adKPZ["ReferenceKinds"]], adKPZ["MustFix"]}, {{"Running", "Running"}, {}}];
test["KPZ: default state = non-references", adKPZ["State"], {lamk}];
test["KPZ: eta matrix l part = identity, no l part of b", {adKPZ["EtaML"], adKPZ["EtaBL"]}, {IdentityMatrix[2], {0, 0}}];
test["KPZ: BetaL free of aq", FreeQ[adKPZ["BetaL"], aq], True];

expectFail["a NonRenormalised reference", AdimensionalEquations[kpzEquations,
   quiet[DeriveDimensions[kpzAnsatz, "References" -> {lamk -> "l", Dk -> "D"}]]], "Adimensional.referencesFlow"];
expectFail["Mode -> Fixed on a reference that appears alone", AdimensionalEquations[kpzEquations,
   quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> {"nu", "Mode" -> "Fixed"}, Dk -> "D"}]]], "Adimensional.referenceMode"];
expectFail["a reference in the state", AdimensionalEquations[kpzEquations, adKPZDims, "State" -> {lamk, nuk}], "Adimensional.state"];
expectFail["dimensions of another ansatz", AdimensionalEquations[kpzEquations,
   quiet[DeriveDimensions[quiet[DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - 2 Dk hbar^2,
       "Couplings" -> {nuk, lamk, Dk}]], "References" -> {nuk -> "nu", Dk -> "D"}]]], "Adimensional.sameAnsatz"];

(* a non-reference coupling in a regulator prefactor (plan D4) *)
adAmpAnsatz = quiet[DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk Ak lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, Ak, lamk, Dk}]];
adAmpSym = quiet[DeclareSymmetries[adAmpAnsatz,
    Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
    Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
    Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
    "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
adAmpEquations[reg_] := Module[{v = quiet[BuildVertices[adAmpAnsatz]], r = quiet[DefineRegulator[adAmpAnsatz, Sequence @@ reg, "Shapes" -> {rn, rd}]], p},
   p = quiet[BuildPropagator[v, r]];
   quiet[AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[v, p, r], adAmpSym]], adAmpSym]]]];
adAmpDims = quiet[DeriveDimensions[adAmpAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}, "Dimensionless" -> {Ak}]];
expectFail["regulator prefactor nuk Ak (not a reference)",
  AdimensionalEquations[adAmpEquations[{{hbar, h} -> nuk Ak q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]}], adAmpDims], "Adimensional.regulatorReferences"];

(* nuk appears only as nuk Ak: a fixed reference, the equation becomes the flow of Ak, eta_nu has a placeholder row *)
adAmp = quiet[AdimensionalEquations[adAmpEquations[{{hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]}], adAmpDims]];
test["amplitude: nuk fixed, Dk running", adAmp["ReferenceKinds"], <|nuk -> "Fixed", Dk -> "Running"|>];
test["amplitude: MUST_FIX = nu", adAmp["MustFix"], {"nu"}];
test["amplitude: beta_A = A eta_nu + loop (linear part A eta_nu)", Expand[adAmp["BetaL"][[First[FirstPosition[adAmp["State"], Ak]]]]], Ak etanu];
test["amplitude: placeholder row for eta_nu", {adAmp["EtaMNL"][[1]], adAmp["EtaBNL"][[1]]}, {{0, 0}, 0}];
test["amplitude: Mode -> Running restores an eta row", quiet[AdimensionalEquations[adAmpEquations[{{hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]}],
     quiet[DeriveDimensions[adAmpAnsatz, "References" -> {nuk -> {"nu", "Mode" -> "Running"}, Dk -> "D"}, "Dimensionless" -> {Ak}]]]]["MustFix"], {}];

(* causal-pole conditions (design/CAUSALITY.md): minimal dimensionless form, its check and its abort *)
test["KPZ causal poles: compact = the atom", adKPZ["CausalPoles"]["Compact"], {AdimAtom[h, hbar, 0, 0]}];
test["KPZ causal poles: explicit 1 + rn(aq^2)", adKPZ["CausalPoles"]["Explicit"], {1 + rn[aq^2]}];
test["KPZ causal poles: stripped nuk k^2 aq^2", adKPZ["CausalPoles"]["Stripped"], {aq^2 k^2 nuk}];
Block[{FRGscalarMSR`Private`causalPolesStrip = Function[{r, x}, {1, 2 Expand[r]}]},
  expectFail["mutated causal-pole strip: wrong factor", AdimensionalEquations[kpzEquations, adKPZDims], "Adimensional.causalPolesMinimal"]];
Block[{FRGscalarMSR`Private`causalPolesStrip = Function[{r, x}, {-1, -Expand[r]}]},
  expectFail["mutated causal-pole strip: negative content", AdimensionalEquations[kpzEquations, adKPZDims], "Adimensional.causalPolesMinimal"]];
adAcausal = With[{q = kpzEquations["LoopMomentum"]}, FRGscalarMSR`Internal`stampStage[KeyDrop[Append[kpzEquations,
      "ValidityConditions" -> <|"Atoms" -> {-nuk q^2 > 0}, "Explicit" -> {-nuk q^2 > 0}|>], "Fingerprint"]]];
expectFail["a causal-pole condition that reduces to a negative constant", AdimensionalEquations[adAcausal, adKPZDims], "Adimensional.causalPoles"];

(* LO branch (design/LOPIPELINE.md step 7); the comparison with BMWflow's equations is in tests/acceptance/KPZLO.wl *)
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
test["LO: functions, etas in FUNCTIONS order, coupling", {kpzloAdim["Pipeline"], kpzloAdim["Functions"], kpzloAdim["Etas"], kpzloAdim["Couplings"], kpzloAdim["MustFix"]},
  {"LO", {fD, fnu}, {etaD, etanu}, {alamk}, {}}];
test["LO: linear part at fixed ap = p/k", Expand[kpzloAdim["Linear"]], Expand[{etaD fD[ap^2] + 2 ap^2 fD'[ap^2], etanu fnu[ap^2] + 2 ap^2 fnu'[ap^2]}]];
test["LO: beta of alamk = DE KPZ linear part", Expand[kpzloAdim["BetaL"] - adKPZ["BetaL"]], {0}];
test["LO: no unit in the loops and eta rows; f and r only at aq^2 and Qhat^2",
  {FreeQ[kpzloAdim /@ {"Loop", "EtaMatrix", "EtaVector"}, k | nuk | Dk | q | p], Union[Cases[kpzloAdim["Loop"], (fD | fnu | rn | rd)[x_] :> x, {0, Infinity}]]},
  {True, Sort[{aq^2, ap^2 + aq^2 + 2 ap aq CosTheta}]}];
test["LO: eta system linear, rows of the eta integrands", {Dimensions[kpzloAdim["EtaMatrix"]], FreeQ[kpzloAdim /@ {"EtaMatrix", "EtaVector"}, etaD | etanu]}, {{2, 2}, True}];
test["LO: measures", {Together[kpzloAdim["Measure"]["Eta"] - LoopMeasure[Dim, aq]], FreeQ[kpzloAdim["Measure"]["Loop"], Sin]}, {0, False}];
adLOnoNR = quiet[loEquationsOf[kpzloTheory, kpzloDensity, {nuk, lamk, Dk}, {fD -> Dk, fnu -> nuk}, Take[kpzloTags, 2], kpzloRegs]];
expectFail["LO: a constant coupling that is not NonRenormalised gets no flow",
  AdimensionalEquations[adLOnoNR, quiet[DeriveDimensions[AnsatzOf[adLOnoNR], "References" -> {nuk -> "nu", Dk -> "D"}]]], "Adimensional.noBeta"];
