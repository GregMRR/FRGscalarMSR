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
