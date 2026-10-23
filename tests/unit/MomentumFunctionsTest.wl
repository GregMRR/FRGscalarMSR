(* ::Package:: *)

(* Momentum-dependent couplings (design/MOMENTUMFUNCTIONS.md): every new check of the hint table, the generic-p chain (oracle at
   generic p, n0 = 0, 1, 2, mixed ansatz with a mass term), the p -> 0 regression against plain KPZ, Dimensions, Adimensional, Related,
   and the export refusal. *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];

mfDensity = hbar (dt[h] - nuk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2;
mfAnsatz[d_, cps_, fns_, o___] := DefineAnsatz[kpzloTheory, "Density" -> d, "Couplings" -> cps, "MomentumFunctions" -> fns, o];
mfTags = Sequence[Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift symmetry"], Neglected[_, "Reason" -> "beyond the ansatz"],
   NonRenormalised[lamk, "Reason" -> "Galilean invariance"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>];

(* ---- Ansatz (A2) ---- *)
test["Ansatz: pipeline set from the functions", kpzloAnsatz["Pipeline"], "MomentumFunctions"];
test["Ansatz: functions stored", kpzloAnsatz["MomentumFunctions"], {Fnu, FD}];
test["Ansatz: function-free stage has no MomentumFunctions key", KeyExistsQ[kpzAnsatz, "MomentumFunctions"], False];
expectFail["Ansatz: Fnu[0] aborts", mfAnsatz[hbar (dt[h] - nuk Fnu[0] lap[h]) - Dk hbar^2, {nuk, Dk}, {Fnu}], "Ansatz.functionBare"];
expectFail["Ansatz: a function in a cubic term aborts", mfAnsatz[hbar (dt[h] - nuk lap[h] - Fnu lamk grad[h].grad[h]) - Dk hbar^2, {nuk, lamk, Dk}, {Fnu}],
  "Ansatz.functionQuadraticOnly"];
expectFail["Ansatz: a function under an operator aborts", mfAnsatz[hbar dt[h] - nuk hbar lap[Fnu h] - Dk hbar^2, {nuk, Dk}, {Fnu}], "Ansatz.functionNotInOperator"];
expectFail["Ansatz: a function also in Couplings aborts", mfAnsatz[mfDensity, {nuk, lamk, Dk, Fnu}, {Fnu, FD}], "Ansatz.functionsDistinct"];
expectFail["Ansatz: a contradicting Pipeline aborts", mfAnsatz[mfDensity, {nuk, lamk, Dk}, {Fnu, FD}, "Pipeline" -> "DerivativeExpansion"], "Ansatz.pipeline"];
expectFail["Ansatz: an undeclared function aborts", mfAnsatz[mfDensity, {nuk, lamk, Dk}, {Fnu}], "Ansatz.knownSymbols"];
expectFail["Ansatz: a function named c (reserved) aborts", mfAnsatz[hbar (dt[h] - nuk c lap[h]) - Dk hbar^2, {nuk, Dk}, {c}], "Ansatz.functions.reserved"];

(* ---- Vertices (A3) ---- *)
test["Vertices: Gamma2_{h hbar} = I w + nuk Fnu(q^2) q^2", Expand[Gamma2Matrix[kpzloVertices, w, q][[1, 2]] - (I w + nuk Fnu[scal[q, q]] scal[q, q])], 0];
test["Vertices: Gamma2_{hbar hbar} = -2 Dk FD(q^2)", Gamma2Matrix[kpzloVertices, w, q][[2, 2]], -2 Dk FD[scal[q, q]]];
Scan[expectFail["Vertices: functional term " <> #[[1]] <> " aborts", BuildVertices[mfAnsatz[#[[2]], #[[3]], {Fnu}]], "Vertices.functionalTerm"] &,
  {{"Fnu^2", hbar (dt[h] - nuk Fnu^2 lap[h]) - Dk hbar^2, {nuk, Dk}},
   {"1/Fnu", hbar (dt[h] - nuk lap[h]/Fnu) - Dk hbar^2, {nuk, Dk}},
   {"nuk lamk Fnu (two couplings)", hbar (dt[h] - nuk lamk Fnu lap[h]) - Dk hbar^2, {nuk, lamk, Dk}},
   {"Fnu at two powers of p^2", hbar (dt[h] - nuk Fnu lap[h] + nuk Fnu lap[lap[h]]) - Dk hbar^2, {nuk, Dk}}}];

(* ---- Regulator (A5) ---- *)
expectFail["Regulator: a function in R aborts", DefineRegulator[kpzloAnsatz, {hbar, h} -> nuk Fnu q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2],
   "Shapes" -> {rn, rd}], "Regulator.noFunctions"];

(* ---- Symmetries (A4, B) ---- *)
mfFunctional = Select[kpzloSymmetries["Classification"], #Class === "Functional" &];
test["Symmetries: two Functional rows (nuk Fnu at n0 = 1, -2 Dk FD at n0 = 0)", {#Pair, #Momentum, #GammaK, #Function} & /@ mfFunctional,
  {{{h, hbar}, 1, nuk Fnu, Fnu}, {{hbar, hbar}, 0, -2 Dk FD, FD}}];
test["Symmetries: the p^0 order of the Fnu slot stays an ordinary Protected row",
  Lookup[SelectFirst[kpzloSymmetries["Classification"], #Pair === {h, hbar} && #Frequency == 0 && #Momentum === 0 &], "Class"], "Protected"];
test["Symmetries: sibling I w slot reported (f^lambda = 1)", MemberQ[kpzloSymmetries["Warnings"][[All, "Check"]], "Symmetries.siblingAtZeroMomentum"], True];
expectFail["Symmetries: a constant at the function's order aborts", DeclareSymmetries[mfAnsatz[hbar (dt[h] - nuk lap[h] - nuk Fnu lap[h]) - Dk hbar^2,
    {nuk, Dk}, {Fnu}], Neglected[_, "Reason" -> "x"], "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>], "Symmetries.functionalSlot"];
expectFail["Symmetries: two functions in one slot abort", DeclareSymmetries[mfAnsatz[hbar (dt[h] - nuk Fnu lap[h] + nuk Gnu lap[lap[h]]) - Dk hbar^2,
    {nuk, Dk}, {Fnu, Gnu}], Neglected[_, "Reason" -> "x"], "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>], "Symmetries.functionalSlot"];
expectFail["Symmetries: a function whose term vanishes in Gamma2 aborts", DeclareSymmetries[mfAnsatz[hbar (dt[h] - nuk lap[h]) - Dk hbar^2 + Fnu div[hbar grad[hbar]],
    {nuk, Dk}, {Fnu}], Neglected[_, "Reason" -> "x"], "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>], "Symmetries.functionSlot"];
expectFail["Symmetries: a function in two slots aborts", DeclareSymmetries[mfAnsatz[hbar (dt[h] - nuk Fnu lap[h]) - Dk Fnu hbar^2,
    {nuk, Dk}, {Fnu}], Neglected[_, "Reason" -> "x"], "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>], "Symmetries.functionSlot"];
expectFail["Symmetries: Neglected at order n0 of a functional slot aborts", DeclareSymmetries[kpzloAnsatz,
   Neglected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "x"], mfTags], "Symmetries.functionalTag"];
expectFail["Symmetries: AllMomenta on a functional slot aborts", DeclareSymmetries[kpzloAnsatz,
   Neglected[{hbar, hbar}, "Frequency" -> 0, "Momentum" -> "AllMomenta", "Reason" -> "x"], mfTags], "Symmetries.contradiction"];
expectFail["Symmetries: NonRenormalised[function] aborts", DeclareSymmetries[kpzloAnsatz, NonRenormalised[Fnu, "Reason" -> "x"], mfTags],
  "Symmetries.nonRenormalisedFunction"];
expectFail["Symmetries: Momentum -> \"All\" (a string) aborts", DeclareSymmetries[kpzloAnsatz, Neglected[_, "Momentum" -> "All", "Reason" -> "x"], mfTags],
  "Symmetries.tagMomentum"];
mfAllSym = quiet[DeclareSymmetries[kpzloAnsatz, Neglected[{h, hbar}, "Frequency" -> 1, "Momentum" -> "AllMomenta", "Reason" -> "LO: f^lambda = 1"], mfTags]];
test["Symmetries: an AllMomenta tag gives one generic-p row", {#Pair, #Frequency, #Class} & /@ Select[mfAllSym["Classification"], #Momentum === "AllMomenta" &],
  {{{h, hbar}, 1, "Neglected"}}];
test["Symmetries: the tagged sibling is no longer reported", Select[mfAllSym["Warnings"], #Check === "Symmetries.siblingAtZeroMomentum" && StringContainsQ[#Detail, "(I w)^1"] &], {}];

(* ---- Projection, FrequencyIntegration (C) ---- *)
test["Projection: functional integrands in invariants only (no vectors)",
  FreeQ[Select[kpzloProjection["Integrands"], #Class === "Functional" &][[All, {"Integrand", "AtZero"}]], _scal], True];
mfFam[n_, ex_, deg_, cv_] := FrequencyIntegralFamilies[n, ex, deg] /. Flatten[Table[FrequencyIntegralC[i, m] -> cv[[i, m + 1]], {i, Length[cv]}, {m, 0, 2}]];
mfNum[n_, ex_, deg_, cv_] := Module[{dl},
   dl[i_, x_] := cv[[i, 1]] + cv[[i, 2]] I x + If[deg == 2, cv[[i, 3]] (I x)^2, 0];
   NIntegrate[x^n/Product[dl[i, x]^ex[[i, 1]] dl[i, -x]^ex[[i, 2]], {i, Length[ex]}], {x, -Infinity, Infinity}, WorkingPrecision -> 30, PrecisionGoal -> 15]];
Scan[test["FrequencyIntegration: two families = NIntegrate, " <> #[[1]], Abs[N[mfFam @@ Rest[#], 30] - mfNum @@ Rest[#]] < 10^-12, True] &,
  {{"degree 1", 0, {{1, 1}, {1, 1}}, 1, {{3/2, 1, 0}, {5/7, 1, 0}}},
   {"degree 1, higher powers", 1, {{2, 1}, {1, 2}}, 1, {{3/2, 2/3, 0}, {5/7, 1, 0}}},
   {"degree 2 overdamped", 2, {{2, 1}, {1, 1}}, 2, {{1/5, 3, 1}, {1/3, 4, 1}}},
   {"degree 2 underdamped", 0, {{1, 1}, {1, 2}}, 2, {{3, 1, 1}, {4, 3/2, 2}}}}];
test["FrequencyIntegration: coincident families = one family with summed powers (finite, E3)",
  Together[(FrequencyIntegralFamilies[1, {{2, 1}, {1, 2}}, 2] /. FrequencyIntegralC[2, m_] :> FrequencyIntegralC[1, m] /. FrequencyIntegralC[1, m_] :> FrequencyIntegralC[m]) -
    FrequencyIntegral[1, 3, 3, 2]], 0];

(* ---- Closure (B2) ---- *)
expectFail["Closure: Protected at order n0 of the Fnu slot (nu(0) does not flow) is false in KPZ and aborts",
  Module[{s = DeclareSymmetries[kpzloAnsatz, Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "false"], mfTags]},
   CheckClosure[IntegrateFrequencies[ProjectOnMonomials[kpzloLoops, s]], s]], "Closure.protectedAtZero"];

(* ---- Assemble, p -> 0 regression (E3 of the plan) ---- *)
test["Assemble: EquationOf[eqs, Fnu] is the functional equation", EquationOf[kpzloEquations, Fnu]["Flowing"], Fnu nuk];
test["Assemble: functional integrands in (q, p, c)", Sort[Intersection[Variables[EquationOf[kpzloEquations, FD]["Integrand"] /. (_rn | _rd | _Fnu | _FD |
       Derivative[__][_][_]) -> 1], {q, p, c}]], Sort[{q, p, c}]];
test["p -> 0 regression: the exact p = 0 row of Fnu (F = 1) = plain KPZ nu equation",
  ExactlyZeroQ[(EquationOf[kpzloEquations, Fnu]["AtZero"]["Integrand"] /. {Fnu -> (1 &), FD -> (1 &)}) - EquationOf[kpzEquations, nuk]["Integrand"]], True];
test["p -> 0 regression: the exact p = 0 row of FD (F = 1) = plain KPZ D equation",
  ExactlyZeroQ[(EquationOf[kpzloEquations, FD]["AtZero"]["Integrand"] /. {Fnu -> (1 &), FD -> (1 &)}) - EquationOf[kpzEquations, Dk]["Integrand"]], True];

(* ---- Oracle at generic p (E2): KPZ LO (n0 = 1, 0), a mixed ansatz with a mass term (n0 = 1 over a declared p^0), a p^4 function (n0 = 2, E15) ---- *)
test["Oracle: KPZ LO functional equations at generic p and p = 0, d = 2, 3", AssociationQ[quiet[RunOracle[kpzloEquations]]], True];
mfPipeline[an_, sym_, regRules_, shapes_] := Module[{v = BuildVertices[an], reg, prop},
   reg = DefineRegulator[an, Sequence @@ regRules, "Shapes" -> shapes]; prop = BuildPropagator[v, reg];
   AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[v, prop, reg], sym]], sym]]];
mfMass = quiet[Module[{an = mfAnsatz[hbar (dt[h] + mk h - nuk Fnu lap[h] - lamk/2 grad[h].grad[h]) - Dk FD hbar^2, {mk, nuk, lamk, Dk}, {Fnu, FD}]},
    mfPipeline[an, DeclareSymmetries[an, Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean"], Neglected[_, "Reason" -> "x"],
      NonRenormalised[lamk, "Reason" -> "Galilean"], "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>],
     {{hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]}, {rn, rd}]]];
test["Oracle: mixed ansatz, mass mk below Fnu (n0 = 1), d = 2", AssociationQ[quiet[RunOracle[mfMass, "Dimensions" -> {2}]]], True];
mfP4 = quiet[Module[{an = mfAnsatz[hbar (dt[h] - nuk lap[h] + nuk Ft lap[lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar^2, {nuk, lamk, Dk}, {Ft}]},
    mfPipeline[an, DeclareSymmetries[an, mfTags], {{hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]}, {rn, rd}]]];
test["Projection: a p^4 function over a declared p^2 (n0 = 2) is regular", EquationOf[mfP4, Ft]["Momentum"], 2];
test["Oracle: p^4 function (n0 = 2, midpoint subtraction), d = 2", AssociationQ[quiet[RunOracle[mfP4, "Dimensions" -> {2}]]], True];

(* ---- Dimensions (A3, D6) ---- *)
mfDims = quiet[DeriveDimensions[kpzloAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]];
test["Dimensions: a function times a reference is dimensionless", {DimensionOf[mfDims, Fnu], DimensionOf[mfDims, FD]}, {1, 1}];
test["Dimensions: dimensionless names aFnu, aFD", {mfDims["Names"][Fnu], mfDims["Names"][FD]}, {aFnu, aFD}];
expectFail["Dimensions: a function as reference aborts", DeriveDimensions[kpzloAnsatz, "References" -> {Fnu -> "nu", Dk -> "D"}], "Dimensions.functionReference"];
expectFail["Dimensions: a function times a non-reference aborts", DeriveDimensions[kpzloAnsatz, "References" -> {lamk -> "l", Dk -> "D"}],
  "Dimensions.functionPrefactor"];
expectFail["Dimensions: Mode Fixed with Equation aborts", DeriveDimensions[kpzloAnsatz, "References" -> {nuk -> {"nu", "Mode" -> "Fixed", "Equation" -> nuk Fnu},
    Dk -> "D"}], "Dimensions.equationModeFixed"];
test["Dimensions: Equation makes the reference Running", quiet[DeriveDimensions[kpzloAnsatz, "References" -> {nuk -> {"nu", "Equation" -> nuk Fnu}, Dk -> "D"}]][
   "ReferenceModes"][nuk], "Running"];
mfFlam = quiet[mfAnsatz[hbar (Flam dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, {nuk, lamk, Dk}, {Flam}]];
expectFail["Dimensions: a bare function with an undetermined dimension aborts (hint: Dimensionless)", DeriveDimensions[mfFlam, "References" -> {nuk -> "nu", Dk -> "D"}],
  "Dimensions.freedoms"];
test["Dimensions: a bare function declared Dimensionless is accepted (E17)", DimensionOf[quiet[DeriveDimensions[mfFlam, "References" -> {nuk -> "nu", Dk -> "D"},
     "Dimensionless" -> {Flam}]], Flam], 1];

(* ---- Adimensional (D3-D7), JuliaExport (E1) ---- *)
mfAdim = quiet[AdimensionalEquations[kpzloEquations, mfDims]];
test["Adimensional: functions in the state, dimensionless names", {mfAdim["State"], mfAdim["MomentumFunctions"]}, {{lamk, Fnu, FD}, {aFnu, aFD}}];
test["Adimensional: linear part of aFnu = etanu aFnu + 2 ap^2 aFnu'", Expand[mfAdim["BetaL"][[2]] - (etanu aFnu[ap^2] + 2 ap^2 Derivative[1][aFnu][ap^2])], 0];
test["Adimensional: references only times functions are Fixed by default", mfAdim["MustFix"], {"nu", "D"}];
mfDimsRun = quiet[DeriveDimensions[kpzloAnsatz, "References" -> {nuk -> {"nu", "Equation" -> nuk Fnu}, Dk -> {"D", "Equation" -> Dk FD}}]];
mfAdimRun = quiet[AdimensionalEquations[kpzloEquations, mfDimsRun]];
kpzAdimRef = quiet[AdimensionalEquations[kpzEquations, quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]]]];
test["Adimensional: eta rows normalised on Fnu, FD (F = 1) = plain KPZ eta rows (C4)",
  ExactlyZeroQ[{(mfAdimRun["EtaBNLExplicit"] - kpzAdimRef["EtaBNLExplicit"]), (mfAdimRun["EtaMNLExplicit"] - kpzAdimRef["EtaMNLExplicit"])} /. {aFnu -> (1 &), aFD -> (1 &)}], True];
mfRel = quiet[AdimensionalEquations[kpzloEquations, mfDimsRun, "Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "time reversal (d = 1)"]}]];
test["Related: time reversal at d = 1 holds; FD eliminated, d = 1, eta_D = eta_nu imposed", {mfRel["State"], mfRel["Dimension"], mfRel["EtaML"]},
  {{lamk, Fnu}, 1, {{1, 0}, {-1, 1}}}];
expectFail["Related: time reversal at d = 2 is rejected (E34)", AdimensionalEquations[kpzloEquations, mfDimsRun,
   "Related" -> {Related[FD -> Fnu, "AtDimension" -> 2, "Reason" -> "false"]}], "Adimensional.related"];
expectFail["Related: relating a function to a coupling aborts", AdimensionalEquations[kpzloEquations, mfDimsRun,
   "Related" -> {Related[FD -> nuk, "Reason" -> "x"]}], "Adimensional.relatedFunctions"];
expectFail["Related: an eliminated function in State aborts", AdimensionalEquations[kpzloEquations, mfDimsRun, "State" -> {lamk, Fnu, FD},
   "Related" -> {Related[FD -> Fnu, "AtDimension" -> 1, "Reason" -> "x"]}], "Adimensional.relatedState"];
expectFail["Related: a cycle aborts", AdimensionalEquations[kpzloEquations, mfDimsRun,
   "Related" -> {Related[FD -> Fnu, "Reason" -> "x"], Related[Fnu -> FD, "Reason" -> "y"]}], "Adimensional.relatedChain"];
expectFail["JuliaExport: functional equations are not exportable", ExportJulia[mfAdim, "Dimension" -> 1, "File" -> FileNameJoin[{$TemporaryDirectory, "mf.jl"}]],
  "JuliaExport.functionalCouplings"];
