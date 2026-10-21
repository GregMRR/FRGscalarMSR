(* ::Package:: *)

(* The target of the Adimensional work: the LPA ansatz of ~/sync_cnrs/Projects/shallow_water/ShallowWaterLPA.nb (mk Psibar dt[H],
   Dmk grad[Hbar].grad[Psibar], forcing Dphik aDphik Psibar^2), with the regulator prefactor of {Hbar, H} changed to nk (plan D18).
   Dphik appears only as Dphik aDphik: eta_phi is fixed at run time (MUST_FIX) and the equation is the flow of aDphik. *)

Get[FileNameJoin[{testRoot, "reference", "JuliaRoundTrip.wl"}]];

lpaTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
lpaAnsatz = quiet[DefineAnsatz[lpaTheory, "Density" ->
      Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]])
    + Psibar (dt[Psi] + mk dt[H] - nk lap[Psi] + gk H + lk/2 grad[Psi].grad[Psi])
    - Dhk grad[Hbar].grad[Hbar] - Dmk grad[Hbar].grad[Psibar] - Dphik aDphik Psibar^2,
    "Couplings" -> {nk, Ak, Hk, lk, gk, mk, Dhk, Dmk, Dphik, aDphik}]];
lpaRegulator = quiet[DefineRegulator[lpaAnsatz,
    {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2],
    "Shapes" -> {arnh, arnphi, armphi}]];
lpaSymmetries = quiet[DeclareSymmetries[lpaAnsatz,
    Protected[{Hbar, _}, "Momentum" -> 0, "Reason" -> "Hbar -> Hbar + c(t) (mass conservation)"],
    Protected[{Psi, _}, "Momentum" -> 0, "Reason" -> "Psi -> Psi + c(t)"],
    Neglected[{H, Psibar}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "subdominant to gk"],
    Neglected[{Psibar, Psibar}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "subdominant to Dphik"],
    Neglected[_, "Frequency" -> 1, "Momentum" -> 1, "Reason" -> "mixed i w q^2 / w^2 q^2 terms, dropped"],
    NonRenormalised[lk, "Reason" -> "Galilean invariance"],
    "TruncationOrders" -> <|"Frequency" -> 1, "Momentum" -> 1|>]];
lpaEquations = Module[{v = quiet[BuildVertices[lpaAnsatz]], p},
   p = quiet[BuildPropagator[v, lpaRegulator]];
   quiet[AssembleEquations[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[v, p, lpaRegulator], lpaSymmetries]], lpaSymmetries]]]];

(* the workflow of AdimensionalLPAsw.nb: dimensions from the ansatz embedded in the (saved) equations *)
lpaFile = FileNameJoin[{$roundTripDir, "lpaEquations.wdx"}];
quiet[SaveStage[lpaEquations, lpaFile]];
lpaLoaded = quiet[LoadStage[lpaFile]];
lpaDims = quiet[DeriveDimensions[AnsatzOf[lpaLoaded], "References" -> {nk -> "nu", Dhk -> "h", Dphik -> "phi"}, "Dimensionless" -> {Ak, aDphik}]];
test["LPA adim: hand dimensions (old notebook conventions)", quiet[CheckDimensions[lpaDims,
    <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, Hk -> k^(Dim/2) Dhk^(1/2) nk^(-1/2), lk -> k^(1 - Dim/2) Dphik^(-1/2) nk^(3/2),
      mk -> k^-1 Dphik^(1/2) Dhk^(-1/2), Dmk -> k^-1 Dphik^(1/2) Dhk^(1/2), "Frequency" -> k^2 nk|>]], {gk, Hk, lk, mk, Dmk, "Frequency"}];
lpaAdim = quiet[AdimensionalEquations[lpaLoaded, lpaDims, "State" -> {gk, lk, Hk, Ak, mk, Dmk, aDphik}]];
test["LPA adim: nk, Dhk running, Dphik fixed", lpaAdim["ReferenceKinds"], <|nk -> "Running", Dhk -> "Running", Dphik -> "Fixed"|>];
test["LPA adim: MUST_FIX = phi", lpaAdim["MustFix"], {"phi"}];
test["LPA adim: state names", lpaAdim["Couplings"], {agk, alk, aHk, Ak, amk, aDmk, aDphik}];
test["LPA adim: linear parts", Expand[lpaAdim["BetaL"]], Expand[{agk (-1 - etah/2 + etanu + etaphi/2), alk (-1 + Dim/2 + 3 etanu/2 - etaphi/2),
     aHk (-Dim/2 + etah/2 - etanu/2), Ak etanu, amk (1 - etah/2 + etaphi/2), aDmk (1 + etah/2 + etaphi/2), aDphik etaphi}]];
test["LPA adim: lk has no loop part (NonRenormalised)", lpaAdim["BetaNL"][[2]], 0];
test["LPA adim: placeholder eta_phi row", {lpaAdim["EtaMNL"][[3]], lpaAdim["EtaBNL"][[3]]}, {{0, 0, 0}, 0}];
juliaRoundTrip["LPA2D", lpaAdim, 2];
