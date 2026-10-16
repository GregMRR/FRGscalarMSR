(* ::Package:: *)

(* LPA acceptance: the LPA ansatz (DE4plus without the Kk, g2k, D2phik terms) derived from scratch by the library equals the
   equations of the old notebooks, exactly (random rationals for every symbol and shape atom):
   (1) tests/reference/legacyLPA/flowsDE4plus_agreed /. DE4toLPA  -- input of LPA_KPZ_thermal_forcing/AdimensionalFlowEQGenericReg.nb,
       already normalised flows, aDhk = ank = aDphik = 1 implicit;
   (2) tests/reference/legacyLPA/flowsLPAold  -- the oldest LPA set (lambda = 1, A = 0, no H-noise regulator), raw RHS with the
       factors of shallow_water_DE4plus.nb section 14.
   Symbol maps (OldtoNew, DimensionlessRegulators) copied from shallow_water_DE4plus.nb; see tests/reference/legacyLPA/README.md. *)

lpaTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
lpaAnsatz = quiet[DefineAnsatz[lpaTheory, "Density" ->
     Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) + Psibar (dt[Psi] - nk ank lap[Psi] + gk H + lk/2 grad[Psi].grad[Psi]) -
     Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2, "Couplings" -> {nk, Ak, Hk, lk, ank, gk, Dhk, aDhk, Dphik, aDphik}]];
lpaVertices = quiet[BuildVertices[lpaAnsatz]];
lpaRegulator = quiet[DefineRegulator[lpaAnsatz, {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2],
    {Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2], "Shapes" -> {arnh, arnphi, armh, armphi}]];
lpaSymmetries = quiet[DeclareSymmetries[lpaAnsatz,
    Protected[{Hbar, _}, "Momentum" -> 0, "Reason" -> "Hbar -> Hbar + c(t) (mass conservation)"],
    Protected[{Psi, _}, "Momentum" -> 0, "Reason" -> "Psi -> Psi + c(t)"],
    Neglected[{H, Psibar}, "Frequency" -> {1, 2}, "Momentum" -> 0, "Reason" -> "truncation, z = 1 power counting"],
    Neglected[{Psibar, Psibar}, "Frequency" -> 2, "Momentum" -> 0, "Reason" -> "coloured forcing, dropped"],
    Neglected[_, "Frequency" -> {1, 2}, "Momentum" -> 1, "Reason" -> "mixed i w q^2 / w^2 q^2 terms, dropped"],
    Neglected[{H, Psibar} | {Hbar, Psibar} | {Psibar, Psibar}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "LPA: g2, K, D2phi dropped"],
    NonRenormalised[lk, "Reason" -> "Galilean invariance"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
lpaPropagator = quiet[BuildPropagator[lpaVertices, lpaRegulator]];
lpaClosure = quiet[CheckClosure[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[lpaVertices, lpaPropagator, lpaRegulator], lpaSymmetries]], lpaSymmetries]];
lpaEquations = quiet[AssembleEquations[lpaClosure]];
test["LPA: 6 equations (A, g, Dh, H, N, Dphi)", Sort[lpaEquations["Equations"][[All, "Flowing"]]], Sort[{nk Ak, gk, Dhk aDhk, Hk lk, nk ank, Dphik aDphik}]];
test["LPA: g2, K, D2phi flows reported as Neglected",
  Sort[#Pair & /@ Select[lpaClosure["Report"], #Class === "Neglected" && #Frequency == 0 && #Momentum == 1 &]], Sort[{{H, Psibar}, {Hbar, Psibar}, {Psibar, Psibar}}]];

(* ours: explicit, d_s R abstract with the script's names (E10), bookkeeping couplings = 1 *)
lpaE10 = {dsR[H, Hbar] -> dsNh, dsR[Hbar, Hbar] -> dsMh, dsR[Psi, Psibar] -> dsNphi, dsR[Psibar, Psibar] -> dsMphi};
lpaOurs[c_] := (EquationOf[lpaEquations, c]["Compact"] /. Normal[lpaEquations["Atoms"]] /. lpaE10) /. {aDhk -> 1, ank -> 1, aDphik -> 1};
(* old notation -> new: OldtoNew + DimensionlessRegulators of shallow_water_DE4plus.nb, H_0 = H_k, dsN\[Phi] = dsNphi *)
lpaOldToNew = {Subscript[N, k] -> nk, Subscript[\[Lambda], k] -> lk, Subscript[K, k] -> Kk, Subscript[A, k] -> Ak, Subscript[D, hk] -> Dhk,
   Subscript[D, \[Phi]k] -> Dphik, Subscript[D, 2 \[Phi]k] -> D2phik, Subscript[g, k] -> gk, Subscript[g, 2 k] -> g2k, Subscript[H, k] -> Hk, Subscript[H, 0] -> Hk,
   Subscript[rn, h][q^2] -> arnh[q^2/k^2], Subscript[rn, \[Phi]][q^2] -> arnphi[q^2/k^2], Subscript[rm, h][q^2] -> armh[q^2/k^2], Subscript[rm, \[Phi]][q^2] -> armphi[q^2/k^2],
   Derivative[1][Subscript[rn, h]][q^2] -> Derivative[1][arnh][q^2/k^2]/k^2, Derivative[1][Subscript[rn, \[Phi]]][q^2] -> Derivative[1][arnphi][q^2/k^2]/k^2,
   Derivative[1][Subscript[rm, h]][q^2] -> Derivative[1][armh][q^2/k^2]/k^2, Derivative[1][Subscript[rm, \[Phi]]][q^2] -> Derivative[1][armphi][q^2/k^2]/k^2,
   dsN\[Phi] -> dsNphi, dsM\[Phi] -> dsMphi};
lpaRef[set_, f_] := Import[FileNameJoin[{testRoot, "reference", "legacyLPA", set, f <> ".wdx"}]] /. lpaOldToNew;
lpaNew = {nk, lk, Kk, Ak, Dhk, Dphik, D2phik, gk, g2k, Hk, k, q, Dim, dsNh, dsMh, dsNphi, dsMphi, arnh, arnphi, armh, armphi};

(* (1) flowsDE4plus_agreed /. DE4toLPA = {Kk, D2phik, g2k -> 0} (DE4toLPA3 of AdimensionalFlowEQGenericReg.nb) *)
Do[Module[{ref = lpaRef["flowsDE4plus_agreed", j[[1]]] /. {Kk -> 0, D2phik -> 0, g2k -> 0}},
   test["LPA vs flowsDE4plus_agreed " <> j[[1]] <> ": old notation fully mapped", Complement[Union[Cases[ref, _Symbol | _Subscript, Infinity]], lpaNew], {}];
   test["LPA vs flowsDE4plus_agreed " <> j[[1]] <> " (exact)", ExactlyZeroQ[lpaOurs[j[[2]]] - ref], True]],
  {j, {{"gEQ", gk}, {"DphiEQ", Dphik aDphik}, {"NEQ", nk ank}, {"DhEQ", Dhk aDhk}, {"ANEQ", nk Ak}}}];

(* (2) flowsLPAold: lambda = 1, A = 0, no H-noise regulator; raw RHS times the LHS factors of shallow_water_DE4plus.nb section 14 *)
Do[Module[{ref = lpaRef["flowsLPAold", j[[1]]]},
   test["LPA vs flowsLPAold " <> j[[1]] <> " (exact, factor " <> ToString[j[[3]], InputForm] <> ")",
    ExactlyZeroQ[(lpaOurs[j[[2]]] /. {Ak -> 0, lk -> 1, armh -> (0 &), dsMh -> 0}) - j[[3]] ref], True]],
  {j, {{"gEQ", gk, 1}, {"DphiEQ", Dphik aDphik, -1/2}, {"NEQ", nk ank, 1/(2 Dim)}, {"DhEQ", Dhk aDhk, -1/(4 Dim)}}}];
