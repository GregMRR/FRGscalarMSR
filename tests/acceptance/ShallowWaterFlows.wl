(* ::Package:: *)

(* Milestone 4 acceptance: the raw RHS of all 9 DE4plus equations equal tests/reference/flowsFAST/*.wdx exactly
   (random rationals for every symbol and every shape atom), first with the abstract d_s R entries mapped to the script's
   placeholders (E10), then with the derived d_s R substituted on both sides. *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

(* {file, pair, momentum order, script LHS factor}: flowsFAST file = raw/factor; our value = coefficient = raw/(2 Dim)^n *)
swJobs = {{"gEQ", {H, Psibar}, 0, 1}, {"DphiaDphiEQ", {Psibar, Psibar}, 0, -2}, {"AnEQ", {H, Hbar}, 1, 2 Dim},
   {"DhaDhEQ", {Hbar, Hbar}, 1, -4 Dim}, {"KEQ", {Hbar, Psibar}, 1, 2 Dim}, {"nanEQ", {Psi, Psibar}, 1, 2 Dim},
   {"HEQ", {Hbar, Psi}, 1, -2 Dim}, {"g2EQ", {H, Psibar}, 1, 2 Dim}, {"DphiD2phiEQ", {Psibar, Psibar}, 1, -4 Dim}};
e10 = {dsR[H, Hbar] -> dsNh, dsR[Hbar, Hbar] -> dsMh, dsR[Psi, Psibar] -> dsNphi, dsR[Psibar, Psibar] -> dsMphi};
sq = LoopMomentumSquared;
(* explicit form: Delta coefficients -> atoms -> couplings and shapes, q^2 -> q^2 *)
explicitFlow[e_] := e /. {DeltaCoefficient[m_][x_] :> swaPropagator["DeltaCoefficients"][x][[m + 1]],
      Derivative[j_][DeltaCoefficient[m_]][x_] :> (D[swaPropagator["DeltaCoefficients"][xs][[m + 1]], {xs, j}] /. xs -> x)} /.
    {(a : GammaAtom[__])[x_] :> swaPropagator["AtomValues"][a][x], Derivative[j_][a : GammaAtom[__]][x_] :> Derivative[j][swaPropagator["AtomValues"][a]][x]} /.
   sq -> q^2;
flowOf[pair_, n_] := SelectFirst[swaIntegrated["Flows"], #Pair === pair && #Frequency == 0 && #Momentum == n &]["Value"];
derivedScript = {dsNh -> dsR[H, Hbar], dsMh -> dsR[Hbar, Hbar], dsNphi -> dsR[Psi, Psibar], dsMphi -> dsR[Psibar, Psibar]} /. DsRRules[swaRegulator, q^2];

Do[Module[{name = j[[1]], ref = Import[FileNameJoin[{testRoot, "reference", "flowsFAST", j[[1]] <> ".wdx"}]], ours},
   ours = explicitFlow[flowOf[j[[2]], j[[3]]]] (2 Dim)^j[[3]]/j[[4]];
   test["flowsFAST " <> name <> " (abstract d_s R, E10 map)", ExactlyZeroQ[(ours /. e10) - ref], True];
   test["flowsFAST " <> name <> " (derived d_s R)", ExactlyZeroQ[(ours /. DsRRules[swaRegulator, q^2]) - (ref /. derivedScript)], True]],
  {j, swJobs}];

test["HEQ vanishes identically", ExactlyZeroQ[explicitFlow[flowOf[{Hbar, Psi}, 1]] /. e10], True];
test["9 declared equations", Count[swaIntegrated["Flows"], KeyValuePattern["Class" -> "Declared"]], 9];
test["Protected, DeclaredNumeric and Causality flows sample to exactly 0",
  Select[swaIntegrated["Flows"], MemberQ[{"Protected", "DeclaredNumeric", "Causality"}, #Class] && Lookup[#, "Samples"] =!= {0, 0, 0} &][[All, {"Pair", "Frequency", "Momentum"}]], {}];
sampleOf[pair_, m_] := First[SelectFirst[swaIntegrated["Flows"], #Pair === pair && #Frequency == m && #Momentum == 0 &]["Samples"]];
test["script section 6b: i w and w^2 in Gamma_{H Psibar}, w^2 in Gamma_{Psibar Psibar} are generated (nonzero)",
  FreeQ[{sampleOf[{H, Psibar}, 1], sampleOf[{H, Psibar}, 2], sampleOf[{Psibar, Psibar}, 2]}, 0, {1}], True];
