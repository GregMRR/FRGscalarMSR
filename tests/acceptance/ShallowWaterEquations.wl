(* ::Package:: *)

(* Milestone 5 acceptance: assembled equations = tests/reference/flowsFAST/*.wdx exactly (derived d_s R on both sides);
   compact LeafCounts comparable to the script's compact.wdx *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

swFiles = {{"gEQ", gk}, {"DphiaDphiEQ", Dphik aDphik}, {"AnEQ", nk Ak}, {"DhaDhEQ", Dhk aDhk}, {"KEQ", Kk}, {"nanEQ", nk ank},
   {"HEQ", Hk lk}, {"g2EQ", g2k}, {"DphiD2phiEQ", Dphik D2phik}};
derivedScriptQ = {dsNh -> dsR[H, Hbar], dsMh -> dsR[Hbar, Hbar], dsNphi -> dsR[Psi, Psibar], dsMphi -> dsR[Psibar, Psibar]} /. DsRRules[swaRegulator, q^2];
scriptCompact = Import[FileNameJoin[{testRoot, "reference", "flowsFAST", "compact.wdx"}]]["equations"];
Do[Module[{name = f[[1]], eq = EquationOf[swaEquations, f[[2]]], ref = Import[FileNameJoin[{testRoot, "reference", "flowsFAST", f[[1]] <> ".wdx"}]]},
   test["Equations " <> name <> " = flowsFAST (exact)", ExactlyZeroQ[eq["Integrand"] - (ref /. derivedScriptQ)], True];
   test["Equations " <> name <> " compact LeafCount <= 1.5 x script (" <> ToString[eq["LeafCount"]["Compact"]] <> " vs " <>
     ToString[LeafCount[scriptCompact[name]]] <> ")", eq["LeafCount"]["Compact"] <= 3/2 LeafCount[scriptCompact[name]] + 10, True]],
  {f, swFiles}];
test["HEQ identically zero", EquationOf[swaEquations, Hk lk]["IdenticallyZero"], True];
test["dsCoupling enter linearly",
  AllTrue[swaEquations["Equations"], Function[e, AllTrue[e["DsCouplings"], ExactlyZeroQ[D[e["Integrand"], {#, 2}], "Seeds" -> {1}] &]]], True];
test["Dim symbolic in every nonzero p^2 equation", AllTrue[Select[swaEquations["Equations"], #Momentum == 1 && ! #IdenticallyZero &], ! FreeQ[#Integrand, Dim] &], True];
test["shapes abstract in q^2/k^2", Union[Cases[swaEquations["Equations"][[All, "Integrand"]], (arnh | armh | arnphi | armphi)[x_] | Derivative[_][arnh | armh | arnphi | armphi][x_] :> x, Infinity]], {q^2/k^2}];
