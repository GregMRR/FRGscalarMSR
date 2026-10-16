(* ::Package:: *)

(* DE4plus acceptance against the old notebooks, exactly (random rationals for every symbol and shape atom):
   (1) tests/reference/legacy/flowsDE4plus_agreed (old notation, normalised flows, aDhk = ank = aDphik = 1 implicit);
   (2) tests/reference/legacy/flowsNEW (shallow_water_DE4plus.nb exports, same normalisation as flowsFAST).
   Not usable as references, see tests/reference/legacy/README.md: agreed HEQ (placeholder symbol HkEQ), agreed g2EQ and DphiD2phiEQ
   (unsubstituted f''; they disagree with this library, flowsNEW and the oracle), flowsNEW DphiD2phiEQ (corrupt, not copied). *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "legacy", "OldNotation.wl"}]];

Do[Module[{ref = legacyRef["flowsDE4plus_agreed", j[[1]]]},
   test["DE4 vs flowsDE4plus_agreed " <> j[[1]] <> ": old notation fully mapped", Complement[Union[Cases[ref, _Symbol | _Subscript, Infinity]], newNotationSymbols], {}];
   test["DE4 vs flowsDE4plus_agreed " <> j[[1]] <> " (exact)", ExactlyZeroQ[(legacyForm[swaEquations, j[[2]]] /. {aDhk -> 1, ank -> 1, aDphik -> 1}) - ref], True]],
  {j, {{"gEQ", gk}, {"DphiEQ", Dphik aDphik}, {"NEQ", nk ank}, {"DhEQ", Dhk aDhk}, {"ANEQ", nk Ak}, {"KEQ", Kk}}}];

Do[test["DE4 vs flowsNEW " <> j[[1]] <> " (exact)", ExactlyZeroQ[legacyForm[swaEquations, j[[2]]] - legacyRef["flowsNEW", j[[1]]]], True],
  {j, {{"gEQ", gk}, {"DphiaDphiEQ", Dphik aDphik}, {"nanEQ", nk ank}, {"DhaDhEQ", Dhk aDhk}, {"AnEQ", nk Ak}, {"KEQ", Kk}, {"HEQ", Hk lk}, {"g2EQ", g2k}}}];
