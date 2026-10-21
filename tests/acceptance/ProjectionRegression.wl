(* ::Package:: *)

(* Invariant of design/MOMENTUMPROJECTION.md: the p^(2n) projection gives, for every n <= 1 monomial of the existing pipelines,
   integrands Expand-identical to the former p-Laplacian operator (kept test-only in tests/reference/OldProjection.wl).
   The declared-d = 2 and LPA pipelines are compared in their own files. *)

Get[FileNameJoin[{testRoot, "reference", "OldProjection.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];
If[! ValueQ[jacProjected], Get[FileNameJoin[{testRoot, "acceptance", "Jacobian2D.wl"}]]];

test["KPZ: new projection = old operator, every integrand", FRGscalarMSR`Private`oldLapMismatches[kpzIntegrated["Projection"]], {}];
test["jac (eps, d = 2): new projection = old operator, every integrand", FRGscalarMSR`Private`oldLapMismatches[jacProjected], {}];
test["shallow water: new projection = old operator, every integrand", FRGscalarMSR`Private`oldLapMismatches[swaProjected], {}];
test["shallow water: integrands compared", Length[swaProjected["Integrands"]] > 0, True];
