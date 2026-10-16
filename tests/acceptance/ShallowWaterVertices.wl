(* ::Package:: *)

(* Milestone 2 acceptance: the shallow-water DE4plus vertices derived from the density equal the hand-written
   Gamma2 / Gamma3 / AnsatzBMW / form* of DE4plus_fast.wl exactly (regulator functions set to 0). *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

test["shallow water: Gamma2 = script Gamma2 /. AnsatzBMW /. form* (regulator off)",
  Expand[Gamma2Matrix[swaVertices, ω, q] - (refGamma2[ω, q] /. AnsatzBMW /. refForms /. regOff)], ConstantArray[0, {4, 4}]];
Do[test["shallow water: Gamma3[" <> ToString[i] <> "] = script Gamma3[" <> ToString[i] <> "] /. AnsatzBMW",
   Expand[Gamma3Matrix[swaVertices, i, {ν, p}, {ω, q}] - (refGamma3[i][ν, p, ω, q] /. AnsatzBMW)], ConstantArray[0, {4, 4}]], {i, 4}];
test["shallow water: no Gamma4", Select[Keys[swaVertices["VertexTable"]], Length[#] == 4 &], {}];
test["shallow water: all vertex checks passed", Count[swaVertices["Checks"], KeyValuePattern["Status" -> "PASS"]], 6];
