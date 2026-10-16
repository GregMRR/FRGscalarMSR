(* ::Package:: *)

(* Milestone 2 acceptance: the shallow-water DE4plus vertices derived from the density equal the hand-written
   Gamma2 / Gamma3 / AnsatzBMW / form* of DE4plus_fast.wl exactly (regulator functions set to 0). *)

Get[FileNameJoin[{testRoot, "reference", "DE4plusReference.wl"}]];

swaTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
swaAnsatz = quiet[DefineAnsatz[swaTheory, "Density" ->
      Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) +
      Psibar (dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi]) -
      Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2 - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar],
    "Couplings" -> {nk, Ak, Hk, lk, ank, gk, g2k, Dhk, aDhk, Dphik, aDphik, D2phik, Kk}]];
swaVertices = quiet[BuildVertices[swaAnsatz]];

refForms = {feta[s_] :> formeta[s], fDh[s_] :> formDh[s], fgamma[s_] :> formgamma[s], fDphi[s_] :> formDphi[s],
   fk[s_] :> formk[s], fl[s_] :> forml[s], fh[s_] :> formh[s]};
regOff = Thread[{arnh, armh, arnphi, armphi} -> (0 &)];

test["shallow water: Gamma2 = script Gamma2 /. AnsatzBMW /. form* (regulator off)",
  Expand[Gamma2Matrix[swaVertices, ω, q] - (refGamma2[ω, q] /. AnsatzBMW /. refForms /. regOff)], ConstantArray[0, {4, 4}]];
Do[test["shallow water: Gamma3[" <> ToString[i] <> "] = script Gamma3[" <> ToString[i] <> "] /. AnsatzBMW",
   Expand[Gamma3Matrix[swaVertices, i, {ν, p}, {ω, q}] - (refGamma3[i][ν, p, ω, q] /. AnsatzBMW)], ConstantArray[0, {4, 4}]], {i, 4}];
test["shallow water: no Gamma4", Select[Keys[swaVertices["VertexTable"]], Length[#] == 4 &], {}];
test["shallow water: all vertex checks passed", Count[swaVertices["Checks"], KeyValuePattern["Status" -> "PASS"]], 6];
