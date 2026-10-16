(* ::Package:: *)

(* Milestone 3 acceptance: the regulator of design section 4 and the automatically inverted propagator reproduce the
   hand-written forms, Delta, A/CC and G2Explicit (proprules) of DE4plus_fast.wl exactly. *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];

test["shallow water: 7 entry atoms", Length[swaPropagator["Atoms"]], 7];
test["shallow water: atoms (Gamma_k + R) = script form* with regulator",
  Expand[Table[(a /. AtomRules[swaPropagator])[scal[q, q]] - (((a /. toScriptAtoms)[scal[q, q]]) /. refForms), {a, swaPropagator["Atoms"]}]],
  ConstantArray[0, 7]];
test["shallow water: Delta = script Delta (abstract atoms)",
  Expand[(swaPropagator["Delta"][ω, scal[q, q]] /. toScriptAtoms) - ΔofS[ω, scal[q, q]]], 0];
test["shallow water: c2 = 1, c1 = 2A, c0 = CC (script ACsub)",
  Expand[(swaPropagator["DeltaCoefficients"][scal[q, q]] /. toScriptAtoms) - {CC, 2 A, 1} /. ACsub], {0, 0, 0}];
test["shallow water: G = script G2Explicit /. AnsatzBMW (Delta abstract on both sides)",
  Map[Together, (PropagatorMatrix[swaPropagator, ω, scal[q, q]] /. toScriptAtoms) - (refG2Explicit = G2Explicit[ω, q] /. AnsatzBMW), {2}],
  ConstantArray[0, {4, 4}]];
test["shallow water: dsR pattern = script dsRmat under the E10 map",
  DsRAbstractMatrix[swaRegulator] /. {dsR[H, Hbar] -> dsNh, dsR[Hbar, Hbar] -> dsMh, dsR[Psi, Psibar] -> dsNphi, dsR[Psibar, Psibar] -> dsMphi},
  dsRmatComplete];
test["shallow water: validity A > 0, CC > 0", (ValidityConditions[swaPropagator][scal[q, q]] /. DeltaCoefficient[m_][_] :> {CC, 2 A, 1}[[m + 1]]),
  {2 A > 0, CC > 0}];
