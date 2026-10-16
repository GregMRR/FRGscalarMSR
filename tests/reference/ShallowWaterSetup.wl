(* ::Package:: *)

(* shallow-water DE4plus pipeline stages of design section 4, shared by the acceptance tests (built once per run) *)

If[! ValueQ[swaRegulator],
  Get[FileNameJoin[{testRoot, "reference", "DE4plusReference.wl"}]];
  swaTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
  swaAnsatz = quiet[DefineAnsatz[swaTheory, "Density" ->
        Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) +
        Psibar (dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi]) -
        Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2 - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar],
      "Couplings" -> {nk, Ak, Hk, lk, ank, gk, g2k, Dhk, aDhk, Dphik, aDphik, D2phik, Kk}]];
  swaVertices = quiet[BuildVertices[swaAnsatz]];
  swaRegulator = quiet[DefineRegulator[swaAnsatz,
     {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2],
     {Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2],
     "Shapes" -> {arnh, arnphi, armh, armphi}]];
  swaPropagator = quiet[BuildPropagator[swaVertices, swaRegulator]];
  swaSymmetries = quiet[DeclareSymmetries[swaAnsatz,
     Protected[{Hbar, _}, "Momentum" -> 0, "Reason" -> "Hbar -> Hbar + c(t) (mass conservation)"],
     Protected[{Psi, _}, "Momentum" -> 0, "Reason" -> "Psi -> Psi + c(t)"],
     Neglected[{H, Psibar}, "Frequency" -> {1, 2}, "Momentum" -> 0, "Reason" -> "truncation, z = 1 power counting"],
     Neglected[{Psibar, Psibar}, "Frequency" -> 2, "Momentum" -> 0, "Reason" -> "coloured forcing, dropped"],
     Neglected[_, "Frequency" -> {1, 2}, "Momentum" -> 1, "Reason" -> "mixed i w q^2 / w^2 q^2 terms, dropped"],
     NonRenormalised[lk, "Reason" -> "Galilean invariance"],
     "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
  swaLoops = quiet[BuildLoopIntegrands[swaVertices, swaPropagator, swaRegulator]];
  swaProjected = quiet[ProjectOnMonomials[swaLoops, swaSymmetries]];
  swaIntegrated = quiet[IntegrateFrequencies[swaProjected]];
  (* script atoms in terms of ours (E10 and the primer's Gamma2 table) *)
  toScriptAtoms = {GammaAtom[H, Hbar, 0] -> feta, GammaAtom[Psi, Psibar, 0] -> fgamma, GammaAtom[Hbar, Hbar, 0] -> (-fDh[#] &),
     GammaAtom[Psibar, Psibar, 0] -> (-fDphi[#] &), GammaAtom[H, Psibar, 0] -> fl, GammaAtom[Hbar, Psi, 0] -> (-fh[#] &),
     GammaAtom[Hbar, Psibar, 0] -> fk, \[CapitalDelta] -> refDelta};
  refForms = {feta[s_] :> formeta[s], fDh[s_] :> formDh[s], fgamma[s_] :> formgamma[s], fDphi[s_] :> formDphi[s],
     fk[s_] :> formk[s], fl[s_] :> forml[s], fh[s_] :> formh[s]};
  regOff = Thread[{arnh, armh, arnphi, armphi} -> (0 &)]];
