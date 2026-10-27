(* ::Package:: *)

(* Dimensions: derived dimensions against hand values (old notebook AdimensionalFlowEQGenericRegLargeScale.nb, section 2),
   KPZ, the LPA ansatz, the 2D jac example, and every abort *)

dimTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];

(* DE4plus: hand values of the old notebook (Dimensionality[x, y, z, u] = k^x Dphik^y Dhk^z nk^u) *)
dimDE4 = quiet[DefineAnsatz[dimTheory, "Density" ->
      Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) +
      Psibar (dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi]) -
      Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2 - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar],
    "Couplings" -> {nk, Ak, Hk, lk, ank, gk, g2k, Dhk, aDhk, Dphik, aDphik, D2phik, Kk}]];
dimsDE4 = quiet[DeriveDimensions[dimDE4, "References" -> {Dphik -> "phi", Dhk -> "h", nk -> "nu"}, "Dimensionless" -> {Ak, ank, aDhk, aDphik}]];
test["DE4plus: old notebook hand values (g, g2, D2phi, K, H, l)",
  quiet[CheckDimensions[dimsDE4, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, g2k -> k^-1 Dphik^(1/2) Dhk^(-1/2) nk, D2phik -> k^-2,
     Kk -> k^-1 Dphik^(1/2) Dhk^(1/2), Hk -> k^(Dim/2) Dhk^(1/2) nk^(-1/2), lk -> k^(1 - Dim/2) Dphik^(-1/2) nk^(3/2)|>]],
  {gk, g2k, D2phik, Kk, Hk, lk}];
test["DE4plus: frequency k^2 nk", quiet[DimensionOf[dimsDE4, "Frequency"]], k^2 nk];
test["DE4plus: linear parts = old betal (d = 2)", Expand[quiet[LinearPart[dimsDE4, #]] /. Dim -> 2] & /@ {gk, g2k, lk, Hk, D2phik, Kk},
  Expand[{-1 - etah/2 + etanu + etaphi/2, 1 - etah/2 + etanu + etaphi/2, (3 etanu)/2 - etaphi/2, -1 + etah/2 - etanu/2, 2, 1 + etah/2 + etaphi/2}]];
test["DE4plus: names", Lookup[dimsDE4["Names"], {gk, Ak, D2phik}], {agk, Ak, aD2phik}];
test["DE4plus: stage chained to the ansatz", dimsDE4["Upstream"], {dimDE4["Fingerprint"]}];

(* LPA (ShallowWaterLPA.nb): mk Psibar dt[H] and Dmk grad[Hbar].grad[Psibar] *)
dimLPA = quiet[DefineAnsatz[dimTheory, "Density" ->
      Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) +
      Psibar (dt[Psi] + mk dt[H] - nk lap[Psi] + gk H + lk/2 grad[Psi].grad[Psi]) -
      Dhk grad[Hbar].grad[Hbar] - Dmk grad[Hbar].grad[Psibar] - Dphik aDphik Psibar^2,
    "Couplings" -> {nk, Ak, Hk, lk, gk, mk, Dhk, Dmk, Dphik, aDphik}]];
dimsLPA = quiet[DeriveDimensions[dimLPA, "References" -> {nk -> "nu", Dhk -> "h", Dphik -> "phi"}, "Dimensionless" -> {Ak, aDphik}]];
test["LPA: hand values of mk, Dmk, gk", quiet[CheckDimensions[dimsLPA,
    <|mk -> Dphik^(1/2) Dhk^(-1/2) k^-1, Dmk -> Dhk^(1/2) Dphik^(1/2) k^-1, gk -> k Dphik^(1/2) Dhk^(-1/2) nk|>]], {mk, Dmk, gk}];

(* KPZ: lambda^2 D / nu^3 has dimension k^(2-d) *)
dimKPZTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
dimKPZ = quiet[DefineAnsatz[dimKPZTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
dimsKPZ = quiet[DeriveDimensions[dimKPZ, "References" -> {nuk -> "nu", Dk -> "D"}]];
test["KPZ: [lambda] = k^(1-d/2) nu^(3/2) D^(-1/2)", quiet[CheckDimensions[dimsKPZ, <|lamk -> k^(1 - Dim/2) nuk^(3/2) Dk^(-1/2), "Frequency" -> k^2 nuk|>]],
  {lamk, "Frequency"}];
test["KPZ: g = lambda^2 D/nu^3 scales as k^(2-d)", PowerExpand[quiet[DimensionOf[dimsKPZ, lamk]]^2 Dk/nuk^3], k^(2 - Dim)];

(* 2D jac term (README example), declared d = 2 *)
dimJacTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]];
dimJac = quiet[DefineAnsatz[dimJacTheory, "Density" ->
      Hbar (dt[H] + lk jac[Psi, H] + kap grad[Psi].grad[H] - nk lap[H]) - dk Hbar^2 + Psibar (dt[Psi] - mk lap[Psi]) - ek Psibar^2,
    "Couplings" -> {lk, kap, nk, dk, mk, ek}]];
dimsJac = quiet[DeriveDimensions[dimJac, "References" -> {nk -> "nu", dk -> "d", ek -> "e"}]];
test["jac: each eps slot counts one gradient ([lk] = [kap])", quiet[DimensionOf[dimsJac, lk]] === quiet[DimensionOf[dimsJac, kap]], True];
test["jac: [mk] = [nk]", quiet[DimensionOf[dimsJac, mk]], nk];
test["jac: exponents numeric for a declared d = 2", FreeQ[dimsJac["Vectors"], Dim], True];

(* aborts *)
expectFail["two references for N = 2", DeriveDimensions[dimDE4, "References" -> {Dphik -> "phi", nk -> "nu"}, "Dimensionless" -> {Ak, ank, aDhk, aDphik}],
  "Dimensions.referenceCount"];
expectFail["amplitudes not declared", DeriveDimensions[dimDE4, "References" -> {Dphik -> "phi", Dhk -> "h", nk -> "nu"}], "Dimensions.freedoms"];
expectFail["dependent references ([mk] = [nk])", DeriveDimensions[dimJac, "References" -> {nk -> "nu", mk -> "m", ek -> "e"}],
  "Dimensions.referencesIndependent"];
expectFail["inconsistent monomials (nuk in hbar lap[h] and hbar h)",
  DeriveDimensions[quiet[DefineAnsatz[dimKPZTheory, "Density" -> hbar (dt[h] - nuk lap[h] + nuk h - lamk/2 grad[h].grad[h]) - Dk hbar^2,
     "Couplings" -> {nuk, lamk, Dk}]], "References" -> {nuk -> "nu", Dk -> "D"}], "Dimensions.consistent"];
expectFail["hand value mismatch", CheckDimensions[dimsKPZ, <|lamk -> k^(1 + Dim/2) nuk^(3/2) Dk^(-1/2)|>], "Dimensions.handCheck"];
expectFail["generated name collides (lamk -> alamk exists)",
  DeriveDimensions[quiet[DefineAnsatz[dimKPZTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk alamk/2 grad[h].grad[h]) - Dk hbar^2,
     "Couplings" -> {nuk, lamk, alamk, Dk}]], "References" -> {nuk -> "nu", Dk -> "D"}, "Dimensionless" -> {alamk}], "Dimensions.names"];
test["Names override resolves the collision",
  quiet[DeriveDimensions[quiet[DefineAnsatz[dimKPZTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk alamk/2 grad[h].grad[h]) - Dk hbar^2,
       "Couplings" -> {nuk, lamk, alamk, Dk}]], "References" -> {nuk -> "nu", Dk -> "D"}, "Dimensionless" -> {alamk}, "Names" -> <|lamk -> lamhat|>]]["Names"][lamk],
  lamhat];
expectFail["dimensionless reference", DeriveDimensions[dimKPZ, "References" -> {nuk -> "nu", Dk -> "D"}, "Dimensionless" -> {Dk}],
  "Dimensions.referenceDimensionless"];
test["AnsatzOf finds the embedded ansatz", AnsatzOf[<|"Stage" -> "X", "Inner" -> <|"Stage" -> "Y", "A" -> dimKPZ|>|>] === dimKPZ, True];

(* LO: functions are dimensionless, so the KPZ LO couplings have the KPZ dimensions (design/LOPIPELINE.md step 7) *)
dimLOTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
dimLOAnsatz = quiet[DefineAnsatz[dimLOTheory, "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
    "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]];
dimLO = quiet[DeriveDimensions[dimLOAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]];
dimKPZDE = quiet[DeriveDimensions[quiet[DefineAnsatz[dimLOTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2,
      "Couplings" -> {nuk, lamk, Dk}]], "References" -> {nuk -> "nu", Dk -> "D"}]];
test["LO: dimensions = those of the DE KPZ ansatz", {dimLO["Vectors"], dimLO["Time"]}, {dimKPZDE["Vectors"], dimKPZDE["Time"]}];
test["LO: rescaling check (fn[f, X] scales as X) and functionReferences on the card",
  SubsetQ[dimLO["Checks"][[All, "Check"]], {"Dimensions.rescaling", "Dimensions.functionReferences"}], True];
test["LO: Pipeline and Functions keys", {dimLO["Pipeline"], dimLO["Functions"]}, {"LO", {fD -> Dk, fnu -> nuk}}];
test["DE: no Functions key", KeyExistsQ[dimKPZDE, "Functions"], False];
expectFail["LO: a function's reference must be a running reference", DeriveDimensions[dimLOAnsatz, "References" -> {nuk -> "nu", Dk -> {"D", "Mode" -> "Fixed"}}],
  "Dimensions.functionReferences"];
