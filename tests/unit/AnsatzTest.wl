(* ::Package:: *)

(* Ansatz: parser, Leibniz expansion, echo table, every assertion *)

swTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]];
swCouplings = {nk, Ak, Hk, lk, ank, gk, g2k, Dhk, aDhk, Dphik, aDphik, D2phik, Kk};
swDensity = Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]]) +
   Psibar (dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi]) -
   Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2 - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar];
sw = quiet[DefineAnsatz[swTheory, "Density" -> swDensity, "Couplings" -> swCouplings]];
table[a_] := Sort[{#Coefficient, MonomialString[#]} & /@ a["Monomials"]];

(* expected table, written by hand from the density (design section 4) *)
test["shallow water: 14 monomials", Length[sw["Monomials"]], 14];
test["shallow water: monomial table", table[sw], Sort[{
    {1, "Hbar dt[H]"}, {-Ak nk, "Hbar lap[H]"}, {Hk lk, "Hbar lap[Psi]"}, {lk, "Hbar grad[H].grad[Psi]"}, {lk, "Hbar H lap[Psi]"},
    {1, "Psibar dt[Psi]"}, {-ank nk, "Psibar lap[Psi]"}, {gk, "Psibar H"}, {-g2k, "Psibar lap[H]"}, {lk/2, "Psibar grad[Psi].grad[Psi]"},
    {-aDhk Dhk, "grad[Hbar].grad[Hbar]"}, {-aDphik Dphik, "Psibar Psibar"}, {-D2phik Dphik, "grad[Psibar].grad[Psibar]"},
    {Kk, "grad[Hbar].grad[Psibar]"}}]];
test["shallow water: orders of Hbar grad[H].grad[Psi]",
  Lookup[SelectFirst[sw["Monomials"], MonomialString[#] === "Hbar grad[H].grad[Psi]" &], {"FieldCount", "ResponseCount", "TimeOrder", "MomentumOrder"}],
  {3, 1, 0, 2}];
test["shallow water: stage chained to theory", sw["Upstream"], {swTheory["Fingerprint"]}];

(* Leibniz expansion *)
parse[th_, d_, c_] := table[quiet[DefineAnsatz[th, "Density" -> d, "Couplings" -> c]]];
test["div[H grad Psi] = grad H.grad Psi + H lap Psi",
  parse[swTheory, Hbar dt[H] + lk Hbar div[H grad[Psi]], {lk}], parse[swTheory, Hbar dt[H] + lk Hbar (grad[H].grad[Psi] + H lap[Psi]), {lk}]];
test["lap of a product", parse[swTheory, Hbar lap[H Psi], {}], Sort[{{1, "Hbar lap[H] Psi"}, {2, "Hbar grad[H].grad[Psi]"}, {1, "Hbar H lap[Psi]"}}]];
test["dt of a product", parse[swTheory, Hbar dt[H Psi], {}], Sort[{{1, "Hbar dt[H] Psi"}, {1, "Hbar H dt[Psi]"}}]];
test["dt dt", parse[swTheory, Hbar dt[dt[H]], {}], {{1, "Hbar dt[dt[H]]"}}];
test["grad of a power", parse[swTheory, Hbar div[grad[H^2]], {}], Sort[{{2, "Hbar grad[H].grad[H]"}, {2, "Hbar H lap[H]"}}]];
test["like monomials combine", parse[swTheory, Hbar dt[H] + 2 grad[Hbar].grad[H] + grad[H].grad[Hbar], {}], Sort[{{1, "Hbar dt[H]"}, {3, "grad[Hbar].grad[H]"}}]];
test["cancelling monomials are dropped", parse[swTheory, Hbar dt[H] + Hbar lap[H] - Hbar div[grad[H]], {}], {{1, "Hbar dt[H]"}}];
test["quartic allowed", Length[parse[swTheory, Hbar dt[H] + u4 Hbar H^3, {u4}]], 2];

(* KPZ, 1 scalar *)
kpzTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
test["KPZ table", parse[kpzTheory, hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, {nuk, lamk, Dk}],
  Sort[{{1, "hbar dt[h]"}, {-nuk, "hbar lap[h]"}, {-lamk/2, "hbar grad[h].grad[h]"}, {-Dk, "hbar hbar"}}]];

(* negative tests: each must abort through the named check *)
bad[d_, c_ : {nk}] := DefineAnsatz[swTheory, "Density" -> d, "Couplings" -> c];
expectFail["System Dt (E7)", bad[Hbar Dt[H]], "Ansatz.noSystemDerivatives"];
expectFail["System Dt of coupling x field (E7)", bad[Hbar Dt[nk H]], "Ansatz.noSystemDerivatives"];
expectFail["dt of a coupling", bad[Hbar dt[nk H]], "Ansatz.couplingsOutsideOperators"];
expectFail["undeclared symbol", bad[Hbar dt[H] + nk zz Hbar lap[H]], "Ansatz.knownSymbols"];
expectFail["explicit t dependence", bad[Hbar dt[H] + nk t Hbar H], "Ansatz.knownSymbols"];
expectFail["no response field (MSR)", bad[Hbar dt[H] + nk H lap[H]], "Ansatz.MSR"];
expectFail["linear term", bad[Hbar dt[H] + nk Hbar], "Ansatz.noLinear"];
expectFail["quintic term", bad[Hbar dt[H] + nk Hbar H^4], "Ansatz.quartic"];
expectFail["uncontracted gradient", bad[Hbar dt[H] + nk Hbar grad[H]], "Ansatz.vectorSum"];
expectFail["vector density", bad[nk Hbar grad[H]], "Ansatz.scalarDensity"];
expectFail["product of two gradients", bad[Hbar dt[H] + nk grad[Hbar] grad[H]], "Ansatz.vectorProduct"];
expectFail["lap of a vector", bad[Hbar dt[H] + nk lap[grad[H]].grad[Hbar]], "Ansatz.vectorRank"];
expectFail["div of a scalar", bad[Hbar dt[H] + nk Hbar div[H]], "Ansatz.vectorRank"];
expectFail["non-polynomial", bad[Hbar dt[H] + nk Hbar Exp[H]], "Ansatz.polynomial"];
expectFail["unknown function of a field", bad[Hbar dt[H] + nk Hbar Sin[H]], "Ansatz.unknownStructure"];
expectFail["negative power of a field", bad[Hbar dt[H] + nk Hbar/H], "Ansatz.polynomial"];
expectFail["operator without field", bad[Hbar dt[H] + Hbar H lap[2]], "Ansatz.operatorOnField"];
expectFail["two-argument operator", bad[Hbar dt[H, H]], "Ansatz.operatorArity"];
expectFail["coupling named like a field", bad[Hbar dt[H], {H}], "Ansatz.couplingsNotFields"];
expectFail["reserved coupling name", bad[Hbar dt[H], {q}], "Ansatz.couplings.reserved"];
expectFail["pipeline 2 not in v1", DefineAnsatz[swTheory, "Density" -> Hbar dt[H], "Couplings" -> {}, "Pipeline" -> "MomentumDependent"], "Ansatz.pipeline"];
expectFail["wrong input stage", DefineAnsatz[sw, "Density" -> Hbar dt[H], "Couplings" -> {}], "Stage.input"];
H = 1;
expectFail["stale theory (a field got a value)", DefineAnsatz[swTheory, "Density" -> Hbar dt[H], "Couplings" -> {}], "Stage.fingerprint"];
ClearAll[H];

(* jac: eps_ij d_i a d_j b in d = 2 (design/LEVICIVITA.md) *)
strip[a_Association] := strip /@ KeyDrop[a, {"Timing", "Fingerprint", "Upstream", "Checks", "Warnings"}]; strip[l_List] := strip /@ l; strip[x_] := x;
Module[{kt = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]], ka},
  ka = quiet[DefineAnsatz[kt, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
  (* content hash (without timings, which enter the fingerprint through the nested Theory stage, and without the card lines),
     the same at commit ba4a80d (before the p^(2n) projection) and after *)
  test["eps-free symbolic-Dim ansatz unchanged by jac and p^(2n) support", Hash[strip[ka], "SHA256", "HexString"],
   "a19bb9df3ce9e7f9f07050ecf3953d40b6a2f6fc816d6bb05d4a206c6a6b1f1c"];
  test["eps-free ansatz has no Epsilon key", KeyExistsQ[ka, "Epsilon"], False]];

jTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]];
jAnsatz[d_, c_ : {lk}] := DefineAnsatz[jTheory, "Density" -> Hbar dt[H] + Psibar dt[Psi] + d, "Couplings" -> c];
jTable[a_] := Sort[{#Coefficient, MonomialString[#]} & /@ a["Monomials"]];

expectFail["jac with a symbolic Dim", DefineAnsatz[swTheory, "Density" -> Hbar dt[H] + lk Hbar jac[H, Psi], "Couplings" -> {lk}],
  "Ansatz.epsilonDimension"];
Module[{a = quiet[jAnsatz[lk Hbar jac[H, Psi]]], m},
  m = SelectFirst[a["Monomials"], #FieldCount == 3 &];
  test["jac parsed: one cubic monomial", Count[a["Monomials"], KeyValuePattern["FieldCount" -> 3]], 1];
  test["jac monomial: coefficient and display", {m["Coefficient"], MonomialString[m]}, {lk, "Hbar jac[H, Psi]"}];
  test["jac monomial: orders", Lookup[m, {"TimeOrder", "MomentumOrder"}], {0, 2}];
  test["jac: Epsilon key", a["Epsilon"], True]];
test["jac[B, A] = -jac[A, B]", Select[jTable[quiet[jAnsatz[lk Hbar jac[Psi, H]]]], StringContainsQ[#[[2]], "jac"] &], {{-lk, "Hbar jac[H, Psi]"}}];
Module[{a = quiet[jAnsatz[lk Hbar (jac[H, Psi] + jac[Psi, H])]]},
  test["Hbar (jac[A,B] + jac[B,A]) cancels", Length[a["Monomials"]], 2];
  test["cancellation reported", MemberQ[a["Warnings"][[All, "Check"]], "Ansatz.cancelled"], True];
  test["no Epsilon key after cancellation", KeyExistsQ[a, "Epsilon"], False]];
test["Hbar jac[Psi, Psi] cancels", Length[quiet[jAnsatz[lk Hbar jac[Psi, Psi]]]["Monomials"]], 2];
test["jac[dt[H], Psi] parsed", Select[jTable[quiet[jAnsatz[lk Hbar jac[dt[H], Psi]]]], StringContainsQ[#[[2]], "jac"] &],
  {{lk, "Hbar jac[dt[H], Psi]"}}];
test["jac of a product: Leibniz", Select[jTable[quiet[jAnsatz[lk Hbar jac[H Psi, Psi]]]], StringContainsQ[#[[2]], "jac"] &],
  {{lk, "Hbar Psi jac[H, Psi]"}}];
expectFail["jac[grad[H], Psi]", jAnsatz[lk Hbar jac[grad[H], Psi]], "Ansatz.vectorRank"];
expectFail["jac[H]", jAnsatz[lk Hbar jac[H]], "Ansatz.operatorArity"];
expectFail["jac[H, nk Psi]", jAnsatz[Hbar jac[H, nk Psi], {nk}], "Ansatz.couplingsOutsideOperators"];
expectFail["jac[H, 1] (argument without a field)", jAnsatz[lk Hbar jac[H, 1]], "Ansatz.operatorOnField"];
expectFail["two jac in one monomial", jAnsatz[lk Hbar jac[H, Psi] jac[Psi, H]], "Ansatz.singleEpsilon"];
expectFail["jac squared (one label, two slots)", jAnsatz[lk Hbar jac[H + Psi, H - Psi]^2], "Ansatz.singleEpsilon"];
expectFail["quadratic jac (total derivative)", jAnsatz[lk jac[Hbar, H]], "Ansatz.jacQuadratic"];

(* p^4 and beyond (design/MOMENTUMPROJECTION.md step 2): powers get fresh dummies (review E14), canonical form with several
   contractions (E9) *)
p4[d_] := Select[parse[swTheory, Hbar dt[H] + d, {nk}], #[[2]] =!= "Hbar dt[H]" &];
test["p^4 monomial accepted", p4[nk Hbar lap[lap[H]]], {{nk, "Hbar lap[lap[H]]"}}];
test["p^4 monomial: MomentumOrder 4",
  Lookup[SelectFirst[quiet[DefineAnsatz[swTheory, "Density" -> Hbar dt[H] + nk Hbar lap[lap[H]], "Couplings" -> {nk}]]["Monomials"], #FieldCount == 2 && #TimeOrder == 0 &],
   "MomentumOrder"], 4];
test["lap[H]^2 parses", Length[p4[nk Hbar lap[H]^2]], 1];
test["lap[H] lap[H] = lap[H]^2", p4[nk Hbar lap[H] lap[H]], p4[nk Hbar lap[H]^2]];
test["(grad Hbar.grad H)^2 parses to one monomial", Length[p4[nk (grad[Hbar].grad[H])^2]], 1];
test["(grad Hbar.grad H)^2 = product written with swapped dummies",
  p4[nk (grad[Hbar].grad[H])^2 - nk (grad[H].grad[Hbar]) (grad[Hbar].grad[H])], {}];
test["lap[lap[H]] = div[grad[lap[H]]]", p4[nk Hbar lap[lap[H]] - nk Hbar div[grad[lap[H]]]], {}];
test["grad[lap[H]].grad[H] = grad[H].grad[lap[H]]", p4[nk Hbar grad[lap[H]].grad[H] - nk Hbar grad[H].grad[lap[H]]], {}];
test["lap[lap[H]] lap[H]: three contractions", Length[p4[nk H Hbar lap[lap[H]] lap[H]]], 1];
test["jac[H, lap[Psi]] = -jac[lap[Psi], H] (two dummies with eps)",
  Select[jTable[quiet[jAnsatz[lk Hbar jac[H, lap[Psi]] + lk Hbar jac[lap[Psi], H]]]], StringContainsQ[#[[2]], "jac"] &], {}];
test["jac[lap[H], lap[H]] = 0", Select[jTable[quiet[jAnsatz[lk Hbar jac[lap[H], lap[H]]]]], StringContainsQ[#[[2]], "jac"] &], {}];
test["jac[H, lap[Psi]] kept", Length[Select[jTable[quiet[jAnsatz[lk Hbar jac[H, lap[Psi]]]]], StringContainsQ[#[[2]], "jac"] &]], 1];

(* LO: fn[f, X] and "Functions" (design/LOPIPELINE.md step 2) *)
loTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
loKPZ = hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar];
loAn = quiet[DefineAnsatz[loTheory, "Density" -> loKPZ, "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]];
test["LO: pipeline inferred from fn", loAn["Pipeline"], "LO"];
test["LO: functions and coupling kinds", {loAn["Functions"], loAn["CouplingKinds"]}, {{fD -> Dk, fnu -> nuk}, <|nuk -> "Function", lamk -> "Constant", Dk -> "Function"|>}];
test["LO: monomial table", table[loAn], Sort[{{1, "hbar dt[h]"}, {-nuk, "hbar fn[fnu, lap[h]]"}, {-lamk/2, "hbar grad[h].grad[h]"}, {-Dk, "hbar fn[fD, hbar]"}}]];
test["LO: fn[f, lap[h]] = lap[fn[f, h]] (f(-Laplacian) commutes with the Laplacian)",
  quiet[DefineAnsatz[loTheory, "Density" -> hbar (dt[h] - nuk lap[fn[fnu, h]] - lamk/2 grad[h].grad[h]) - Dk fn[fD, hbar] hbar,
     "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]]["Monomials"] === loAn["Monomials"], True];
test["LO: MonomialFunction", Sort[DeleteCases[MonomialFunction /@ loAn["Monomials"], None]], {fD, fnu}];
test["DE: an ansatz without fn has no LO keys", {sw["Pipeline"], KeyExistsQ[sw, "Functions"], KeyExistsQ[sw, "CouplingKinds"]}, {"DerivativeExpansion", False, False}];
expectFail["LO: explicit Pipeline must match the inferred one", DefineAnsatz[loTheory, "Density" -> loKPZ, "Couplings" -> {nuk, lamk, Dk},
   "Functions" -> {fD -> Dk, fnu -> nuk}, "Pipeline" -> "DerivativeExpansion"], "Ansatz.pipeline"];
expectFail["LO: Functions without fn", DefineAnsatz[loTheory, "Density" -> hbar dt[h], "Couplings" -> {}, "Functions" -> {fnu -> nuk}], "Ansatz.functionReference"];
expectFail["LO: undeclared function", DefineAnsatz[loTheory, "Density" -> loKPZ, "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fnu -> nuk}], "Ansatz.functionReference"];
expectFail["LO: coefficient is not number x reference", DefineAnsatz[loTheory, "Density" -> loKPZ, "Couplings" -> {nuk, lamk, Dk},
   "Functions" -> {fD -> nuk, fnu -> Dk}], "Ansatz.functionReference"];
expectFail["LO: fn on a product", DefineAnsatz[loTheory, "Density" -> loKPZ - nuk fn[fnu, h h] hbar, "Couplings" -> {nuk, lamk, Dk},
   "Functions" -> {fD -> Dk, fnu -> nuk}], "Ansatz.fnOneField"];
expectFail["LO: fn in a vertex", DefineAnsatz[loTheory, "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 fn[fl, h] h) - Dk hbar fn[fD, hbar],
   "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk, fl -> lamk}], "Ansatz.fnBilinear"];
expectFail["LO: function times i w", DefineAnsatz[loTheory, "Density" -> hbar (dt[fn[fZ, h]] - nuk fn[fnu, lap[h]]) - Dk hbar fn[fD, hbar],
   "Couplings" -> {nuk, Dk, Zk}, "Functions" -> {fD -> Dk, fnu -> nuk, fZ -> Zk}], "Ansatz.functionTimesOmega"];
expectFail["LO: a DE monomial next to the functions (mixed ansatz)", DefineAnsatz[loTheory, "Density" -> loKPZ + tauk hbar lap[lap[h]],
   "Couplings" -> {nuk, lamk, Dk, tauk}, "Functions" -> {fD -> Dk, fnu -> nuk}], "Ansatz.mixedNotYet"];
