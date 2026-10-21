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
expectFail["p^4 monomial (E2)", bad[Hbar dt[H] + nk Hbar lap[lap[H]]], "Ansatz.momentumOrder"];
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
strip[a_Association] := strip /@ KeyDrop[a, {"Timing", "Fingerprint", "Upstream"}]; strip[l_List] := strip /@ l; strip[x_] := x;
Module[{kt = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]], ka},
  ka = quiet[DefineAnsatz[kt, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}]];
  (* content hash (without timings, which enter the fingerprint through the nested Theory stage) recorded at commit 1fbb12e *)
  test["eps-free symbolic-Dim ansatz unchanged by jac support", Hash[strip[ka], "SHA256", "HexString"],
   "2821c78a100ebada39ecf39e511201ef4a3991401eacd155c3852c890d69cee4"];
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
