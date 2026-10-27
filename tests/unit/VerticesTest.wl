(* ::Package:: *)

(* Vertices: Fourier transform, symmetrisation, leg convention, checks *)

vTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
vAnsatz[d_, c_] := quiet[DefineAnsatz[vTheory, "Density" -> d, "Couplings" -> c]];
kpz = quiet[BuildVertices[vAnsatz[hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, {nuk, lamk, Dk}]]];

test["KPZ Gamma2", Gamma2Matrix[kpz, ω, q], {{0, I ω + nuk scal[q, q]}, {-I ω + nuk scal[q, q], -2 Dk}}];
test["KPZ Gamma3 h h hbar (symmetrised: factor 2)", VertexValue[kpz, {h, h, hbar}, {{ω1, q1}, {ω2, q2}}], lamk scal[q1, q2]];
test["KPZ Gamma3 leg order: hbar first, last leg = -(sum)", VertexValue[kpz, {hbar, h, h}, {{ω1, q1}, {ω2, q2}}],
  Expand[-lamk scal[q2, q1 + q2]]];
test["KPZ Gamma3 matrix entry", Gamma3Matrix[kpz, 1, {ν, p}, {ω, q}][[1, 2]], lamk scal[p, q]];
test["KPZ vertex symmetric in identical legs", Expand[VertexValue[kpz, {h, h, hbar}, {{ω1, q1}, {ω2, q2}}] - VertexValue[kpz, {h, h, hbar}, {{ω2, q2}, {ω1, q1}}]], 0];
test["KPZ absent vertex is 0", VertexValue[kpz, {h, h, h}, {{ω1, q1}, {ω2, q2}}], 0];
test["KPZ: no warnings", kpz["Warnings"], {}];

quartic = quiet[BuildVertices[vAnsatz[hbar dt[h] + u4 hbar h^3 + nuk hbar dt[dt[h]], {u4, nuk}]]];
test["Gamma4 h h h hbar = 3! u4", VertexValue[quartic, {h, h, h, hbar}, {{0, q1}, {0, q2}, {0, q3}}], 6 u4];
test["Gamma4 warning on the card", MemberQ[quartic["Warnings"][[All, "Check"]], "Vertices.gamma4"], True];
test["dt dt gives (I w)^2", Gamma2Matrix[quartic, ω, q][[1, 2]], I ω - nuk ω^2];

(* E1: the time-mirrored or anti-damped theory has its free pole in the lower half plane *)
expectFail["anti-damping sign aborts (E1)", BuildVertices[vAnsatz[hbar (dt[h] + nuk lap[h]), {nuk}]], "Vertices.retardedPole[h]"];
expectFail["time-mirrored convention aborts (E1)", BuildVertices[vAnsatz[hbar (-dt[h] - nuk lap[h]), {nuk}]], "Vertices.retardedPole[h]"];
test["no damping: convention check is a warning", MemberQ[quiet[BuildVertices[vAnsatz[hbar dt[h] - Dk hbar^2, {Dk}]]]["Warnings"][[All, "Check"]],
   "Vertices.retardedPole[h]"], True];

(* a deliberately wrong vertex: monomials that do not match the density (simulated parser bug) *)
Module[{a = vAnsatz[hbar (dt[h] - nuk lap[h]), {nuk}], wrong},
  wrong = FRGscalarMSR`Internal`stampStage[KeyDrop[MapAt[ReplacePart[#, {2, "Coefficient"} -> 2 nuk] &, a, "Monomials"], "Fingerprint"]];
  expectFail["deliberately wrong vertex aborts (plane-wave check)", BuildVertices[wrong], "Vertices.planeWave"]];
expectFail["wrong input stage", BuildVertices[vTheory], "Stage.input"];
expectFail["VertexValue needs n-1 legs", VertexValue[kpz, {h, h, hbar}, {{ω, q}}], "Vertices.legCount"];

(* jac (design/LEVICIVITA.md F3): eps_ij d_i A d_j B at legs (l, m) -> -wedge[q_l, q_m] *)
jvTheory = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> FixedDimension[2]]];
jv = quiet[BuildVertices[quiet[DefineAnsatz[jvTheory, "Density" -> Hbar dt[H] + Psibar dt[Psi] + c3 Hbar jac[H, Psi], "Couplings" -> {c3}]]]];
test["jac: Gamma3_{Hbar H Psi} = -c wedge[q_H, q_Psi] (by hand)", VertexValue[jv, {Hbar, H, Psi}, {{ω1, q1}, {ω2, q2}}], -c3 wedge[q2, -q1 - q2]];
test["jac: Gamma3 leg order", VertexValue[jv, {H, Psi, Hbar}, {{ω1, q1}, {ω2, q2}}], -c3 wedge[q1, q2]];
test["jac: Gamma2 unaffected (wedge[q, -q] = 0)", Gamma2Matrix[jv, ω, q], {{0, I ω, 0, 0}, {-I ω, 0, 0, 0}, {0, 0, 0, I ω}, {0, 0, -I ω, 0}}];
test["jac: 2D plane-wave check on the card", StringContainsQ[SelectFirst[jv["Checks"], #Check === "Vertices.planeWave" &]["Detail"], "d = 2 components"], True];

(* LO: fn[f, X] gives f of the momentum of the leg of X; the plane-wave route sees the same (design/LOPIPELINE.md step 3) *)
vLOTheory = quiet[DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]];
vLO = quiet[BuildVertices[DefineAnsatz[vLOTheory, "Density" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],
     "Couplings" -> {nuk, lamk, Dk}, "Functions" -> {fD -> Dk, fnu -> nuk}]]];
test["LO: Gamma2 = {{0, I w + nuk q^2 fnu(q^2)}, {-I w + nuk q^2 fnu(q^2), -2 Dk fD(q^2)}}", Expand[Gamma2Matrix[vLO, ω, q] /. scal[q, q] -> qq2],
  Expand[{{0, I ω + nuk qq2 fnu[qq2]}, {-I ω + nuk qq2 fnu[qq2], -2 Dk fD[qq2]}}]];
test["LO: the Gamma3 vertex is unchanged (k-only constant)", VertexValue[vLO, {hbar, h, h}, {{0, p1}, {0, p2}}] === VertexValue[kpz, {hbar, h, h}, {{0, p1}, {0, p2}}], True];
test["LO: card", Intersection[vLO["Checks"][[All, "Check"]], {"Vertices.planeWave", "Vertices.gamma2Functions", "Vertices.gamma2Polynomial"}],
  {"Vertices.gamma2Functions", "Vertices.planeWave"}];
expectFail["LO: a Gamma2 coefficient that is not reference x F (internal check, forged)",
  FRGscalarMSR`Private`loGamma2Check[{{0, nuk fnu[scal[q, q]] + lamk fnu[scal[q, q]]}}, ω, q, {fnu -> nuk}], "Vertices.gamma2Functions"];
