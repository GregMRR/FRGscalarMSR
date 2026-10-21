(* ::Package:: *)

(* Core: FRGAssert, ExactlyZeroQ, scal, stage validation *)

test["FRGAssert pass returns the check name", FRGAssert["x.y", 1 + 1 == 2, Null, ""], "x.y"];
expectFail["FRGAssert failure aborts and records the name", FRGAssert["x.fails", 1 == 2, {1, 2}, "hint"], "x.fails"];

test["ExactlyZeroQ identity", ExactlyZeroQ[(a + b)^2 - a^2 - 2 a b - b^2], True];
test["ExactlyZeroQ rational identity", ExactlyZeroQ[1/(a - b) - 1/(a + b) - 2 b/(a^2 - b^2)], True];
test["ExactlyZeroQ non-zero", ExactlyZeroQ[a - b], False];
test["ExactlyZeroQ randomises Dim (E12)", ExactlyZeroQ[Dim - 3], False];
test["ExactlyZeroQ same shape atom, same value", ExactlyZeroQ[r[q^2/k^2] - r[q^2/k^2]], True];
test["ExactlyZeroQ r and r' independent (E12)", ExactlyZeroQ[r'[y0] - r[y0]], False];
test["ExactlyZeroQ r' and r'' independent (E12)", ExactlyZeroQ[r''[y0] - r'[y0]], False];
test["ExactlyZeroQ two shapes independent", ExactlyZeroQ[r1[y0] - r2[y0]], False];
test["ExactlyZeroQ Pi and I are constants", ExactlyZeroQ[I Pi a - Pi I a], True];
(* derivative atoms are seen as wholes (design/MOMENTUMPROJECTION.md step 0, E1) *)
test["ExactlyZeroQ sees GammaAtom'", ExactlyZeroQ[Derivative[1][GammaAtom[hb, h, 0]][QQ]], False];
test["ExactlyZeroQ sees DeltaCoefficient''", ExactlyZeroQ[Derivative[2][DeltaCoefficient[1]][QQ]], False];
test["ExactlyZeroQ GammaAtom' and GammaAtom independent", ExactlyZeroQ[Derivative[1][GammaAtom[hb, h, 0]][QQ] - GammaAtom[hb, h, 0][QQ]], False];
test["ExactlyZeroQ chain rule through an atom", ExactlyZeroQ[D[GammaAtom[a, b, 0][x^2], x] - 2 x Derivative[1][GammaAtom[a, b, 0]][x^2]], True];
expectFail["ExactlyZeroQ refuses a non-numeric leftover", ExactlyZeroQ["text"], "ExactlyZeroQ.numeric"];

test["scal Orderless", scal[b, a] === scal[a, b], True];
test["scal linear over sums", scal[a + b, c], scal[a, c] + scal[b, c]];
test["scal linear over numbers", scal[2 a, -3 b], -6 scal[a, b]];
test["scal not linear over couplings (E14)", FreeQ[scal[nk a, b], scal[a, b]], True];
test["scal of zero", scal[0, a], 0];

Module[{th = quiet[DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]]},
  expectFail["validateStage wrong stage", FRGscalarMSR`Internal`validateStage[th, "Ansatz"], "Stage.input"];
  expectFail["validateStage tampered stage", FRGscalarMSR`Internal`validateStage[Append[th, "Dimension" -> d2], "Theory"], "Stage.fingerprint"];
  test["validateStage untouched stage", FRGscalarMSR`Internal`validateStage[th, "Theory"] === th, True]];

(* wedge: eps_ij a_i b_j in d = 2 (design/LEVICIVITA.md F1, F2) *)
test["wedge antisymmetric", wedge[b, a] + wedge[a, b], 0];
test["wedge canonical order (not Orderless)", wedge[b, a], -wedge[a, b]];
test["wedge[a, a] = 0", wedge[a, a], 0];
test["wedge linear in slot 1", wedge[a + 2 b, c], wedge[a, c] + 2 wedge[b, c]];
test["wedge linear in slot 2", wedge[a, -b + 3 c], -wedge[a, b] + 3 wedge[a, c]];
test["wedge of zero", {wedge[0, a], wedge[a, 0]}, {0, 0}];
test["wedge of a difference with itself", wedge[a + b, a + b], 0];
test["wedge not linear over couplings", FreeQ[wedge[nk a, b], wedge[a, b]], True];
Module[{vecs, num, e},
  SeedRandom[5]; vecs = Thread[{a, b, c, d} -> RandomInteger[{-9, 9}, {4, 2}]];
  num[x_] := x /. {scal[u_, v_] :> (u /. vecs) . (v /. vecs), wedge[u_, v_] :> Module[{s = u /. vecs, t = v /. vecs}, s[[1]] t[[2]] - s[[2]] t[[1]]]};
  e = 3 wedge[a, b] wedge[c, d] + nk wedge[a, c]^2 wedge[b, d] - wedge[a + b, c]^2 scal[a, d] + wedge[a, d]^3 wedge[b, c];
  test["reduceWedge2D: at most one wedge per term", Max[Exponent[# /. wedge[__] -> ww, ww] & /@ List @@ FRGscalarMSR`Internal`reduceWedge2D[e]], 1];
  test["reduceWedge2D = original on explicit 2-vectors", Expand[num[FRGscalarMSR`Internal`reduceWedge2D[e]] - num[e]], 0];
  test["F2 on explicit 2-vectors", Expand[num[wedge[a, b] wedge[c, d]] - num[scal[a, c] scal[b, d] - scal[a, d] scal[b, c]]], 0]];
test["displayString wedge", FRGscalarMSR`Internal`displayString[wedge[q1, q2] + scal[q1, q1]], "q1^2 + wedge(q1,q2)"];
