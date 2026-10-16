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
