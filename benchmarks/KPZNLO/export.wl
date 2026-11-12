(* ::Package:: *)

(* export.wl -- writes the two NLO BMWflow modules (contract v4) into out/ (git-ignored): out/KPZ_NLO.py (f_lambda = 1: KCW LO =
   BMWflow ModelKPZ "NLO") and out/KPZLambda_NLO.py (flowing f_lambda: KCW NLO), and out/KPZLambda_NLO_points.json: the Mathematica
   values of both modules' loop rows (computed with 30 digits, written as doubles) at 5 random points, 2 of them at w = 0, for BMWflow's tests (tests/fixtures/).
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZNLO/export.wl *)

kpznloDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[kpznloDir, 2], "FRGscalarMSR.wl"}]];
Block[{Print}, Get[FileNameJoin[{ParentDirectory[kpznloDir], #, "Pipeline.wl"}]] & /@ {"KPZLO", "KPZLOLambda"};
  Get[FileNameJoin[{kpznloDir, "Pipeline.wl"}]]];
ExportPython[kpznloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{kpznloDir, "out", "KPZ_NLO.py"}]];
ExportPython[kpzlamnloAdim, "G" -> alamk^2, "Model" -> "KPZLambda", "File" -> FileNameJoin[{kpznloDir, "out", "KPZLambda_NLO.py"}]];

(* reference points: every loop argument (dyadic rationals), the rows of the NLO module and, at w = 0, of the LO module (KPZLOLambda) *)
kpznloPoints = Module[{parts = PythonParts[kpzlamnloAdim, alamk^2], partsLO = PythonParts[kpzlamAdim, alamk^2], dy, tup, pts},
   SeedRandom[20261112];
   dy[] := RandomInteger[{16, 127}]/128; tup[n_] := Table[dy[], n];
   pts = Table[With[{p = dy[], q = dy[], c = RandomInteger[{-120, 120}]/128},
      <|"p" -> p, "w" -> If[j <= 2, 0, 4 dy[]], "q" -> q, "c" -> c, "Q2" -> p^2 + q^2 + 2 p q c, "fp" -> tup[3], "fq" -> tup[3], "fQ" -> tup[3],
       "rq" -> tup[2], "rq_" -> -tup[2], "rQ" -> tup[2], "eta" -> tup[2], "couplings" -> tup[1], "d" -> 1|>], {j, 5}];
   Function[pt, Module[{val},
      val[e_] := N[e /. {PythonAtom[a_, i_, b_, j_] :> pt[a][[i + 1]] + pt[b][[j + 1]], PythonVariable[s_, i_] :> pt[s][[i + 1]], PythonVariable[s_] :> pt[s]}, 30];
      Join[N[pt], <|"loop" -> val[parts["loop"]], "loopLO" -> If[pt["w"] == 0, val[partsLO["loop"]], Null]|>]]] /@ pts];
Export[FileNameJoin[{kpznloDir, "out", "KPZLambda_NLO_points.json"}], N[kpznloPoints], "RawJSON"];
Print["wrote out/KPZ_NLO.py, out/KPZLambda_NLO.py, out/KPZLambda_NLO_points.json"];
