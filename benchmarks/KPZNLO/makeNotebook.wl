(* ::Package:: *)

(* makeNotebook.wl -- writes KPZNLO.nb (this folder) from the cell list below.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZNLO/makeNotebook.wl
   Check that it runs headless:
         math -script ~/Physics/Sources/FRGscalarMSR/examples/runNotebook.wl ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZNLO/KPZNLO.nb
   Edit the cells here, not in the .nb. The pipeline cells mirror Pipeline.wl. *)

dir = DirectoryName[$InputFileName];
title[s_] := Cell[s, "Title"];
section[s_] := Cell[s, "Section"];
text[s_] := Cell[s, "Text"];
input[s_] := Cell[s, "Input"];

cells = {
   title["Benchmark: KPZ at NLO (functions of p and of the external frequency w)"],
   text["Kloss, Canet, Wschebor, PRE 86, 051124 (2012), Eqs. 43, 44, 49: the functions f_D, f_nu, f_lambda depend on (w, p) on the \
left-hand side; in the loop every function is at zero frequency, and w enters only through G(omega + w, Q). The f_D row is the even part \
in w of d_s Gamma(0,2), f_nu the even part of d_s Gamma(1,1)/(nu p^2), f_lambda the odd part of d_s Gamma(1,1) divided by i w \
(design/NLOPIPELINE.md). With f_lambda = 1 this is KCW's LO, which BMWflow calls \"NLO\"."],
   input["Get[\"~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl\"];\nkpznloDir = ExpandFileName[\"~/Physics/Sources/FRGscalarMSR/benchmarks/KPZNLO\"];\n\
Block[{Print}, Get[FileNameJoin[{ParentDirectory[kpznloDir], #, \"Pipeline.wl\"}]] & /@ {\"KPZLO\", \"KPZLOLambda\"}];"],

   section["1. NLO equations"],
   text["The LO loops of benchmarks/KPZLO and KPZLOLambda, with \"ExternalFrequency\" -> True. The card shows the NLO checks: every row at \
w = 0 equals its LO row; Gamma(hbar, hbar) is even in w; the f_lambda row vanishes at p = 0 and the p^0 order of the f_nu row vanishes, \
for every w."],
   input["kpzlamnloEquations = BuildLOEquations[kpzlamLoops, kpzlamSymmetries, \"ExternalFrequency\" -> True];"],
   input["RunOracle[kpzlamnloEquations];"],
   input["kpzlamnloAdim = AdimensionalEquations[kpzlamnloEquations, kpzlamDims];\nkpzlamnloAdim[\"Linear\"]"],
   input["kpznloEquations = BuildLOEquations[kpzloLoops, kpzloSymmetries, \"ExternalFrequency\" -> True];\nkpznloAdim = AdimensionalEquations[kpznloEquations, kpzloDims];"],

   section["2. References and checks"],
   text["Generic f_lambda and w against BMWflow's NLO kernels with f_lambda slots (Reference.wl); the w = 0 slice against the LO \
benchmarks; the linear part against KCW Eq. 56; the oracle at w != 0."],
   input["Get[FileNameJoin[{kpznloDir, #}]] & /@ {\"Reference.wl\", \"Checks.wl\"};\nkpznloResults = kpznloChecks[];\nkpznloSummary[kpznloResults]"],
   input["If[AllTrue[kpznloResults, #Pass &], \"all checks pass\", Abort[]]"],

   section["3. Export (contract v4)"],
   input["ExportPython[kpzlamnloAdim, \"G\" -> alamk^2, \"Model\" -> \"KPZLambda\", \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZLambda_NLO.py\"}]]"]};

Export[FileNameJoin[{dir, "KPZNLO.nb"}], Notebook[cells, StyleDefinitions -> "Default.nb"]];
