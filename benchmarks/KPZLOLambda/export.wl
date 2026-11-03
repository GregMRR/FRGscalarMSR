(* ::Package:: *)

(* export.wl -- writes the two BMWflow modules compared by runBMWflow.py into out/ (git-ignored):
   out/KPZ_LO.py (f_lambda = 1, benchmarks/KPZLO) and out/KPZLambda_LO.py (flowing f_lambda, Pipeline.wl here), contract v3.
   Run:  math -script ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/export.wl *)

kpzlamDir = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[kpzlamDir, 2], "FRGscalarMSR.wl"}]];
Block[{Print}, Get[FileNameJoin[{ParentDirectory[kpzlamDir], "KPZLO", "Pipeline.wl"}]]; Get[FileNameJoin[{kpzlamDir, "Pipeline.wl"}]]];
ExportPython[kpzloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{kpzlamDir, "out", "KPZ_LO.py"}]];
ExportPython[kpzlamAdim, "G" -> alamk^2, "Model" -> "KPZLambda", "File" -> FileNameJoin[{kpzlamDir, "out", "KPZLambda_LO.py"}]];
