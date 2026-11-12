(* ::Package:: *)

(* KPZ NLO stages shared by unit and acceptance tests: benchmarks/KPZNLO/Pipeline.wl (kpznlo*: f_lambda = 1; kpzlamnlo*: flowing
   f_lambda), on top of the KPZLO and KPZLOLambda pipelines. Built once per run. *)

Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
If[! ValueQ[kpzlamAdim], quiet[Get[FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZLOLambda", "Pipeline.wl"}]]]];
If[! ValueQ[kpzlamnloAdim], quiet[Get[FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZNLO", "Pipeline.wl"}]]]];
