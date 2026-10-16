(* ::Package:: *)

(* Stages: save / load round trip, fingerprint and staleness checks *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

stageFile = FileNameJoin[{$TemporaryDirectory, "FRGscalarMSR-test-" <> ToString[$ProcessID] <> ".wdx"}];
quiet[SaveStage[kpzClosure, stageFile]];
loaded = quiet[LoadStage[stageFile, "Expect" -> {kpzTheory, kpzAnsatz, kpzRegulator, kpzSymmetries}]];
test["round trip: same stage", loaded === kpzClosure, True];
test["round trip: equations usable downstream", Length[quiet[AssembleEquations[loaded]]["Equations"]], 2];

(* a different ansatz (an extra coupling) makes the saved stage stale *)
otherAnsatz = quiet[DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2 + mk hbar h, "Couplings" -> {nuk, lamk, Dk, mk}]];
expectFail["stale LoadStage", LoadStage[stageFile, "Expect" -> {otherAnsatz}], "Stages.stale"];
Export[stageFile, Append[kpzClosure, "Equations" -> {}], "WDX"];
expectFail["modified file", LoadStage[stageFile], "Stage.fingerprint"];
expectFail["missing file", LoadStage[stageFile <> ".none"], "Stages.file"];
Quiet[DeleteFile[stageFile]];
