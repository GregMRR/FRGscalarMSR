(* ::Package:: *)

(* Cross-repo guard of the LO handoff (docs/handoff/BMWflow_LO.md, design/LOPIPELINE.md section 7): export the KPZ LO stage and run
   BMWflow's own acceptance test on that file (ModelGenerated on it vs the hand-written ModelKPZ: integrals, eta, beta, a 2000-step
   run). BMWflow is at BMWFLOW_PATH or ~/Physics/Sources/BMWflow, with a venv in .venv (scipy, numba, pytest). Without them the test
   is SKIP (a failure with --require-python). *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];

bmwRoot = Replace[Environment["BMWFLOW_PATH"], $Failed -> FileNameJoin[{$HomeDirectory, "Physics", "Sources", "BMWflow"}]];
bmwPython = FileNameJoin[{bmwRoot, ".venv", "bin", "python"}];
bmwTest = FileNameJoin[{bmwRoot, "tests", "test_model_generated.py"}];

If[! (FileExistsQ[bmwPython] && FileExistsQ[bmwTest] && RunProcess[{bmwPython, "-c", "import pytest, scipy, numba"}, "ExitCode"] === 0),
  skip["BMWflow acceptance on a fresh KPZ LO export", "BMWflow with .venv (pytest, scipy, numba) not found at " <> bmwRoot, "python"],
  Module[{dir = FileNameJoin[{$TemporaryDirectory, "FRGscalarMSR-bmwflow-" <> ToString[$ProcessID]}], file, out},
   file = quiet[ExportPython[kpzloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{dir, "KPZ_LO.py"}]]];
   out = RunProcess[{bmwPython, "-m", "pytest", "-q", "-p", "no:cacheprovider", bmwTest, "--module", file}, ProcessDirectory -> bmwRoot];
   If[out["ExitCode"] =!= 0, Print["    pytest: ", StringTake[out["StandardOutput"], -UpTo[1500]]]];
   test["BMWflow acceptance on a fresh KPZ LO export (ModelGenerated = ModelKPZ, d = 1, 2)", out["ExitCode"], 0];
   test["the fresh export equals BMWflow's vendored tests/fixtures/KPZ_LO.py (else: copy it there)",
    Import[file, "Text"] === Import[FileNameJoin[{bmwRoot, "tests", "fixtures", "KPZ_LO.py"}], "Text"], True]]];
