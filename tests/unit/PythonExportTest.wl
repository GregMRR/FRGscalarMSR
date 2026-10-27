(* ::Package:: *)

(* PythonExport (design/LOPIPELINE.md step 8): the KPZ LO module, its constants, the Python round trip (SKIP without python3 +
   numpy), and every check (the contract checks on stages forged from the KPZ one, as tests/acceptance/KPZAdimensional.wl does). *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "PythonRoundTrip.wl"}]];

pyDir = FileNameJoin[{$pythonDir, "unit"}];
pyFile = quiet[ExportPython[kpzloAdim, "Model" -> "KPZ", "File" -> FileNameJoin[{pyDir, "KPZ_LO.py"}]]];
pyText = Import[pyFile, "Text"];
test["KPZ module: contract constants", AllTrue[{"FUNCTIONS = (\"fD\", \"fnu\")", "ETAS = (\"etaD\", \"etanu\")", "COUPLINGS = (\"alamk\",)",
    "SHAPES = (\"rn\", \"rd\")", "P2_DIVIDED = (0, 1)", "APPROXIMATION = \"LO\"", "MODEL = \"KPZ\""}, StringContainsQ[pyText, #] &], True];
test["KPZ module: contract signatures", AllTrue[{"def linear(p, f, fprime, eta, couplings, d):", "def loop(p, q, c, fq, fQ, rq, rq_, rQ, eta, couplings, d):",
    "def eta_matrix(q, f, fprime, r, r_, couplings, d):", "def eta_vector(q, f, fprime, r, r_, couplings, d):", "def beta(couplings, eta, d):"},
   StringContainsQ[pyText, #] &], True];
test["KPZ module: Python syntax only (no ^, no Mathematica brackets)", StringFreeQ[StringRiffle[Select[StringSplit[pyText, "\n"], ! StringStartsQ[#, "#"] &], "\n"], "^" | "Derivative" | "Power[" | "Sqrt["], True];
test["KPZ module: linear = eta f + p f' and beta of alamk", {Expand[PythonParts[kpzloAdim]["linear"]], Expand[PythonParts[kpzloAdim]["beta"]]},
  {Expand[{PythonVariable["eta", 0] PythonVariable["f", 0] + PythonVariable["p"] PythonVariable["fprime", 0],
     PythonVariable["eta", 1] PythonVariable["f", 1] + PythonVariable["p"] PythonVariable["fprime", 1]}],
   Expand[{PythonVariable["couplings", 0] (PythonVariable["d"]/2 - 1 + 3 PythonVariable["eta", 1]/2 - PythonVariable["eta", 0]/2)}]}];
test["printer: powers, rationals, sqrt", FRGscalarMSR`Private`py[x^2 + 3/2 y^-1 + Sqrt[z] - 2 Pi], "(-2.0)*math.pi + x**2 + (3.0/2.0)/(y) + np.sqrt(z)"];

pythonRoundTrip["KPZ LO", kpzloAdim, "KPZ"];

(* the checks, each on a stage that breaks it *)
pyForge[changes_] := FRGscalarMSR`Internal`stampStage[KeyDrop[Join[kpzloAdim, changes], "Fingerprint"]];
pyOut = FileNameJoin[{pyDir, "bad.py"}];
expectFail["a DE stage is refused", ExportPython[quiet[AdimensionalEquations[kpzEquations, quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]]]],
   "Model" -> "KPZ", "File" -> pyOut], "PythonExport.pipeline"];
expectFail["ETAS not in FUNCTIONS order", ExportPython[pyForge[<|"Etas" -> {etanu, etaD}|>], "Model" -> "KPZ", "File" -> pyOut], "PythonExport.etaPerFunction"];
expectFail["two couplings", ExportPython[pyForge[<|"Couplings" -> {alamk, blamk}, "BetaL" -> {0, 0}|>], "Model" -> "KPZ", "File" -> pyOut], "PythonExport.singleCoupling"];
expectFail["a row divided by p^4", ExportPython[pyForge[<|"P2Divided" -> {0, 2}|>], "Model" -> "KPZ", "File" -> pyOut], "PythonExport.p2DividedAtMostOne"];
expectFail["an unknown model", ExportPython[kpzloAdim, "Model" -> "Burgers", "File" -> pyOut], "PythonExport.modelTable"];
pyRenamed = quiet[Module[{eq = loEquationsOf[kpzloTheory, kpzloDensity /. lamk -> lamz, {nuk, lamz, Dk}, {fD -> Dk, fnu -> nuk},
       kpzloTags /. lamk -> lamz, kpzloRegs]}, AdimensionalEquations[eq, DeriveDimensions[AnsatzOf[eq], "References" -> {nuk -> "nu", Dk -> "D"}]]]];
expectFail["KPZ with a renamed coupling (alamz) does not match the KPZ row", ExportPython[pyRenamed, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.modelTable"];
expectFail["a loop measure in a kernel", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] Gamma[Dim/2]|>], "Model" -> "KPZ", "File" -> pyOut], "PythonExport.noMeasure"];
expectFail["a unit left in a kernel", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] k|>], "Model" -> "KPZ", "File" -> pyOut], "PythonExport.knownSymbols"];
expectFail["F' in the loop (not in the contract)", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] fnu'[aq^2]|>], "Model" -> "KPZ", "File" -> pyOut],
  "PythonExport.knownSymbols"];
expectFail["a System function without translation", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] Log[aq]|>], "Model" -> "KPZ", "File" -> pyOut],
  "PythonExport.systemHeads"];
