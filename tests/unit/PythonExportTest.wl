(* ::Package:: *)

(* PythonExport (design/LOPIPELINE.md step 8): the KPZ LO module, its constants, the Python round trip (SKIP without python3 +
   numpy), and every check (the contract checks on stages forged from the KPZ one, as tests/acceptance/KPZAdimensional.wl does). *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "PythonRoundTrip.wl"}]];

pyDir = FileNameJoin[{$pythonDir, "unit"}];
pyFile = quiet[ExportPython[kpzloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{pyDir, "KPZ_LO.py"}]]];
pyText = Import[pyFile, "Text"];
test["KPZ module: contract constants", AllTrue[{"CONTRACT_VERSION = 2", "FUNCTIONS = (\"fD\", \"fnu\")", "ETAS = (\"etaD\", \"etanu\")", "COUPLINGS = (\"g\",)",
    "SHAPES = (\"rn\", \"rd\")", "P2_DIVIDED = (0, 1)", "FIXED_DIMENSION = None", "APPROXIMATION = \"LO\"", "MODEL = \"KPZ\""}, StringContainsQ[pyText, #] &], True];
test["KPZ module: exporting twice gives the same file (no date in the header)",
  Import[quiet[ExportPython[kpzloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{pyDir, "KPZ_LO_again.py"}]]], "Text"] === pyText, True];
test["KPZ module: contract signatures", AllTrue[{"def linear(p, f, fprime, eta, couplings, d):", "def loop(p, q, c, Q2, fq, fQ, rq, rq_, rQ, eta, couplings, d):",
    "def eta_matrix(q, f, fprime, r, r_, couplings, d):", "def eta_vector(q, f, fprime, r, r_, couplings, d):", "def beta(couplings, eta, d):"},
   StringContainsQ[pyText, #] &], True];
test["KPZ module: Python syntax only (no ^, no Mathematica brackets)", StringFreeQ[StringRiffle[Select[StringSplit[pyText, "\n"], ! StringStartsQ[#, "#"] &], "\n"], "^" | "Derivative" | "Power[" | "Sqrt["], True];
test["KPZ module: linear = eta f + p f' and beta of alamk", {Expand[PythonParts[kpzloAdim]["linear"]], Expand[PythonParts[kpzloAdim]["beta"]]},
  {Expand[{PythonVariable["eta", 0] PythonVariable["f", 0] + PythonVariable["p"] PythonVariable["fprime", 0],
     PythonVariable["eta", 1] PythonVariable["f", 1] + PythonVariable["p"] PythonVariable["fprime", 1]}],
   Expand[{PythonVariable["couplings", 0] (PythonVariable["d"]/2 - 1 + 3 PythonVariable["eta", 1]/2 - PythonVariable["eta", 0]/2)}]}];
test["KPZ module in g = alamk^2: beta_g = g (d - 2 - eta_D + 3 eta_nu), BMWflow's g_rhs", Expand[PythonParts[kpzloAdim, alamk^2]["beta"]],
  Expand[{PythonVariable["couplings", 0] (PythonVariable["d"] - 2 + 3 PythonVariable["eta", 1] - PythonVariable["eta", 0])}]];
test["KPZ module in g: the kernels are the alamk kernels at alamk^2 = g", ExactlyZeroQ[
   (Values[KeyDrop[PythonParts[kpzloAdim, alamk^2], "beta"]] /. PythonVariable["couplings", 0] -> PythonVariable["couplings", 0]^2) -
    Values[KeyDrop[PythonParts[kpzloAdim], "beta"]]], True];
test["printer: a row free of the arrays broadcasts", FRGscalarMSR`Private`pyRowString[#, "q"] & /@ {0, 3/2, PythonVariable["eta", 0], PythonVariable["q"]},
  {"0.0*q", "(3.0/2.0) + 0.0*q", "eta[0] + 0.0*q", "q"}];
test["printer: powers, rationals, sqrt", FRGscalarMSR`Private`py[x^2 + 3/2 y^-1 + Sqrt[z] - 2 Pi], "(-2.0)*math.pi + x**2 + (3.0/2.0)/(y) + np.sqrt(z)"];

pythonRoundTrip["KPZ LO", kpzloAdim, "KPZ", alamk^2];

(* the checks, each on a stage that breaks it *)
pyForge[changes_] := FRGscalarMSR`Internal`stampStage[KeyDrop[Join[kpzloAdim, changes], "Fingerprint"]];
pyOut = FileNameJoin[{pyDir, "bad.py"}];
expectFail["a DE stage is refused", ExportPython[quiet[AdimensionalEquations[kpzEquations, quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]]]],
   "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.pipeline"];
expectFail["ETAS not in FUNCTIONS order", ExportPython[pyForge[<|"Etas" -> {etanu, etaD}|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.etaPerFunction"];
expectFail["two couplings", ExportPython[pyForge[<|"Couplings" -> {alamk, blamk}, "BetaL" -> {0, 0}|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.singleCoupling"];
expectFail["a row divided by p^4", ExportPython[pyForge[<|"P2Divided" -> {0, 2}|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.p2DividedAtMostOne"];
expectFail["no \"G\"", ExportPython[kpzloAdim, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.gMap"];
expectFail["\"G\" with a symbol other than the coupling", ExportPython[kpzloAdim, "G" -> alamk^2 Dk, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.gMap"];
expectFail["\"G\" -> alamk^3 leaves g^(2/3) in the kernels", ExportPython[kpzloAdim, "G" -> alamk^3, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.gEven"];
expectFail["no \"Model\" label", ExportPython[kpzloAdim, "G" -> alamk^2, "File" -> pyOut], "PythonExport.model"];
pyRenamed = quiet[Module[{eq = loEquationsOf[kpzloTheory, kpzloDensity /. lamk -> lamz, {nuk, lamz, Dk}, {fD -> Dk, fnu -> nuk},
       kpzloTags /. lamk -> lamz, kpzloRegs]}, AdimensionalEquations[eq, DeriveDimensions[AnsatzOf[eq], "References" -> {nuk -> "nu", Dk -> "D"}]]]];
expectFail["KPZ with a renamed coupling (alamz): \"G\" -> alamk^2 names another coupling", ExportPython[pyRenamed, "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.gMap"];
expectFail["a loop measure in a kernel", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] Gamma[Dim/2]|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.noMeasure"];
expectFail["a unit left in a kernel", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] k|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut], "PythonExport.knownSymbols"];
expectFail["F' in the loop (not in the contract)", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] fnu'[aq^2]|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut],
  "PythonExport.knownSymbols"];
expectFail["a System function without translation", ExportPython[pyForge[<|"Loop" -> kpzloAdim["Loop"] Log[aq]|>], "G" -> alamk^2, "Model" -> "KPZ", "File" -> pyOut],
  "PythonExport.systemHeads"];
