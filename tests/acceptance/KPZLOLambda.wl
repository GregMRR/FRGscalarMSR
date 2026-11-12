(* ::Package:: *)

(* Suite hook of benchmarks/KPZLOLambda (design/LOFLAMBDA.md): KPZ at LO with a flowing f_lambda. Acceptance 1: with f_lambda = 1 the
   equations of fD and fnu equal the KPZLO benchmark's exactly; the same through the Python printer (numpy only, SKIP without it):
   the exported KPZLambda module at f_lambda = 1, f_lambda' = 0 equals the exported KPZ module at random points. The benchmark's
   logic lives in benchmarks/KPZLOLambda; this file wraps its checks in test[...]. *)

kpzlamBenchDir = FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZLOLambda"}];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "PythonRoundTrip.wl"}]];
If[! ValueQ[kpzlamDEEquations], quiet[Get[FileNameJoin[{kpzlamBenchDir, "Pipeline.wl"}]]]];
Get[FileNameJoin[{kpzlamBenchDir, #}]] & /@ {"Reference.wl", "Checks.wl"};

Scan[test["KPZ LO + f_lambda benchmark: " <> #Name <> " [" <> #Paper <> "]", #Pass, True] &, kpzlamChecks[]];

(* the comparison discriminates: a 1/100 change of f_lambda at the external leg in the fnu row is detected *)
test["KPZ LO + f_lambda benchmark: a perturbed fnu row is rejected",
  ExactlyZeroQ[((kpzlamAdim["Loop"][[2]] /. flam[ap^2] -> 101/100 flam[ap^2]) /. kpzlamOne) - kpzloAdim["Loop"][[2]]], False];

pythonRoundTrip["KPZ LO + f_lambda", kpzlamAdim, "KPZLambda", alamk^2];

(* f_lambda = 1 through the printer: KPZLambda_LO.py (fD, fnu rows) = KPZ_LO.py at random points, d in {2, 3/2, 1} *)
If[$pythonExe === None, skip["KPZ LO + f_lambda: Python module at f_lambda = 1 = KPZ module", "python3 with numpy not found", "python"],
  Module[{dir = FileNameJoin[{$pythonDir, "flam1"}], fK, fL, script, out},
   fK = quiet[ExportPython[kpzloAdim, "G" -> alamk^2, "Model" -> "KPZ", "File" -> FileNameJoin[{dir, "KPZ_LO.py"}]]];
   fL = quiet[ExportPython[kpzlamAdim, "G" -> alamk^2, "Model" -> "KPZLambda", "File" -> FileNameJoin[{dir, "KPZLambda_LO.py"}]]];
   script = StringRiffle[{"import importlib.util, numpy as np",
      "def load(f):", "    s = importlib.util.spec_from_file_location('m' + str(abs(hash(f))), f); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m",
      "K = load(" <> ToString[fK, InputForm] <> "); L = load(" <> ToString[fL, InputForm] <> ")",
      "assert L.ETAS == ('etaD', 'etanu', None) and L.FUNCTIONS == ('fD', 'fnu', 'flam') and L.P2_DIVIDED == (0, 1, 0) and L.CONTRACT_VERSION == 4",
      "rng = np.random.default_rng(7); worst = 0.0",
      "for d in (2.0, 1.5, 1.0):",
      "    for _ in range(20):",
      "        p, q = rng.uniform(0.1, 2.0, 2); c = rng.uniform(-1, 1); Q2 = p*p + q*q + 2*p*q*c",
      "        f2 = tuple(rng.uniform(0.5, 2.0, 2)); fq2 = tuple(rng.uniform(0.5, 2.0, 2)); fQ2 = tuple(rng.uniform(0.5, 2.0, 2)); fp2 = tuple(rng.uniform(0.5, 2.0, 2))",
      "        fpr2 = tuple(rng.uniform(-1, 1, 2)); r = tuple(rng.uniform(0.1, 2.0, 2)); r_ = tuple(-rng.uniform(0.1, 2.0, 2)); rQ = tuple(rng.uniform(0.1, 2.0, 2))",
      "        eta = tuple(rng.uniform(-1, 1, 2)); g = (rng.uniform(0.5, 5.0),)",
      "        pairs = [(K.loop(p, q, c, Q2, fp2, fq2, fQ2, r, r_, rQ, eta, g, d), L.loop(p, q, c, Q2, fp2 + (1.0,), fq2 + (1.0,), fQ2 + (1.0,), r, r_, rQ, eta, g, d)[:2]),",
      "                 (K.linear(p, f2, fpr2, eta, g, d), L.linear(p, f2 + (1.0,), fpr2 + (0.0,), eta, g, d)[:2]),",
      "                 (K.eta_vector(q, f2, fpr2, r, r_, g, d), L.eta_vector(q, f2 + (1.0,), fpr2 + (0.0,), r, r_, g, d)),",
      "                 (np.ravel(K.eta_matrix(q, f2, fpr2, r, r_, g, d)), np.ravel(L.eta_matrix(q, f2 + (1.0,), fpr2 + (0.0,), r, r_, g, d))),",
      "                 (K.beta(g, eta, d), L.beta(g, eta, d))]",
      "        for a, b in pairs:",
      "            a = np.array(a, dtype=float); b = np.array(b, dtype=float)",
      "            worst = max(worst, float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-300))))",
      "print(worst)"}, "\n"];
   Export[FileNameJoin[{dir, "flam1.py"}], script, "Text"];
   out = RunProcess[{$pythonExe, "-I", FileNameJoin[{dir, "flam1.py"}]}];
   If[out["ExitCode"] =!= 0, Print["    python stderr: ", StringTake[out["StandardError"], UpTo[800]]]];
   test["KPZ LO + f_lambda: Python module at f_lambda = 1 = KPZ module (loop, linear, eta, beta; rel. 1e-13)",
    out["ExitCode"] === 0 && ToExpression[StringReplace[StringTrim[out["StandardOutput"]], "e" -> "*^"]] < 10^-13, True]]];
