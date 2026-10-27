(* ::Package:: *)

(* Python round trip of an exported LO module (design/LOPIPELINE.md section 7). pythonRoundTrip[name, adim, model] exports adim,
   imports the module in python3 (numpy), evaluates every function (linear, loop, eta_matrix, eta_vector, beta) at a few points,
   including a non-integer d, and compares with the same parts evaluated exactly in Mathematica (relative 1e-13). The point
   values are dyadic rationals, exact in binary floating point. A second call passes numpy arrays: loop must broadcast
   elementwise. Without python3 + numpy the checks are reported SKIP (a failure with --require-python). *)

$pythonExe = Quiet[Check[If[RunProcess[{"python3", "-c", "import numpy"}, "ExitCode"] === 0, "python3", None], None]];
$pythonDir = FileNameJoin[{$TemporaryDirectory, "FRGscalarMSR-python-" <> ToString[$ProcessID]}];
Quiet[CreateDirectory[$pythonDir], CreateDirectory::eexist];

fromPython[s_String] := ToExpression[StringReplace[s, {"e" -> "*^", "nan" -> "Indeterminate", "inf" -> "Infinity", "[" -> "{", "]" -> "}"}]];
pyLiteral[x_] := ToString[N[x], CForm];
pyTupleOf[l_List] := "(" <> StringRiffle[pyLiteral /@ l, ", "] <> If[Length[l] == 1, ",", ""] <> ")";

(* one random dyadic point: every argument of the contract *)
pythonPoint[adim_, seed_, d_] := Module[{nf = Length[adim["Functions"]], ns = Length[adim["Shapes"]], ne = Length[adim["Etas"]], nc = Length[adim["Couplings"]], dy, tup},
   SeedRandom[seed];
   dy[] := RandomInteger[{16, 127}]/128;
   tup[n_] := Table[dy[], n];
   <|"p" -> dy[], "q" -> dy[], "c" -> RandomInteger[{-120, 120}]/128, "fq" -> tup[nf], "fQ" -> tup[nf], "rq" -> tup[ns], "rq_" -> -tup[ns], "rQ" -> tup[ns],
    "f" -> tup[nf], "fprime" -> -tup[nf], "r" -> tup[ns], "r_" -> -tup[ns], "eta" -> tup[ne], "couplings" -> tup[nc], "d" -> d|>];

pythonValue[e_, pt_] := N[e /. {PythonVariable[s_, i_] :> pt[s][[i + 1]], PythonVariable[s_] :> pt[s]}, 30];

pythonArgs = <|"linear" -> {"p", "f", "fprime", "eta", "couplings", "d"}, "loop" -> {"p", "q", "c", "fq", "fQ", "rq", "rq_", "rQ", "eta", "couplings", "d"},
   "eta_matrix" -> {"q", "f", "fprime", "r", "r_", "couplings", "d"}, "eta_vector" -> {"q", "f", "fprime", "r", "r_", "couplings", "d"},
   "beta" -> {"couplings", "eta", "d"}|>;

pythonRoundTrip[name_String, adim_, model_String] := Module[{file, parts, pts, mma, script, out, py, rel, call, bro},
   If[$pythonExe === None, skip[name <> ": Python round trip", "python3 with numpy not found", "python"]; Return[Null]];
   file = quiet[ExportPython[adim, "Model" -> model, "File" -> FileNameJoin[{$pythonDir, name, model <> "_LO.py"}]]];
   parts = PythonParts[adim];
   pts = {pythonPoint[adim, 1, 2], pythonPoint[adim, 2, 3/2], pythonPoint[adim, 3, 1], pythonPoint[adim, 4, 5/2]};
   mma = Flatten[Table[pythonValue[parts[fn], pt], {pt, pts}, {fn, Keys[pythonArgs]}]];
   call[fn_, pt_] := "M." <> fn <> "(" <> StringRiffle[If[ListQ[pt[#]], pyTupleOf[pt[#]], pyLiteral[pt[#]]] & /@ pythonArgs[fn], ", "] <> ")";
   (* broadcasting: loop with p an array of two values = the two scalar calls *)
   bro = With[{pt = First[pts], p2 = 37/128}, "np.allclose(np.array(M.loop(np.array([" <> pyLiteral[pt["p"]] <> ", " <> pyLiteral[p2] <> "]), " <>
       StringRiffle[If[ListQ[pt[#]], pyTupleOf[pt[#]], pyLiteral[pt[#]]] & /@ Rest[pythonArgs["loop"]], ", "] <> ")), np.array([" <>
       call["loop", pt] <> ", " <> call["loop", Append[pt, "p" -> p2]] <> "]).T, rtol=1e-15, atol=0)"];
   script = StringRiffle[{"import importlib.util, json, numpy as np",
      "spec = importlib.util.spec_from_file_location('m', " <> ToString[file, InputForm] <> ")",
      "M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)",
      "flat = lambda x: [float(v) for v in np.ravel(np.array(x, dtype=float))]",
      "v = []",
      Sequence @@ Flatten[Table["v += flat(" <> call[fn, pt] <> ")", {pt, pts}, {fn, Keys[pythonArgs]}]],
      "print(json.dumps(v))",
      "print(bool(" <> bro <> "))"}, "\n"];
   Export[FileNameJoin[{$pythonDir, name, "roundtrip.py"}], script, "Text"];
   out = RunProcess[{$pythonExe, "-I", FileNameJoin[{$pythonDir, name, "roundtrip.py"}]}];
   test[name <> ": python imports and evaluates the exported module", out["ExitCode"], 0];
   If[out["ExitCode"] =!= 0, Print["    python stderr: ", StringTake[out["StandardError"], UpTo[800]]]; Return[Null]];
   {py, bro} = StringSplit[StringTrim[out["StandardOutput"]], "\n"];
   py = fromPython[py];
   rel = Max[MapThread[If[#2 == 0, Abs[#1], Abs[(#1 - #2)/#2]] &, {py, mma}]];
   test[name <> ": Python = Mathematica for every function at 4 points, d in {2, 3/2, 1, 5/2} (rel. 1e-13)", Length[py] == Length[mma] && rel < 10^-13, True];
   test[name <> ": loop broadcasts over a numpy array of p", bro, "True"]];
