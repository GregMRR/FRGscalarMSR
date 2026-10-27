(* ::Package:: *)

(* PythonExport.wl -- ExportPython: write an LO Adimensional stage as a Python module for BMWflow's generic adapter
   (contract: docs/handoff/BMWflow_LO.md = design/LOPIPELINE.md section 6). Pure numpy + math, elementwise, no state, no
   value-dependent branching. The loop measure is not in the module: the adapter applies it.
   Objective, assumptions and checks: docs/MODULES.md section "PythonExport". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

ExportPython::usage = "ExportPython[adim, \"Model\" -> \"KPZ\", \"File\" -> \".../eqs/KPZ_LO.py\"] writes the LO Adimensional stage as \
a Python module (constants FUNCTIONS, ETAS, COUPLINGS, SHAPES, P2_DIVIDED, APPROXIMATION, MODEL; functions linear, loop, eta_matrix, \
eta_vector, beta; docs/handoff/BMWflow_LO.md) and returns its path. \"Model\" names the row of the BMWflow model table (PythonModelTable).";
PythonParts::usage = "PythonParts[adim] gives the exported expressions with the Python argument names as placeholders \
PythonVariable[\"fq\", i] etc. (Mathematica expressions, as ExportPython prints them): <|\"linear\" -> ..., \"loop\" -> ..., ...|>.";
PythonVariable::usage = "PythonVariable[name] or PythonVariable[name, i] is the Python argument name (name[i], 0-based) in PythonParts.";
PythonModelTable::usage = "PythonModelTable is the BMWflow model table (data): MODEL -> <|\"Coupling\", \"G\", \"BetaG\", \"Shapes\", \
\"ShapeBinding\"|>: how the BMWflow adapter maps the exported coupling to its scalar g, and how it binds the regulator shapes.";

Begin["`Private`"];

SetAttributes[PythonVariable, Protected];

(* the model table of design/LOPIPELINE.md section 6 (data; BMWflow's adapter implements the same rows) *)
PythonModelTable = <|
   "KPZ" -> <|"Coupling" -> "alamk", "G" -> "g = alamk^2", "BetaG" -> "beta_g = 2 alamk beta_alamk", "Shapes" -> {"rn", "rd"},
     "ShapeBinding" -> "BMWflow's single r binds rn and rd"|>|>;

(* ---------- printer: Power -> **, numbers as floats, math.pi; nothing but + - * / ** and np.sqrt ---------- *)
pyNumber[n_Integer] := ToString[n] <> ".0";
pyNumber[r_Rational] := "(" <> pyNumber[Numerator[r]] <> "/" <> pyNumber[Denominator[r]] <> ")";
pyAtomQ[e_] := MatchQ[e, _Symbol | _PythonVariable | _Integer | _Rational | Pi | E];
pyParen[e_] := If[pyAtomQ[e] && ! (MatchQ[e, _Integer | _Rational] && e < 0), py[e], "(" <> py[e] <> ")"];
py[n_Integer] := pyNumber[n];
py[r_Rational] := pyNumber[r];
py[Pi] := "math.pi";
py[E] := "math.e";
py[PythonVariable[s_String]] := s;
py[PythonVariable[s_String, i_Integer]] := s <> "[" <> ToString[i] <> "]";
py[s_Symbol] := StringReplace[SymbolName[s], "$" -> "_"];
py[e_Plus] := StringRiffle[py /@ List @@ e, " + "];
py[e_Times] := Module[{fs = List @@ e, num, den},
   num = Select[fs, ! MatchQ[#, Power[_, _?Negative]] &]; den = Cases[fs, Power[b_, x_?Negative] :> Power[b, -x]];
   If[den === {}, StringRiffle[pyParen /@ num, "*"],
    If[num === {}, "1.0", StringRiffle[pyParen /@ num, "*"]] <> "/(" <> StringRiffle[pyParen /@ den, "*"] <> ")"]];
py[Power[b_, 1/2]] := "np.sqrt(" <> py[b] <> ")";
py[Power[b_, -1/2]] := "1.0/np.sqrt(" <> py[b] <> ")";
py[Power[b_, n_Integer]] := If[n < 0, "1.0/" <> pyParen[Power[b, -n]], pyParen[b] <> "**" <> ToString[n]];
py[Power[b_, x_]] := pyParen[b] <> "**(" <> py[x] <> ")";
py[e_] := FRGAssert["PythonExport.printable", False, e, "internal: no Python translation for this expression"];

(* ---------- the exported expressions with Python placeholders ---------- *)

(* functions f at aq^2 / at Qhat^2, shapes at aq^2 (and r' = dr/d(q^2)) / at Qhat^2; F' = d/d(p^2) -> f'/(2 momentum), f' = d/d|momentum| *)
pyLoopRules[fs_, shapes_, aq_, ap_, c_] := With[{q2 = aq^2, Q2 = Expand[ap^2 + aq^2 + 2 ap aq c]}, Join[
    Flatten[MapIndexed[{#1[q2] -> PythonVariable["fq", #2[[1]] - 1], #1[Q2] -> PythonVariable["fQ", #2[[1]] - 1]} &, fs]],
    Flatten[MapIndexed[{#1[q2] -> PythonVariable["rq", #2[[1]] - 1], Derivative[1][#1][q2] -> PythonVariable["rq_", #2[[1]] - 1],
        #1[Q2] -> PythonVariable["rQ", #2[[1]] - 1]} &, shapes]],
    {ap -> PythonVariable["p"], aq -> PythonVariable["q"], c -> PythonVariable["c"]}]];
pyEtaRules[fs_, shapes_, aq_] := Join[
   Flatten[MapIndexed[{#1[aq^2] -> PythonVariable["f", #2[[1]] - 1], Derivative[1][#1][aq^2] -> PythonVariable["fprime", #2[[1]] - 1]/(2 PythonVariable["q"])} &, fs]],
   Flatten[MapIndexed[{#1[aq^2] -> PythonVariable["r", #2[[1]] - 1], Derivative[1][#1][aq^2] -> PythonVariable["r_", #2[[1]] - 1]} &, shapes]],
   {aq -> PythonVariable["q"]}];
pyLinearRules[fs_, ap_] := Join[
   Flatten[MapIndexed[{#1[ap^2] -> PythonVariable["f", #2[[1]] - 1], Derivative[1][#1][ap^2] -> PythonVariable["fprime", #2[[1]] - 1]/(2 PythonVariable["p"])} &, fs]],
   {ap -> PythonVariable["p"]}];

PythonParts[adim_Association] := Block[{$FRGModule = "PythonExport"},
  Module[{fs = adim["Functions"], shapes = adim["Shapes"], aq = adim["LoopMomentum"], ap = adim["ExternalMomentum"], c = adim["Cosine"], common},
   common = Join[Thread[adim["Etas"] -> (PythonVariable["eta", #] & /@ Range[0, Length[adim["Etas"]] - 1])],
     Thread[adim["Couplings"] -> (PythonVariable["couplings", #] & /@ Range[0, Length[adim["Couplings"]] - 1])],
     If[adim["Dimension"] === 2, {}, {adim["Dimension"] -> PythonVariable["d"]}]];
   <|"linear" -> Expand[adim["Linear"] /. pyLinearRules[fs, ap] /. common],
    "loop" -> (adim["Loop"] /. pyLoopRules[fs, shapes, aq, ap, c] /. common),
    "eta_matrix" -> (adim["EtaMatrix"] /. pyEtaRules[fs, shapes, aq] /. common),
    "eta_vector" -> (adim["EtaVector"] /. pyEtaRules[fs, shapes, aq] /. common),
    "beta" -> (adim["BetaL"] /. common)|>]];

pyTuple[l_List] := "(" <> StringRiffle[l, ", "] <> If[Length[l] == 1, ",", ""] <> ")";
pyStrings[l_List] := pyTuple[("\"" <> ToString[#] <> "\"") & /@ l];
pyDef[name_, args_, body_] := "def " <> name <> "(" <> StringRiffle[args, ", "] <> "):\n    return " <> body <> "\n";
pyRows[l_List] := "(\n        " <> StringRiffle[If[# === 0, "0.0", py[#]] & /@ l, ",\n        "] <> If[Length[l] == 1, ",", ""] <> ")";

Options[ExportPython] = {"Model" -> None, "File" -> None};

ExportPython[adim_, opts : OptionsPattern[]] := Block[{$FRGModule = "PythonExport"},
  Module[{model = OptionValue["Model"], file = OptionValue["File"], parts, checks = {}, warns = {}, row, all, fs, cps, shapes, etaNames, bad, units, header,
    defs, allowed},
   FRGAssert["PythonExport.pipeline", AssociationQ[adim] && Lookup[adim, "Pipeline", None] === "LO", If[AssociationQ[adim], Lookup[adim, "Pipeline", None], adim],
    "ExportPython writes LO stages (functions of p^2) for BMWflow; a derivative-expansion stage goes to ExportJulia (FlowRG)"];
   validateStage[adim, "Adimensional"];
   FRGAssert["PythonExport.options", MatchQ[{opts}, {(Rule | RuleDelayed)["Model" | "File", _] ...}], {opts}, "options: \"Model\" -> name, \"File\" -> path"];
   FRGAssert["PythonExport.file", StringQ[file] && StringEndsQ[file, ".py"], file, "\"File\" -> \".../eqs/<Model>_LO.py\""];
   fs = adim["Functions"]; cps = adim["Couplings"]; shapes = adim["Shapes"]; etaNames = adim["EtaNames"];

   AppendTo[checks, assertCheck["PythonExport.etaPerFunction", Length[adim["Etas"]] == Length[fs] &&
       adim["Etas"] === (adim["Dimensions"]["Etas"][#] & /@ (fs /. adim["FunctionReferences"])), {fs, adim["Etas"]},
     "one eta per function, ETAS in FUNCTIONS order (eta_X belongs to the reference of FUNCTIONS[X])",
     "ETAS " <> pyStrings["eta" <> # & /@ etaNames] <> " <-> FUNCTIONS " <> pyStrings[fs]]];
   AppendTo[checks, assertCheck["PythonExport.singleCoupling", Length[cps] == 1, cps,
     "BMWflow's evolution state has one scalar g: exactly one exported coupling (n_c > 1 is in docs/FUTURE.md)", "one coupling: " <> ToString[cps]]];
   AppendTo[checks, assertCheck["PythonExport.p2DividedAtMostOne", Max[adim["P2Divided"]] <= 1, Thread[fs -> adim["P2Divided"]],
     "a row divided by p^4 or more: the theta-quadrature error divided by p1^4 breaks BMWflow's tolerance (docs/FUTURE.md)",
     "P2_DIVIDED " <> ToString[adim["P2Divided"]] <> " <= 1"]];
   row = Lookup[PythonModelTable, model, Missing[]];
   AppendTo[checks, assertCheck["PythonExport.modelTable", ! MissingQ[row] && {row["Coupling"]} === (ToString /@ cps) && Sort[row["Shapes"]] === Sort[ToString /@ shapes],
     <|"model" -> model, "row" -> row, "COUPLINGS" -> cps, "SHAPES" -> shapes|>,
     "\"Model\" must name a row of PythonModelTable whose coupling is the exported one and whose shapes are the regulator shapes (as a set)",
     "MODEL " <> ToString[model] <> ": " <> If[MissingQ[row], "", row["G"] <> ", " <> row["ShapeBinding"]]]];

   parts = PythonParts[adim];
   all = Values[parts];
   AppendTo[checks, assertCheck["PythonExport.noMeasure", FreeQ[all, Gamma | Sin | Power[Pi | 2 | 2 Pi, x_ /; ! FreeQ[x, PythonVariable["d"]]]],
     Cases[all, _Gamma | _Sin | Power[Pi | 2 | 2 Pi, x_ /; ! FreeQ[x, PythonVariable["d"]]], {0, Infinity}],
     "internal: the loop measure (v_d, v_{d-1}, Gamma(d), the (2 Pi)^d) must not be in the module: the adapter applies it",
     "no v_d, v_{d-1}, Gamma(d), (2 Pi)^d in the kernels (the adapter applies the measure)"]];
   units = Join[adim["Dimensions"]["Basis"], {adim["Equations"]["LoopMomentum"], adim["Equations"]["ExternalMomentum"], adim["LoopMomentum"], adim["ExternalMomentum"],
      adim["Cosine"]}, If[adim["Dimension"] === 2, {}, {adim["Dimension"]}]];
   bad = Union[Cases[all /. _PythonVariable -> 1, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> True]];
   AppendTo[checks, assertCheck["PythonExport.knownSymbols", bad === {}, bad,
     "a symbol is left after the Python argument substitution: a unit, an unexpected function atom (F' in loop, r' at Q, r''), or a coupling outside COUPLINGS",
     "only the contract arguments (p, q, c, fq, fQ, rq, rq_, rQ, f, fprime, r, r_, eta, couplings, d)"]];
   bad = Complement[Union[Cases[all /. _PythonVariable -> pyPlaceholder, h_Symbol[___] /; Context[h] === "System`" :> h, {0, Infinity}, Heads -> True]],
     {Plus, Times, Power, Rational, List}];
   AppendTo[checks, assertCheck["PythonExport.systemHeads", bad === {}, bad, "a System function without a Python translation is left", "only + * ** (np.sqrt)"]];

   header = StringRiffle[{
      "# Generated by FRGscalarMSR ExportPython on " <> DateString[{"Year", "-", "Month", "-", "Day", " ", "Time"}] <> ". Do not hand-edit: regenerate.",
      "# Source: Adimensional stage " <> StringTake[adim["Fingerprint"], 8] <> " (LO), model " <> ToString[model] <> ".",
      "# Contract: FRGscalarMSR docs/handoff/BMWflow_LO.md. Pure numpy + math, elementwise, no state.",
      "# Dimensionless; p, q in units of the RG scale; f_X(0) = 1; fprime = d f/d|momentum|; r_ = dr/d(q^2); every rate is d_s x, s = ln(kappa/Lambda).",
      "# The loop measure is NOT included: I_X(p) = v_{d-1}/(2 pi) Int dq q^(d-1) Int dtheta sin^(d-2)(theta) loop_X; eta: <.> = v_d Int dq q^(d-1) (.).",
      "# loop_X is already divided by p^(2 P2_DIVIDED[X]) (non-finite at p = 0 for those rows: the adapter extrapolates).",
      "# Model table row " <> ToString[model] <> ": " <> If[MissingQ[row], "", row["G"] <> ", " <> row["BetaG"] <> "; " <> row["ShapeBinding"] <> "."]}, "\n"];
   defs = {
     pyDef["linear", {"p", "f", "fprime", "eta", "couplings", "d"}, pyRows[parts["linear"]]],
     pyDef["loop", {"p", "q", "c", "fq", "fQ", "rq", "rq_", "rQ", "eta", "couplings", "d"}, pyRows[parts["loop"]]],
     pyDef["eta_matrix", {"q", "f", "fprime", "r", "r_", "couplings", "d"}, "(\n        " <> StringRiffle[pyRows /@ parts["eta_matrix"], ",\n        "] <>
       If[Length[parts["eta_matrix"]] == 1, ",", ""] <> ")"],
     pyDef["eta_vector", {"q", "f", "fprime", "r", "r_", "couplings", "d"}, pyRows[parts["eta_vector"]]],
     pyDef["beta", {"couplings", "eta", "d"}, pyRows[parts["beta"]]]};
   Quiet[CreateDirectory[DirectoryName[file]], CreateDirectory::eexist];
   Export[file, StringRiffle[Join[{header, "", "import math", "", "import numpy as np", "",
       "FUNCTIONS = " <> pyStrings[fs], "ETAS = " <> pyStrings["eta" <> # & /@ etaNames], "COUPLINGS = " <> pyStrings[cps], "SHAPES = " <> pyStrings[shapes],
       "P2_DIVIDED = " <> pyTuple[ToString /@ adim["P2Divided"]], "APPROXIMATION = \"LO\"", "MODEL = \"" <> ToString[model] <> "\"", "", ""},
      Riffle[defs, "\n"]], "\n"], "Text"];
   AppendTo[warns, infoCheck["INFO", "PythonExport.adapter", "the adapter applies the measure, the p = 0 extrapolation of P2_DIVIDED rows and the model-table map " <>
      If[MissingQ[row], "", row["G"]]]];
   Print["ExportPython: wrote ", file, " (", Length[fs], " functions, ", Length[cps], " coupling, model ", model, ")"];
   Scan[Print[StringPadRight[#Status, 11], #Check, If[#Detail === "", "", "   " <> ToString[#Detail]]] &, Join[checks, warns]];
   file]];

End[];
EndPackage[];
