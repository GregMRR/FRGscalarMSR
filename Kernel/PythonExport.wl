(* ::Package:: *)

(* PythonExport.wl -- ExportPython: write an LO Adimensional stage as a Python module for BMWflow's generic adapter
   (contract: docs/handoff/BMWflow_LO.md = design/LOPIPELINE.md section 6). Pure numpy + math, elementwise, no state, no
   value-dependent branching. The loop measure is not in the module: the adapter applies it.
   Objective, assumptions and checks: docs/MODULES.md section "PythonExport". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

ExportPython::usage = "ExportPython[adim, \"G\" -> alamk^2, \"Model\" -> \"KPZ\", \"File\" -> \".../eqs/KPZ_LO.py\"] writes the LO \
Adimensional stage as a Python module (constants CONTRACT_VERSION, FUNCTIONS, ETAS, COUPLINGS, SHAPES, P2_DIVIDED, FIXED_DIMENSION, \
APPROXIMATION, MODEL; functions linear, loop, eta_matrix, eta_vector, beta; docs/handoff/BMWflow_LO.md) and returns its path. \"G\" is \
BMWflow's scalar coupling g in terms of the single exported coupling: every kernel is written in g (COUPLINGS = (\"g\",)) and beta is \
d_s g. \"Model\" is a label (MODEL).";
PythonParts::usage = "PythonParts[adim] gives the exported expressions with the Python argument names as placeholders \
PythonVariable[\"fq\", i] etc. (Mathematica expressions, as ExportPython prints them): <|\"linear\" -> ..., \"loop\" -> ..., ...|>. \
PythonParts[adim, G] writes them in g = G: PythonVariable[\"couplings\", 0] is then g and \"beta\" is {d_s g}.";
PythonAtom::usage = "PythonAtom[\"fQ\", i, \"rQ\", j] in PythonParts is the propagator atom fQ[i] + rQ[j] (printed as that sum): function i and the \
regulator shape j of its Gamma2 entry. The loop is written in atoms so that large r(q), r(Q) do not cancel numerically.";
PythonVariable::usage = "PythonVariable[name] or PythonVariable[name, i] is the Python argument name (name[i], 0-based) in PythonParts.";

Begin["`Private`"];

SetAttributes[{PythonVariable, PythonAtom}, Protected];

(* version of the contract docs/handoff/BMWflow_LO.md; BMWflow's ModelGenerated refuses any other *)
pyContractVersion = 2;
pyHashPlaceholder = "@@HASH@@";

(* ---------- printer: Power -> **, numbers as floats, math.pi; nothing but + - * / ** and np.sqrt ---------- *)
pyNumber[n_Integer] := ToString[n] <> ".0";
pyNumber[r_Rational] := "(" <> pyNumber[Numerator[r]] <> "/" <> pyNumber[Denominator[r]] <> ")";
pyAtomQ[e_] := MatchQ[e, _Symbol | _PythonVariable | _PythonAtom | _Integer | _Rational | Pi | E];
pyParen[e_] := If[pyAtomQ[e] && ! (MatchQ[e, _Integer | _Rational] && e < 0), py[e], "(" <> py[e] <> ")"];
py[n_Integer] := pyNumber[n];
py[r_Rational] := pyNumber[r];
py[Pi] := "math.pi";
py[E] := "math.e";
py[PythonVariable[s_String]] := s;
py[PythonVariable[s_String, i_Integer]] := s <> "[" <> ToString[i] <> "]";
py[PythonAtom[f_String, i_Integer, r_String, j_Integer]] := "(" <> f <> "[" <> ToString[i] <> "] + " <> r <> "[" <> ToString[j] <> "])";
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

(* the coupling in terms of g = G (pyG stands for g until the placeholders are put in); Missing[] when g = G has no solution *)
pyGRule[c_, gexpr_] := Module[{sol = Quiet[Solve[pyG == gexpr, c]]}, If[MatchQ[sol, {{c -> _} ..}], Last[sol][[1]], Missing[]]];

(* the loop, conditioned for BMWflow's grids. (1) Atoms: f -> A - r for each function and the shape of its Gamma2 entry's regulator
   (pyAtomPairs), at q and at Q: large r(q), r(Q) then stay inside A = f + r instead of cancelling between expanded terms. (2) Q2 = p^2 + q^2
   + 2 p q c, BMWflow's Q^2 (where fQ and rQ are evaluated): p^(2n) -> (Q2 - q^2 - 2 p q c)^n, numerators re-expanded, until each term
   has p at most linearly; near Q -> 0 the cancellation of p^2 + 2 p q c + q^2 against large atoms is then exact. Both are identities
   (PythonExport.conditionedLoop checks it exactly). *)
pyAtomPairs[adim_] := Module[{eqs = adim["Equations"]["Equations"], reg = EmbeddedStage[adim["Equations"], "Regulator"]["Entries"], fs = adim["Functions"],
    shapes = adim["Shapes"]},
   DeleteCases[Function[r, Module[{sh = Union[Cases[Lookup[reg, Key[r["Pair"]], 0], s_[_] /; MemberQ[shapes, s] :> s, {0, Infinity}]]},
       If[Length[sh] == 1, {First[FirstPosition[fs, r["Function"]]] - 1, First[FirstPosition[shapes, First[sh]]] - 1}, Nothing]]] /@ eqs, Nothing]];
pyConditioned[e_, pairs_] := Module[{P = PythonVariable["p"], q = PythonVariable["q"], c = PythonVariable["c"], Q2 = PythonVariable["Q2"], atoms, step},
   atoms = Flatten[{PythonVariable["fq", #1] -> PythonAtom["fq", #1, "rq", #2] - PythonVariable["rq", #2],
        PythonVariable["fQ", #1] -> PythonAtom["fQ", #1, "rQ", #2] - PythonVariable["rQ", #2]} & @@@ pairs];
   step[x_] := Module[{t = Together[x //. Power[P, n_Integer /; n >= 2] :> P^Mod[n, 2] (Q2 - q^2 - 2 P q c)^Quotient[n, 2]]},
     Expand[Numerator[t]]/Denominator[t]];
   FixedPoint[step, e /. atoms, 10]];
pyUnconditioned[e_] := e /. {PythonAtom[f_, i_, r_, j_] :> PythonVariable[f, i] + PythonVariable[r, j],
    PythonVariable["Q2"] -> PythonVariable["p"]^2 + PythonVariable["q"]^2 + 2 PythonVariable["p"] PythonVariable["q"] PythonVariable["c"]};

PythonParts[adim_Association] := PythonParts[adim, None];
PythonParts[adim_Association, gexpr_, conditioned_ : True] := Block[{$FRGModule = "PythonExport"},
  Module[{fs = adim["Functions"], shapes = adim["Shapes"], aq = adim["LoopMomentum"], ap = adim["ExternalMomentum"], c = adim["Cosine"], common, cRule,
    cps = adim["Couplings"], betaL = adim["BetaL"]},
   If[gexpr =!= None,
    cRule = {pyGRule[First[cps], gexpr]};
    betaL = {Expand[(D[gexpr, First[cps]] First[betaL]) /. cRule]}; cps = {pyG},
    cRule = {}];
   common = Join[Thread[adim["Etas"] -> (PythonVariable["eta", #] & /@ Range[0, Length[adim["Etas"]] - 1])],
     Thread[cps -> (PythonVariable["couplings", #] & /@ Range[0, Length[cps] - 1])],
     If[adim["Dimension"] === 2, {}, {adim["Dimension"] -> PythonVariable["d"]}]];
   <|"linear" -> Expand[adim["Linear"] /. cRule /. pyLinearRules[fs, ap] /. common],
    "loop" -> (If[conditioned, pyConditioned[#, pyAtomPairs[adim]] &, Identity] /@ (adim["Loop"] /. cRule /. pyLoopRules[fs, shapes, aq, ap, c]) /. common),
    "eta_matrix" -> (adim["EtaMatrix"] /. cRule /. pyEtaRules[fs, shapes, aq] /. common),
    "eta_vector" -> (adim["EtaVector"] /. cRule /. pyEtaRules[fs, shapes, aq] /. common),
    "beta" -> (betaL /. common)|>]];

pyTuple[l_List] := "(" <> StringRiffle[l, ", "] <> If[Length[l] == 1, ",", ""] <> ")";
pyStrings[l_List] := pyTuple[("\"" <> ToString[#] <> "\"") & /@ l];
pyDef[name_, args_, body_] := "def " <> name <> "(" <> StringRiffle[args, ", "] <> "):\n    return " <> body <> "\n";
(* a row free of the array arguments (a constant, or eta and g only) gets "+ 0.0*arg" so that it has the arrays' shape *)
pyScalarArgs = {"eta", "couplings", "d"};
pyRowString[e_, arg_] := Which[arg === None, If[e === 0, "0.0", py[e]], e === 0, "0.0*" <> arg,
   FreeQ[e, PythonVariable[s_ /; ! MemberQ[pyScalarArgs, s], ___]], py[e] <> " + 0.0*" <> arg, True, py[e]];
pyRows[l_List, arg_ : None] := "(\n        " <> StringRiffle[pyRowString[#, arg] & /@ l, ",\n        "] <> If[Length[l] == 1, ",", ""] <> ")";

Options[ExportPython] = {"G" -> None, "Model" -> None, "File" -> None};

ExportPython[adim_, opts : OptionsPattern[]] := Block[{$FRGModule = "PythonExport"},
  Module[{model = OptionValue["Model"], file = OptionValue["File"], gexpr = OptionValue["G"], parts, checks = {}, warns = {}, all, fs, cps, shapes, etaNames,
    bad, units, header, body, defs, gSyms, fixedDim},
   FRGAssert["PythonExport.pipeline", AssociationQ[adim] && Lookup[adim, "Pipeline", None] === "LO", If[AssociationQ[adim], Lookup[adim, "Pipeline", None], adim],
    "ExportPython writes LO stages (functions of p^2) for BMWflow; a derivative-expansion stage goes to ExportJulia (FlowRG)"];
   validateStage[adim, "Adimensional"];
   FRGAssert["PythonExport.options", MatchQ[{opts}, {(Rule | RuleDelayed)["G" | "Model" | "File", _] ...}], {opts},
    "options: \"G\" -> expression of g in the coupling, \"Model\" -> label, \"File\" -> path"];
   FRGAssert["PythonExport.model", StringQ[model] && StringMatchQ[model, RegularExpression["[A-Za-z][A-Za-z0-9_]*"]], model, "\"Model\" -> \"KPZ\": a label, written as MODEL"];
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
   gSyms = Union[Cases[{gexpr}, s_Symbol /; Context @@ {s} =!= "System`", {0, Infinity}, Heads -> True]];
   AppendTo[checks, assertCheck["PythonExport.gMap", gexpr =!= None && gSyms === cps && ! ExactlyZeroQ[D[gexpr, First[cps]]] && ! MissingQ[pyGRule[First[cps], gexpr]],
     <|"G" -> gexpr, "COUPLINGS" -> cps|>,
     "\"G\" -> BMWflow's g as an invertible expression of the exported coupling alone (KPZ: \"G\" -> alamk^2)",
     "g = " <> ToString[gexpr, InputForm]]];
   parts = PythonParts[adim, gexpr];
   all = Values[parts];
   bad = Union[Cases[all, Power[b_, x_] /; ! IntegerQ[x] && ! FreeQ[b, PythonVariable["couplings", _]], {0, Infinity}]];
   AppendTo[checks, assertCheck["PythonExport.gEven", bad === {}, bad,
     "a kernel or beta is not a rational function of g = G (a radical of g is left): BMWflow's state is g, choose G accordingly",
     "every kernel and beta rational in g"]];
   AppendTo[checks, assertCheck["PythonExport.conditionedLoop", ExactlyZeroQ[pyUnconditioned[parts["loop"]] - PythonParts[adim, gexpr, False]["loop"]], "loop",
     "internal: the loop in atoms and Q2 must equal the loop in f, r, p, q, c", "loop in atoms (f + r) and Q2 = p^2 + q^2 + 2 p q c (exact identity), atoms " <>
      StringRiffle[(ToString[adim["Functions"][[#1 + 1]]] <> " + " <> ToString[adim["Shapes"][[#2 + 1]]]) & @@@ pyAtomPairs[adim], ", "]]];
   AppendTo[checks, assertCheck["PythonExport.noMeasure", FreeQ[all, Gamma | Sin | Power[Pi | 2 | 2 Pi, x_ /; ! FreeQ[x, PythonVariable["d"]]]],
     Cases[all, _Gamma | _Sin | Power[Pi | 2 | 2 Pi, x_ /; ! FreeQ[x, PythonVariable["d"]]], {0, Infinity}],
     "internal: the loop measure (v_d, v_{d-1}, Gamma(d), the (2 Pi)^d) must not be in the module: the adapter applies it",
     "no v_d, v_{d-1}, Gamma(d), (2 Pi)^d in the kernels (the adapter applies the measure)"]];
   units = Join[adim["Dimensions"]["Basis"], {adim["Equations"]["LoopMomentum"], adim["Equations"]["ExternalMomentum"], adim["LoopMomentum"], adim["ExternalMomentum"],
      adim["Cosine"]}, If[adim["Dimension"] === 2, {}, {adim["Dimension"]}]];
   bad = Union[Cases[all /. (_PythonVariable | _PythonAtom) -> 1, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> True]];
   AppendTo[checks, assertCheck["PythonExport.knownSymbols", bad === {}, bad,
     "a symbol is left after the Python argument substitution: a unit, an unexpected function atom (F' in loop, r' at Q, r''), or a coupling outside COUPLINGS",
     "only the contract arguments (p, q, c, Q2, fq, fQ, rq, rq_, rQ, f, fprime, r, r_, eta, couplings, d)"]];
   bad = Complement[Union[Cases[all /. (_PythonVariable | _PythonAtom) -> pyPlaceholder, h_Symbol[___] /; Context[h] === "System`" :> h, {0, Infinity}, Heads -> True]],
     {Plus, Times, Power, Rational, List}];
   AppendTo[checks, assertCheck["PythonExport.systemHeads", bad === {}, bad, "a System function without a Python translation is left", "only + * ** (np.sqrt)"]];

   fixedDim = If[IntegerQ[adim["Dimension"]], ToString[adim["Dimension"]], "None"];
   header = StringRiffle[{
      "# Generated by FRGscalarMSR ExportPython. Do not hand-edit: regenerate.",
      "# Source: LO Adimensional stage, model " <> model <> "; content hash " <> pyHashPlaceholder <> " (of everything below the header).",
      "# Contract: FRGscalarMSR docs/handoff/BMWflow_LO.md, version " <> ToString[pyContractVersion] <> ". Pure numpy + math, elementwise, no state.",
      "# Dimensionless; p, q in units of the RG scale; f_X(0) = 1; fprime = d f/d|momentum|; r_ = dr/d(q^2); every rate is d_s x, s = ln(kappa/Lambda).",
      "# The loop measure is NOT included: I_X(p) = v_{d-1}/(2 pi) Int dq q^(d-1) Int dtheta sin^(d-2)(theta) loop_X; eta: <.> = v_d Int dq q^(d-1) (.).",
      "# loop_X is already divided by p^(2 P2_DIVIDED[X]) (non-finite at p = 0 for those rows: the adapter extrapolates).",
      "# loop takes Q2 = p^2 + q^2 + 2 p q c, the Q^2 at which fQ and rQ are evaluated; it is written in Q2 and in the atoms (f + r): well conditioned as Q -> 0.",
      "# The coupling is BMWflow's g = " <> ToString[gexpr, InputForm] <> "; beta is d_s g."}, "\n"];
   defs = {
     pyDef["linear", {"p", "f", "fprime", "eta", "couplings", "d"}, pyRows[parts["linear"], "p"]],
     pyDef["loop", {"p", "q", "c", "Q2", "fq", "fQ", "rq", "rq_", "rQ", "eta", "couplings", "d"}, pyRows[parts["loop"], "p"]],
     pyDef["eta_matrix", {"q", "f", "fprime", "r", "r_", "couplings", "d"}, "(\n        " <> StringRiffle[pyRows[#, "q"] & /@ parts["eta_matrix"], ",\n        "] <>
       If[Length[parts["eta_matrix"]] == 1, ",", ""] <> ")"],
     pyDef["eta_vector", {"q", "f", "fprime", "r", "r_", "couplings", "d"}, pyRows[parts["eta_vector"], "q"]],
     pyDef["beta", {"couplings", "eta", "d"}, pyRows[parts["beta"]]]};
   Quiet[CreateDirectory[DirectoryName[file]], CreateDirectory::eexist];
   body = StringRiffle[Join[{"import math", "", "import numpy as np", "",
       "CONTRACT_VERSION = " <> ToString[pyContractVersion],
       "FUNCTIONS = " <> pyStrings[fs], "ETAS = " <> pyStrings["eta" <> # & /@ etaNames], "COUPLINGS = (\"g\",)", "SHAPES = " <> pyStrings[shapes],
       "P2_DIVIDED = " <> pyTuple[ToString /@ adim["P2Divided"]], "FIXED_DIMENSION = " <> fixedDim, "APPROXIMATION = \"LO\"", "MODEL = \"" <> model <> "\"", "", ""},
      Riffle[defs, "\n"]], "\n"];
   Export[file, StringReplace[header, pyHashPlaceholder -> StringTake[IntegerString[Hash[body, "SHA256"], 16, 64], 12]] <> "\n\n" <> body, "Text"];
   AppendTo[warns, infoCheck["INFO", "PythonExport.adapter", "the adapter applies the measure, the p = 0 extrapolation of P2_DIVIDED rows and binds each regulator shape by name"]];
   Print["ExportPython: wrote ", file, " (", Length[fs], " functions, ", Length[cps], " coupling, g = ", ToString[gexpr, InputForm], ", model ", model, ")"];
   Scan[Print[StringPadRight[#Status, 11], #Check, If[#Detail === "", "", "   " <> ToString[#Detail]]] &, Join[checks, warns]];
   file]];

End[];
EndPackage[];
