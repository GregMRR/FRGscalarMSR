(* ::Package:: *)

(* JuliaExport.wl -- ExportJulia: write the Adimensional stage, at one integer space dimension, as a FlowRG equations file
   (~/Physics/Sources/FlowRG.jl/docs/equations-interface.md): COUPLINGS, ETAS, SHAPES, MUST_FIX, betal, betanl, etaMl, etaMnl,
   etabl, etabnl, with positional arguments; plus causal_poles (design/CAUSALITY.md), an extra definition FlowRG does not call.
   Objective, assumptions and checks: docs/MODULES.md section "JuliaExport". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

ExportJulia::usage = "ExportJulia[adim, \"Dimension\" -> d, \"File\" -> path] writes the FlowRG equations file for space dimension d and \
returns its path. \"Form\" -> \"Compact\" (default: dimensionless atoms as local variables) or \"Explicit\". Shape derivatives are written \
Darmh(aq^2), D2armh(aq^2). causal_poles(couplings..., aq) gives the causal-pole conditions: each entry must be > 0 at every aq.";
JuliaParts::usage = "JuliaParts[adim, d] gives the parts of the equations file at dimension d (Mathematica expressions, shapes and atoms \
still symbolic): <|\"betal\" -> ..., \"betanl\" -> ..., ...|>, as ExportJulia writes them.";

Begin["`Private`"];

(* Julia name of a dimensionless atom: AdimAtom[H, Hbar, 0, 1] -> atom_H_Hbar_0_1 *)
atomName[AdimAtom[a__]] := "atom_" <> StringRiffle[ToString /@ {a}, "_"];

(* shape derivatives as applied functions: arnh'[y] -> Darnh[y], arnh''[y] -> D2arnh[y] *)
derivativeRule[shapes_] := Derivative[n_Integer][f : Alternatives @@ shapes][y_] :>
   (FRGAssert["JuliaExport.shapeDerivative", n <= 2, Derivative[n][f],
     "only r, r' and r'' are exported: p^4 (and higher) flows need r''' / r'''' and a FlowRG interface extension (docs/FUTURE.md)"];
    Symbol[Context[f] <> If[n == 1, "D", "D2"] <> SymbolName[f]][y]);

toJulia[e_] := Module[{s},
   s = ToString[e /. {Pi -> juliaPi, E -> juliaE}, InputForm, PageWidth -> Infinity];
   s = StringReplace[s, {"[" -> "(", "]" -> ")", "juliaPi" -> "pi", "juliaE" -> "MathConstants.e",
      "Sqrt(" -> "sqrt(", "Exp(" -> "exp(", "Log(" -> "log(", "Abs(" -> "abs("}];
   StringReplace[s, RegularExpression["FRGscalarMSR`([A-Za-z]+`)*"] -> ""]];

symTuple[l_List] := "(" <> StringRiffle[":" <> # & /@ l, ", "] <> If[Length[l] == 1, ",", ""] <> ")";
tupleString[exprs_List] := "(" <> StringRiffle[exprs, ",\n     "] <> If[Length[exprs] == 1, ",", ""] <> ")";

(* the parts at dimension d: Dim -> d, shape derivatives renamed; loop parts include the loop measure *)
JuliaParts[adim_Association, d_Integer, form_String : "Compact"] := Block[{$FRGModule = "JuliaExport"},
  Module[{dimSym = adim["Dimension"], aq = adim["LoopMomentum"], w, sub, nl, mnl, bnl},
   sub[e_] := (e /. If[dimSym === 2, {}, {dimSym -> d}]) /. derivativeRule[adim["Shapes"]];
   w = adim["Measure"] /. If[dimSym === 2, {}, {dimSym -> d}];
   {nl, mnl, bnl} = If[form === "Compact", adim /@ {"BetaNL", "EtaMNL", "EtaBNL"}, adim /@ {"BetaNLExplicit", "EtaMNLExplicit", "EtaBNLExplicit"}];
   <|"betal" -> sub[adim["BetaL"]], "betanl" -> sub[w nl], "etaMl" -> Flatten[adim["EtaML"]], "etaMnl" -> sub[w Flatten[mnl]],
    "etabl" -> adim["EtaBL"], "etabnl" -> sub[w bnl],
    "causalPoles" -> sub[adim["CausalPoles"][If[form === "Compact", "Compact", "Explicit"]]],
    "atoms" -> If[form === "Compact", sub[Normal[adim["Atoms"]]], {}], "measure" -> w|>]];

Options[ExportJulia] = {"Dimension" -> None, "File" -> None, "Form" -> "Compact"};

ExportJulia[adim_, opts : OptionsPattern[]] := Block[{$FRGModule = "JuliaExport"},
  Module[{t0 = AbsoluteTime[], d = OptionValue["Dimension"], file = OptionValue["File"], form = OptionValue["Form"], dimSym, parts, cps, etaSyms,
    etaNames, shapes, aq, units, all, bad, usedShapes, d2Q, header, defs, def, dims, checks = {}, warns = {}, atomKeys, placeholder, rename, jl},
   FRGAssert["JuliaExport.pipeline", NoneTrue[{adim}, AssociationQ[#] && Lookup[#, "Pipeline", None] === "LO" &], "LO", "FlowRG files are for the derivative expansion; export an LO stage with ExportPython (BMWflow)"];
   validateStage[adim, "Adimensional"];
   FRGAssert["JuliaExport.options", MatchQ[{opts}, {(Rule | RuleDelayed)["Dimension" | "File" | "Form", _] ...}], {opts},
    "options: \"Dimension\" -> d, \"File\" -> path, \"Form\" -> \"Compact\" | \"Explicit\""];
   dimSym = adim["Dimension"];
   FRGAssert["JuliaExport.dimension", IntegerQ[d] && d >= 1, d, "\"Dimension\" -> d, a positive integer (one file per d, FlowRG rule 5)"];
   FRGAssert["JuliaExport.declaredDimension", dimSym =!= 2 || d == 2, d, "the theory declares d = 2 (FixedDimension[2]): export with \"Dimension\" -> 2"];
   FRGAssert["JuliaExport.file", StringQ[file] && StringEndsQ[file, ".jl"], file, "\"File\" -> \".../EquationsND.jl\""];
   FRGAssert["JuliaExport.form", MemberQ[{"Compact", "Explicit"}, form], form, "\"Form\" -> \"Compact\" or \"Explicit\""];
   dims = adim["Dimensions"];
   parts = JuliaParts[adim, d, form];
   cps = adim["Couplings"]; etaSyms = adim["Etas"]; etaNames = adim["EtaNames"]; aq = adim["LoopMomentum"];
   shapes = adim["Shapes"];

   (* check 7: no unit, no dimensionful loop momentum, no symbolic dimension; only declared symbols *)
   units = Join[dims["Basis"], {adim["Equations"]["LoopMomentum"]}, If[dimSym === 2, {}, {dimSym}]];
   all = Values[KeyDrop[parts, "measure"]];
   AppendTo[checks, assertCheck["JuliaExport.noUnits", FreeQ[all, Alternatives @@ units], Select[units, ! FreeQ[all, #] &],
     "internal: a unit, q or the symbolic dimension is left", "free of k, references, q" <> If[dimSym === 2, "", " and " <> ToString[dimSym]]]];
   usedShapes = Select[shapes, ! FreeQ[all, #] || ! FreeQ[all, Symbol[Context[#] <> "D" <> SymbolName[#]]] || ! FreeQ[all, Symbol[Context[#] <> "D2" <> SymbolName[#]]] &];
   bad = Complement[Union[Cases[all /. _AdimAtom -> 1, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> True]],
     cps, etaSyms, {aq}, shapes, Symbol[Context[#] <> "D" <> SymbolName[#]] & /@ shapes, Symbol[Context[#] <> "D2" <> SymbolName[#]] & /@ shapes];
   AppendTo[checks, assertCheck["JuliaExport.knownSymbols", bad === {}, bad, "internal: undeclared symbols in the exported expressions",
     "only couplings, etas, aq, shapes and atoms"]];
   bad = Complement[Union[Cases[all /. _AdimAtom -> 1, h_Symbol[___] /; Context[h] === "System`" :> h, {0, Infinity}, Heads -> True]],
     {Plus, Times, Power, Rational, Complex, List, Rule, Sqrt, Exp, Log, Abs}];
   AppendTo[checks, assertCheck["JuliaExport.systemHeads", bad === {}, bad, "a System function without a Julia translation is left", "only + * ^ and sqrt/exp/log/abs"]];
   AppendTo[checks, assertCheck["JuliaExport.linearNoMomentum", FreeQ[parts /@ {"betal", "etaMl", "etabl"}, aq], parts["betal"],
     "internal: an l part depends on aq", "l parts free of aq"]];
   d2Q = ! FreeQ[all, s_Symbol /; StringStartsQ[SymbolName[s], "D2"] && MemberQ[SymbolName /@ shapes, StringDrop[SymbolName[s], 2]]];
   If[d2Q, AppendTo[warns, infoCheck["WARN", "JuliaExport.secondDerivatives", "the file uses D2<shape> (r''): FlowRG must bind second derivatives (docs/handoff/FlowRG-MUST_FIX.md, task 2)"]]];
   If[adim["MustFix"] =!= {}, AppendTo[warns, infoCheck["INFO", "JuliaExport.mustFix", "MUST_FIX = " <> symTuple[adim["MustFix"]] <>
       ": load with fixed_etas = (" <> StringRiffle[# <> " = <value>" & /@ adim["MustFix"], ", "] <> ",)"]]];

   (* one Julia definition; compact loop parts compute the atoms they use first *)
   (* atoms are written through placeholder symbols atomJ<i>Z, renamed to atom_a_b_m_j in the final string *)
   atomKeys = parts["atoms"][[All, 1]];
   placeholder = AssociationThread[atomKeys -> (Symbol["FRGscalarMSR`Private`atomJ" <> ToString[#] <> "Z"] & /@ Range[Length[atomKeys]])];
   rename = Thread[(SymbolName /@ Values[placeholder]) -> (atomName /@ atomKeys)];
   jl[e_] := StringReplace[toJulia[e /. a_AdimAtom :> placeholder[a]], rename];
   def[name_, argSyms_, exprs_, atomsQ_] := Module[{as = If[atomsQ, Select[parts["atoms"], ! FreeQ[exprs, #[[1]]] &], {}], body},
     body = tupleString[If[# === 0, "0.0", jl[#]] & /@ exprs];
     If[as === {},
      name <> "(" <> StringRiffle[ToString /@ argSyms, ", "] <> ") =\n    " <> body,
      "function " <> name <> "(" <> StringRiffle[ToString /@ argSyms, ", "] <> ")\n" <>
       StringJoin[("    " <> atomName[#[[1]]] <> " = " <> jl[#[[2]]] <> "\n") & /@ as] <> "    return " <> body <> "\nend"]];
   defs = {
     def["betal", Join[cps, etaSyms], parts["betal"], False],
     If[AllTrue[parts["betanl"], # === 0 &], Nothing, def["betanl", Join[cps, etaSyms, {aq}], parts["betanl"], True]],
     def["etaMl", cps, parts["etaMl"], False],
     If[AllTrue[parts["etaMnl"], # === 0 &], Nothing, def["etaMnl", Append[cps, aq], parts["etaMnl"], True]],
     If[AllTrue[parts["etabl"], # === 0 &], Nothing, def["etabl", cps, parts["etabl"], False]],
     If[AllTrue[parts["etabnl"], # === 0 &], Nothing, def["etabnl", Append[cps, aq], parts["etabnl"], True]],
     If[parts["causalPoles"] === {}, Nothing, def["causal_poles", Append[cps, aq], parts["causalPoles"], True]]};
   If[parts["causalPoles"] === {}, AppendTo[warns, infoCheck["INFO", "JuliaExport.causalPoles", "causal poles trivially satisfied: no causal_poles written"]]];

   header = StringRiffle[{
      "# Generated by FRGscalarMSR ExportJulia on " <> DateString[{"Year", "-", "Month", "-", "Day", " ", "Time"}] <> ". Do not hand-edit: regenerate.",
      "# Source: Adimensional stage " <> StringTake[adim["Fingerprint"], 8] <> ", space dimension d = " <> ToString[d] <> ", form " <> form <> ".",
      "# FlowRG equations-file interface (FlowRG.jl/docs/equations-interface.md).",
      "# References: " <> StringRiffle[(ToString[#] <> " (eta" <> dims["EtaNames"][#] <> ", " <> adim["ReferenceKinds"][#] <> ")") & /@ adim["References"], ", "] <>
       ";  d_s ref = -eta ref.",
      "# Loop parts (betanl, etaMnl, etabnl) include the measure " <> toJulia[parts["measure"]] <> " = S_d aq^(d-1)/(2 pi)^d.",
      If[parts["causalPoles"] === {}, Nothing, "# causal_poles(couplings..., aq): each entry must be > 0 at every aq (causal poles of the propagator; " <>
        "references > 0 assumed). Not called by FlowRG."],
      If[adim["MustFix"] === {}, Nothing, "# MUST_FIX: these etas have a placeholder row eta = 0; give their value with load_model(...; fixed_etas = (name = value,))."],
      If[d2Q, "# Uses second shape derivatives D2<shape>(y) = r''(y): needs a FlowRG version that binds them.", Nothing]}, "\n"];
   Quiet[CreateDirectory[DirectoryName[file]], CreateDirectory::eexist];
   Export[file, StringRiffle[Join[{header, "",
       "const COUPLINGS = " <> symTuple[ToString /@ cps],
       "const ETAS = " <> symTuple[etaNames],
       "const SHAPES = " <> symTuple[ToString /@ usedShapes]},
      If[adim["MustFix"] === {}, {}, {"const MUST_FIX = " <> symTuple[adim["MustFix"]]}], {""}, Riffle[defs, ""]], "\n"] <> "\n", "Text"];
   Print["ExportJulia: wrote ", file, " (", Length[cps], " couplings, ", Length[etaSyms], " etas, d = ", d, ", ", form, ")"];
   Scan[Print[StringPadRight[#Status, 11], #Check, If[#Detail === "", "", "   " <> ToString[#Detail]]] &, Join[checks, warns]];
   file]];

End[];
EndPackage[];
