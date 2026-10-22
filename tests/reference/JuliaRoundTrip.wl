(* ::Package:: *)

(* Julia round trip of an exported equations file (plan D14). juliaRoundTrip[name, adim, d] exports adim at dimension d,
   evaluates every part of the file in plain Julia (smooth test shapes r_i(y) = 1/(i + y)^2 and their derivatives, defined next to
   the file, at a random rational state, etas and aq) and compares with the same parts evaluated in Mathematica (relative 1e-11).
   Then, if FlowRG is available and the file needs no second shape derivative, load_model must accept it and give finite
   etas and beta (Litim; fixed etas from MUST_FIX set to 0.1). Without julia or FlowRG the checks are reported SKIP. *)

(* julia runs with minimal compilation: the generated functions are huge and are evaluated only a few times here *)
$juliaFlags = {"--startup-file=no", "--compile=min", "-O0"};
$juliaExe = Quiet[Check[If[RunProcess[{"julia", "--version"}, "ExitCode"] === 0, "julia", None], None]];
$flowRGPath = Replace[Environment["FLOWRG_PATH"], $Failed -> FileNameJoin[{$HomeDirectory, "Physics", "Sources", "FlowRG.jl"}]];
$roundTripDir = FileNameJoin[{$TemporaryDirectory, "FRGscalarMSR-roundtrip-" <> ToString[$ProcessID]}];
Quiet[CreateDirectory[$roundTripDir], CreateDirectory::eexist];

testShape[i_] := Function[y, 1/(i + y)^2];
(* a test shape (or a derivative) as a Julia definition *)
juliaShapeDef[name_String, f_] := "    " <> name <> "(y) = " <> StringReplace[ToString[f[Symbol["y"]], InputForm], {"[" -> "(", "]" -> ")"}];
fromJulia[s_String] := ToExpression[StringReplace[s, {"e" -> "*^", "NaN" -> "Indeterminate", "Inf" -> "Infinity"}]];

juliaRoundTrip[name_String, adim_, d_Integer] := Module[{file, parts, shapes, shapeRules, pt, cps, etas, aq, mma, script, out, jl, rel, d2Q, fixed, fl},
   If[$juliaExe === None, skip[name <> ": Julia round trip", "julia not found"]; Return[Null]];
   file = quiet[ExportJulia[adim, "Dimension" -> d, "File" -> FileNameJoin[{$roundTripDir, name, "Equations" <> ToString[d] <> "D.jl"}]]];
   parts = JuliaParts[adim, d];
   shapes = adim["Shapes"]; cps = adim["Couplings"]; etas = adim["Etas"]; aq = adim["LoopMomentum"];
   (* shape i and its derivatives, in Mathematica and as Julia definitions *)
   shapeRules = Flatten[Table[With[{f = testShape[i], s = shapes[[i]]},
       {s -> f, Symbol[Context[s] <> "D" <> SymbolName[s]] -> Function[y, Evaluate[D[f[y], y]]],
        Symbol[Context[s] <> "D2" <> SymbolName[s]] -> Function[y, Evaluate[D[f[y], {y, 2}]]]}], {i, Length[shapes]}]];
   SeedRandom[17];
   pt = Thread[Join[cps, etas, {aq}] -> RandomInteger[{1, 40}, Length[cps] + Length[etas] + 1]/RandomInteger[{41, 97}, Length[cps] + Length[etas] + 1]];
   mma = N[Flatten[(parts /@ {"betal", "betanl", "etaMl", "etaMnl", "etabl", "etabnl", "causalPoles"}) /. parts["atoms"] /. shapeRules /. (#[[1]] -> Rationalize[N[#[[2]]], 0] & /@ pt)], 30];
   script = StringRiffle[Join[{"module M"},
      juliaShapeDef[ToString[#[[1]]], #[[2]]] & /@ shapeRules,
      {"    include(" <> ToString[file, InputForm] <> ")", "end",
       "c = (" <> StringRiffle[ToString[N[#], CForm] & /@ (cps /. pt), ", "] <> ",)",
       "e = (" <> StringRiffle[ToString[N[#], CForm] & /@ (etas /. pt), ", "] <> ",)",
       "aq = " <> ToString[N[aq /. pt], CForm],
       "z(f, k, args...) = isdefined(M, f) ? collect(getfield(M, f)(args...)) : zeros(k)",
       "N, n = " <> ToString[Length[cps]] <> ", " <> ToString[Length[etas]],
       "v = vcat(collect(M.betal(c..., e...)), z(:betanl, N, c..., e..., aq), collect(M.etaMl(c...)), z(:etaMnl, n*n, c..., aq), z(:etabl, n, c...), z(:etabnl, n, c..., aq), z(:causal_poles, " <> ToString[Length[parts["causalPoles"]]] <> ", c..., aq))",
       "println(join(repr.(v), \",\"))"}], "\n"];
   Export[FileNameJoin[{$roundTripDir, name, "roundtrip.jl"}], script, "Text"];
   out = RunProcess[Join[{$juliaExe}, $juliaFlags, {FileNameJoin[{$roundTripDir, name, "roundtrip.jl"}]}]];
   test[name <> ": julia evaluates the exported file", out["ExitCode"], 0];
   If[out["ExitCode"] =!= 0, Print["    julia stderr: ", StringTake[out["StandardError"], UpTo[600]]]; Return[Null]];
   jl = fromJulia["{" <> StringTrim[out["StandardOutput"]] <> "}"];
   rel = Max[MapThread[Abs[#1 - #2]/Max[1, Abs[#2]] &, {jl, mma}]];
   test[name <> ": Julia = Mathematica for every part (rel. 1e-11)", Length[jl] == Length[mma] && rel < 10^-11, True];
   (* FlowRG *)
   d2Q = StringContainsQ[Import[file, "Text"], RegularExpression["\\bD2\\w+\\("]];
   Which[
    ! DirectoryQ[$flowRGPath], skip[name <> ": FlowRG load_model", "FlowRG not found at " <> $flowRGPath],
    d2Q, skip[name <> ": FlowRG load_model", "the file uses D2<shape> (r''): needs FlowRG task 2 of docs/handoff/FlowRG-MUST_FIX.md"],
    True,
    fixed = If[adim["MustFix"] === {}, "(;)", "(" <> StringRiffle[# <> " = 0.1" & /@ adim["MustFix"], ", "] <> ",)"];
    fl = RunProcess[Join[{$juliaExe}, $juliaFlags, {"--project=" <> $flowRGPath, "-e",
       "using FlowRG; m = load_model(" <> ToString[file, InputForm] <> "; regulator = Litim, fixed_etas = " <> fixed <> "); " <>
        "u = fill(0.3, length(m.names)); h = values(etas(m, u)); b = beta(m, u); println(all(isfinite, h) && all(isfinite, b))"}]];
    test[name <> ": FlowRG load_model accepts the file, finite etas and beta", {fl["ExitCode"], StringTrim[fl["StandardOutput"]]}, {0, "true"}]]];
