(* ::Package:: *)

(* Adimensional + JuliaExport acceptance on shallow-water DE4plus, against the hand-built export of the old notebook
   ManualCodeOld/DE4plus_largscale_forcing_LLM/AdimensionalFlowEQGenericRegLargeScale.nb -> eqs/Equations2D.jl (copied to
   tests/reference/legacy/DE4plusEquations2D.jl). Same state order and names; ETAS by name (old (nu, phi, h), here (phi, h, nu)).
   Known differences, fixed in advance:
   - the old file has no (2 Pi)^-d in the measure: every loop part here = old loop part x (2 Pi)^-2 (plan D17, review E7);
   - the old file multiplies d_s M_h by eta_phi instead of eta_h (plan D9): compared at eta_h = eta_phi, and the eta system as
     (M eta, b) at such etas. *)

Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];
Get[FileNameJoin[{testRoot, "reference", "JuliaRoundTrip.wl"}]];

de4Dims = quiet[DeriveDimensions[swaAnsatz, "References" -> {Dphik -> {"phi", "Mode" -> "Running"}, Dhk -> {"h", "Mode" -> "Running"},
     nk -> {"nu", "Mode" -> "Running", "Equation" -> nk ank}}, "Dimensionless" -> {Ak, ank, aDhk, aDphik}]];
de4Adim = quiet[AdimensionalEquations[swaEquations, de4Dims, "State" -> {gk, g2k, lk, Hk, aDphik, D2phik, Kk, Ak, ank, aDhk}]];

test["DE4 adim: state names = old COUPLINGS", de4Adim["Couplings"], {agk, ag2k, alk, aHk, aDphik, aD2phik, aKk, Ak, ank, aDhk}];
test["DE4 adim: betal = old betal at d = 2", Expand[de4Adim["BetaL"] /. Dim -> 2],
  Expand[{agk (-1 - etah/2 + etanu + etaphi/2), ag2k (1 - etah/2 + etanu + etaphi/2), alk ((3 etanu)/2 - etaphi/2), aHk (-1 + etah/2 - etanu/2),
    aDphik etaphi, aD2phik (2 + etaphi), aKk (1 + etah/2 + etaphi/2), Ak etanu, ank etanu, aDhk etah}]];
test["DE4 adim: all references running, nothing to fix", {Values[de4Adim["ReferenceKinds"]], de4Adim["MustFix"]}, {{"Running", "Running", "Running"}, {}}];
test["DE4 adim: causal poles (design/CAUSALITY.md section 4)", Expand[de4Adim["CausalPoles"]["Explicit"]],
  Expand[{Ak + ank + arnh[aq^2] + arnphi[aq^2], agk aHk alk + aq^2 (ag2k aHk alk + (Ak + arnh[aq^2]) (ank + arnphi[aq^2]))}]];

(* old vs new in Julia, same smooth test shapes, random state, eta_h = eta_phi *)
de4Compare[] := Module[{newFile, oldFile, dir, shapes, defs, script, out, v},
   If[$juliaExe === None, skip["DE4 adim: old vs new export", "julia not found"]; Return[Null]];
   dir = FileNameJoin[{$roundTripDir, "DE4old"}]; Quiet[CreateDirectory[dir], CreateDirectory::eexist];
   newFile = quiet[ExportJulia[de4Adim, "Dimension" -> 2, "File" -> FileNameJoin[{dir, "Equations2D.jl"}]]];
   oldFile = FileNameJoin[{testRoot, "reference", "legacy", "DE4plusEquations2D.jl"}];
   shapes = de4Adim["Shapes"];
   defs = StringRiffle[Flatten[Table[With[{f = testShape[i], s = ToString[shapes[[i]]]},
       {juliaShapeDef[s, f], juliaShapeDef["D" <> s, Function[y, Evaluate[D[f[y], y]]]], juliaShapeDef["D2" <> s, Function[y, Evaluate[D[f[y], {y, 2}]]]]}],
      {i, Length[shapes]}]], "\n"] <> "\n";
   script = "module New\n" <> defs <> "    include(" <> ToString[newFile, InputForm] <> ")\nend\nmodule Old\n" <> defs <>
     "    include(" <> ToString[oldFile, InputForm] <> ")\nend\n" <> "
c = (0.31, 0.17, 0.43, 0.29, 0.61, 0.23, 0.37, 0.53, 0.47, 0.71)
nu, phi = 0.13, 0.27; h = phi                          # eta_h = eta_phi neutralises the old d_s M_h bug
aq = 0.41
r = (2pi)^2                                            # the old file has no (2 Pi)^-d
newEta, oldEta = (phi, h, nu), (nu, phi, h)
p = [3, 1, 2]                                          # old row/column i = new row/column p[i]
bl = maximum(abs.(collect(New.betal(c..., newEta...)) .- collect(Old.betal(c..., oldEta...))))
bnl = maximum(abs.(r .* collect(New.betanl(c..., newEta..., aq)) .- collect(Old.betanl(c..., oldEta..., aq))) ./ max.(1, abs.(collect(Old.betanl(c..., oldEta..., aq)))))
Mn = reshape(collect(New.etaMnl(c..., aq)), 3, 3)'      # row by row
Mo = reshape(collect(Old.etaMnl(c..., aq)), 3, 3)'
bn = collect(New.etabnl(c..., aq)); bo = collect(Old.etabnl(c..., aq))
etaVec = [nu, phi, h]
eta = r .* (Mn * collect(newEta))[p] .- Mo * etaVec
bb = r .* bn[p] .- bo
println(join(string.([bl, bnl, maximum(abs.(eta) ./ max.(1, abs.(Mo * etaVec))), maximum(abs.(bb) ./ max.(1, abs.(bo)))]), \",\"))
";
   Export[FileNameJoin[{dir, "compare.jl"}], script, "Text"];
   out = RunProcess[Join[{$juliaExe}, $juliaFlags, {FileNameJoin[{dir, "compare.jl"}]}]];
   test["DE4 adim: julia runs old and new files", out["ExitCode"], 0];
   If[out["ExitCode"] =!= 0, Print["    julia stderr: ", StringTake[out["StandardError"], UpTo[800]]]; Return[Null]];
   v = fromJulia["{" <> StringTrim[out["StandardOutput"]] <> "}"];
   test["DE4 adim: betal = old (exact up to rounding)", v[[1]] < 10^-13, True];
   test["DE4 adim: betanl = old x (2 Pi)^-2 at eta_h = eta_phi", v[[2]] < 10^-10, True];
   test["DE4 adim: eta system M eta = old x (2 Pi)^-2 at eta_h = eta_phi", v[[3]] < 10^-10, True];
   test["DE4 adim: eta system b = old x (2 Pi)^-2", v[[4]] < 10^-10, True]];
de4Compare[];
juliaRoundTrip["DE4plus2D", de4Adim, 2];
