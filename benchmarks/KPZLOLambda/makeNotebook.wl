(* ::Package:: *)

(* makeNotebook.wl -- writes KPZLOLambda.nb (this folder) from the cell list below.
   Run:  math -script ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/makeNotebook.wl
   Check that it runs headless:
         math -script ~/sources/FRGscalarMSR/examples/runNotebook.wl ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/KPZLOLambda.nb
   Edit the cells here, not in the .nb. The pipeline cells mirror Pipeline.wl; the last section checks that they build the same stage. *)

dir = DirectoryName[$InputFileName];
title[s_] := Cell[s, "Title"];
section[s_] := Cell[s, "Section"];
text[s_] := Cell[s, "Text"];
input[s_] := Cell[s, "Input"];

cells = {
   title["Benchmark: KPZ at LO with a flowing f_lambda"],
   text["The Kloss-Canet-Wschebor LO ansatz with three functions of p^2: f_lambda multiplies the Galilean-covariant combination \
dt h - lam/2 (grad h)^2, so it sits at i w in Gamma(1,1) and in the cubic vertex, at the momentum of the hbar leg. Its flow is the \
i nu coefficient of d_s Gamma(1,1)(nu, p): d/d(i nu) before the frequency integral (design/LOFLAMBDA.md)."],
   input["Get[\"~/sources/FRGscalarMSR/FRGscalarMSR.wl\"];\nkpzlamDir = ExpandFileName[\"~/sources/FRGscalarMSR/benchmarks/KPZLOLambda\"];"],

   section["1. Ansatz"],
   text["\"Functions\" -> {..., flam -> 1}: f_lambda has reference 1. Its f(0) = 1 is the Galilean Ward identity (the Protected i w tag), \
so it has no eta; BuildLOEquations checks that its flow vanishes exactly at p = 0."],
   input["kpzlamTheory = DefineTheory[\"Scalars\" -> {h}, \"Dimension\" -> Dim];\nkpzlamAnsatz = DefineAnsatz[kpzlamTheory,\n  \"Density\" -> fn[flam, hbar] (dt[h] - lamk/2 grad[h].grad[h]) - nuk hbar fn[fnu, lap[h]] - Dk hbar fn[fD, hbar],\n  \"Couplings\" -> {nuk, lamk, Dk}, \"Functions\" -> {fD -> Dk, fnu -> nuk, flam -> 1}];"],
   input["kpzlamSymmetries = DeclareSymmetries[kpzlamAnsatz,\n  Protected[{hbar, h}, \"Frequency\" -> 1, \"Momentum\" -> 0, \"Reason\" -> \"Galilean Ward identity: f_lambda(0) = 1\"],\n  Protected[{h, hbar}, \"Frequency\" -> 0, \"Momentum\" -> 0, \"Reason\" -> \"h -> h + c (shift symmetry)\"],\n  NonRenormalised[lamk, \"Reason\" -> \"Galilean invariance\"]];"],

   section["2. Vertices, propagator, loops, LO equations"],
   text["The propagator now has i w f_lambda(q^2): Delta(w) = c0 + c1 (i w) with c1 = f_lambda (the causal-pole assumption is f_lambda > 0, on \
the ExternalMomentum card). The f_lambda row has double poles at Q (d/d nu of 1/Delta(w + nu, Q))."],
   input["kpzlamVertices = BuildVertices[kpzlamAnsatz];\nkpzlamRegulator = DefineRegulator[kpzlamAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], \"Shapes\" -> {rn, rd}];\nkpzlamPropagator = BuildPropagator[kpzlamVertices, kpzlamRegulator];"],
   input["kpzlamLoops = BuildLoopIntegrands[kpzlamVertices, kpzlamPropagator, kpzlamRegulator];\nkpzlamEquations = BuildLOEquations[kpzlamLoops, kpzlamSymmetries];"],
   input["RunOracle[kpzlamEquations];"],
   input["kpzlamDims = DeriveDimensions[kpzlamAnsatz, \"References\" -> {nuk -> \"nu\", Dk -> \"D\"}];\nkpzlamAdim = AdimensionalEquations[kpzlamEquations, kpzlamDims];\n{kpzlamAdim[\"Functions\"], kpzlamAdim[\"HasEta\"], kpzlamAdim[\"Etas\"], kpzlamAdim[\"Linear\"]}"],

   section["3. References and checks"],
   text["f_lambda = 1 against benchmarks/KPZLO; generic f_lambda against the Kloss et al. kernels (Reference.wl); the derivative expansion \
at f_lambda = 1 + x p^2 against the Taylor coefficients of the three rows (sign and normalisation of the d/d(i nu) row)."],
   input["Block[{Print}, Get[FileNameJoin[{ParentDirectory[kpzlamDir], \"KPZLO\", \"Pipeline.wl\"}]]; Get[FileNameJoin[{kpzlamDir, \"Pipeline.wl\"}]]];\nGet[FileNameJoin[{kpzlamDir, #}]] & /@ {\"Reference.wl\", \"Checks.wl\"};\nkpzlamResults = kpzlamChecks[];\nkpzlamSummary[kpzlamResults]"],
   input["If[AllTrue[kpzlamResults, #Pass &], \"all checks pass\", Abort[]]"],

   section["4. Export (contract v3)"],
   input["ExportPython[kpzlamAdim, \"G\" -> alamk^2, \"Model\" -> \"KPZLambda\", \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZLambda_LO.py\"}]]"]};

Export[FileNameJoin[{dir, "KPZLOLambda.nb"}], Notebook[cells, StyleDefinitions -> "Default.nb"]];
