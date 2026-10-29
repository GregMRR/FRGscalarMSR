(* ::Package:: *)

(* makeNotebook.wl -- writes KPZLO.nb (this folder) from the cell list below.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZLO/makeNotebook.wl
   Check that it runs headless:
         math -script ~/Physics/Sources/FRGscalarMSR/examples/runNotebook.wl ~/Physics/Sources/FRGscalarMSR/benchmarks/KPZLO/KPZLO.nb
   Edit the cells here, not in the .nb. The pipeline cells mirror Pipeline.wl; section 6 checks that they build the same stage. *)

dir = DirectoryName[$InputFileName];
title[s_] := Cell[s, "Title"];
section[s_] := Cell[s, "Section"];
text[s_] := Cell[s, "Text"];
input[s_] := Cell[s, "Input"];

cells = {
   title["Benchmark: KPZ at LO (functions of p^2)"],
   text["Reference: the KPZ LO flow of BMWflow, written in notes/LO_flow_equations.tex and notes/kpz_notation.tex (Kloss, Canet, \
Wschebor, PRE 86, 051124 (2012), with f_lambda = 1 and the functions at zero frequency) and coded by hand in flow/models/model_kpz.py. \
This notebook derives the same flow with FRGscalarMSR from the ansatz alone and compares it with the notes, exactly, at every \
external momentum p, loop momentum q, angle and dimension d."],
   text["Why this benchmark: it is the first test of the LO pipeline (design/LOPIPELINE.md): couplings that are functions of p^2, a loop \
at generic external momentum whose frequency integral has two pole families (at q and at Q = p + q), the anomalous dimensions from \
f(0) = 1 as the p -> 0 limit of the same loop, and the dimensionless flow at fixed p/k."],
   input["Get[\"~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl\"];\nkpzloDir = ExpandFileName[\"~/Physics/Sources/FRGscalarMSR/benchmarks/KPZLO\"];"],

   section["1. Model and ansatz"],
   text["Gamma_k = Int hbar { dt h - nu fnu(-Laplacian) lap h - lam/2 (grad h)^2 } - D hbar fD(-Laplacian) hbar, so \
Gamma(1,1) = i w + nu p^2 fnu(p^2) and Gamma(0,2) = -2 D fD(p^2). fn[f, X] is f(-Laplacian) acting on X; \"Functions\" ties each \
function to its reference coupling and states f(0) = 1. The fn in the density selects the LO pipeline."],
   input["kpzloTheory = DefineTheory[\"Scalars\" -> {h}, \"Dimension\" -> Dim];"],
   input["kpzloAnsatz = DefineAnsatz[kpzloTheory,\n  \"Density\" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],\n  \"Couplings\" -> {nuk, lamk, Dk}, \"Functions\" -> {fD -> Dk, fnu -> nuk}];"],

   section["2. Symmetries, vertices, regulator, propagator, loops"],
   text["In LO the symmetries are declarations: the Galilean Ward identity (the 1 of i w) is checked at p = 0, the shift symmetry by the \
vanishing of the p^0 order of Gamma(1,1) (so the fnu row can be divided by p^2). The regulators are those of kpz_notation.tex with two \
independent shapes, R_nu = nu q^2 rn and R_D = D rd; the notes use one shape r, substituted in the comparison."],
   input["kpzloSymmetries = DeclareSymmetries[kpzloAnsatz,\n  Protected[{hbar, h}, \"Frequency\" -> 1, \"Momentum\" -> 0, \"Reason\" -> \"Galilean Ward identity (the 1 in front of dt h)\"],\n  Protected[{h, hbar}, \"Frequency\" -> 0, \"Momentum\" -> 0, \"Reason\" -> \"h -> h + c (shift symmetry)\"],\n  NonRenormalised[lamk, \"Reason\" -> \"Galilean invariance\"]];"],
   input["kpzloVertices = BuildVertices[kpzloAnsatz];"],
   input["kpzloRegulator = DefineRegulator[kpzloAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], \"Shapes\" -> {rn, rd}];"],
   input["kpzloPropagator = BuildPropagator[kpzloVertices, kpzloRegulator];"],
   input["kpzloLoops = BuildLoopIntegrands[kpzloVertices, kpzloPropagator, kpzloRegulator];"],

   section["3. LO equations, oracle, dimensionless form"],
   text["BuildLOEquations integrates the loop over the frequency by residues at the poles of Delta(w, q^2) and Delta(w, Q^2), divides \
each Gamma2 entry by its prefactor (nuk p^2 for fnu, -2 Dk for fD) and takes the p -> 0 limit for eta. RunOracle compares the residue \
result with a numerical frequency integral at a random point. AdimensionalEquations writes F(p^2) = f(p^2/k^2) and the flow at fixed \
p/k, which gives the linear part eta f + 2 p^2 f'(p^2)."],
   input["kpzloEquations = BuildLOEquations[kpzloLoops, kpzloSymmetries];"],
   input["RunOracle[kpzloEquations, \"ZeroFlows\" -> True];"],
   input["kpzloDims = DeriveDimensions[kpzloAnsatz, \"References\" -> {nuk -> \"nu\", Dk -> \"D\"}];\nkpzloAdim = AdimensionalEquations[kpzloEquations, kpzloDims];"],
   input["{kpzloAdim[\"Functions\"], kpzloAdim[\"Etas\"], kpzloAdim[\"Linear\"], kpzloAdim[\"BetaL\"], kpzloAdim[\"P2Divided\"]}"],

   section["4. The reference (LO_flow_equations.tex, kpz_notation.tex)"],
   text["Reference.wl transcribes the kernels A_a, A_c, A_d, the loop integrands of I_D and I_nu (with I_nu's 1/p^2), the eta \
matrix and vector, and d_s g = g (d - 2 - eta_D + 3 eta_nu), in the notes' variables ph, qh, ch = cos(theta), l_q = q^2 (fnu + r), \
k_q = fD + r, g = lam^2 D/nu^3 k^(d-2)."],
   input["Get[FileNameJoin[{kpzloDir, \"Reference.wl\"}]];\nColumn[{\"D\" -> Short[kpzloRefLoop[\"D\"], 4], \"nu\" -> Short[kpzloRefLoop[\"nu\"], 4]}, Spacings -> 2]"],

   section["5. Comparison"],
   text["Checks.wl writes the library's stage in the notes' variables (ap -> ph, aq -> qh, CosTheta -> ch, rn = rd = r, alamk^2 = g) and \
compares with ExactlyZeroQ: random rationals for every symbol and every function atom f(x), f'(x), r(x), r'(x), so the test is pointwise, \
not after integration, and holds at symbolic d."],
   input["Get[FileNameJoin[{kpzloDir, \"Checks.wl\"}]];\n{ExactlyZeroQ[kpzloLibrary[\"Loop\"][[1]] - kpzloRefLoop[\"D\"]], ExactlyZeroQ[kpzloLibrary[\"Loop\"][[2]] - kpzloRefLoop[\"nu\"]]}"],
   input["{ExactlyZeroQ[kpzloLibrary[\"EtaMatrix\"] - kpzloRefEta[\"Matrix\"]], ExactlyZeroQ[kpzloLibrary[\"EtaVector\"] - kpzloRefEta[\"Vector\"]]}"],

   section["6. Summary and the Python export"],
   text["The same checks as benchmarks/KPZLO/run.wl and the suite hook tests/acceptance/KPZLO.wl. The next cell rebuilds the stage from \
Pipeline.wl (cards hidden) and checks that it gives the same loops, so that this notebook and the scripted benchmark cannot drift apart. \
The last cell writes the BMWflow module (docs/handoff/BMWflow_LO.md) to a temporary directory."],
   input["kpzloResults = kpzloChecks[];\nkpzloSummary[kpzloResults]"],
   input["If[AllTrue[kpzloResults, #Pass &], \"all checks pass\", Abort[]]"],
   input["nbLoops = kpzloAdim[\"Loop\"];\nBlock[{kpzloAdim, Print}, Get[FileNameJoin[{kpzloDir, \"Pipeline.wl\"}]]; If[kpzloAdim[\"Loop\"] === nbLoops, \"same stage as Pipeline.wl\", Abort[]]]"],
   input["ExportPython[kpzloAdim, \"G\" -> alamk^2, \"Model\" -> \"KPZ\", \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZ_LO.py\"}]]"]};

Export[FileNameJoin[{dir, "KPZLO.nb"}], Notebook[cells, StyleDefinitions -> "Default.nb"]];
