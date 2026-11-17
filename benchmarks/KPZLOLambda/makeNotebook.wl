(* ::Package:: *)

(* makeNotebook.wl -- writes KPZLOLambda.nb (this folder) from the cell list below.
   Run:  math -script ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/makeNotebook.wl
   Check that it runs headless:
         math -script ~/sources/FRGscalarMSR/examples/runNotebook.wl ~/sources/FRGscalarMSR/benchmarks/KPZLOLambda/KPZLOLambda.nb
   Edit the cells here, not in the .nb. The pipeline cells mirror Pipeline.wl; section 8 checks that they build the same stage. *)

dir = DirectoryName[$InputFileName];
title[s_] := Cell[s, "Title"];
section[s_] := Cell[s, "Section"];
text[s_] := Cell[s, "Text"];
input[s_] := Cell[s, "Input"];

cells = {
   title["Benchmark: KPZ at LO with a flowing f_lambda"],
   text["Model: the Kloss-Canet-Wschebor LO ansatz (PRE 86, 051124 (2012)) with three functions of p^2. BMWflow's hand-written \
KPZ model sets f_lambda = 1; here f_lambda flows. This notebook derives the flow of f_D, f_nu and f_lambda with FRGscalarMSR from \
the ansatz alone, checks it three independent ways, and exports it to BMWflow (contract v3, docs/handoff/BMWflow_LO.md), where it \
is integrated to the fixed point in d = 1, 2, 3."],
   text["Why this benchmark: it is the first LO model with a function times i w and inside a vertex (design/LOFLAMBDA.md). That needs \
(i) the flow of an i w coefficient, d/d(i nu) of d_s Gamma(1,1)(nu, p) taken before the frequency integral, which gives double poles; \
(ii) a function evaluated at the external momentum p as well as at q and Q = p + q; (iii) a function with reference 1, whose \
f(0) = 1 is a Ward identity, so it has no eta."],
   input["Get[\"~/sources/FRGscalarMSR/FRGscalarMSR.wl\"];\nkpzlamDir = ExpandFileName[\"~/sources/FRGscalarMSR/benchmarks/KPZLOLambda\"];"],

   section["1. Model and ansatz"],
   text["Gamma_k = Int hbar f_lambda(-Laplacian) { dt h - lam/2 (grad h)^2 } - nu hbar f_nu(-Laplacian) lap h - D hbar f_D(-Laplacian) hbar, \
so Gamma(1,1) = i w f_lambda(p^2) + nu p^2 f_nu(p^2), Gamma(0,2) = -2 D f_D(p^2), and the cubic vertex is lam f_lambda(k_hbar^2) (k1.k2), \
with f_lambda at the momentum of its hbar leg. Galilean invariance ties the f_lambda of the vertex to the f_lambda of dt h, which is \
why f_lambda multiplies the whole combination. fn[f, X] is f(-Laplacian) acting on X."],
   text["\"Functions\" -> {fD -> Dk, fnu -> nuk, flam -> 1}: f_D and f_nu are normalised by their couplings (f(0) = 1 defines D and nu, \
and gives eta_D, eta_nu); f_lambda has reference 1, i.e. no coupling and no eta."],
   input["kpzlamTheory = DefineTheory[\"Scalars\" -> {h}, \"Dimension\" -> Dim];"],
   input["kpzlamAnsatz = DefineAnsatz[kpzlamTheory,\n  \"Density\" -> fn[flam, hbar] (dt[h] - lamk/2 grad[h].grad[h]) - nuk hbar fn[fnu, lap[h]] - Dk hbar fn[fD, hbar],\n  \"Couplings\" -> {nuk, lamk, Dk}, \"Functions\" -> {fD -> Dk, fnu -> nuk, flam -> 1}];"],

   section["2. Symmetries, vertices, regulator, propagator, loops"],
   text["The Galilean Ward identity says the i w coefficient of Gamma(1,1) at p = 0 is not renormalised: f_lambda(0) stays 1. It is \
declared as a Protected i w tag and checked at p = 0 when the LO equations are built (ExternalMomentum.ward). The shift symmetry \
removes the p^0 order of Gamma(1,1), so the f_nu row can be divided by p^2. lam is not renormalised (Galilean invariance)."],
   input["kpzlamSymmetries = DeclareSymmetries[kpzlamAnsatz,\n  Protected[{hbar, h}, \"Frequency\" -> 1, \"Momentum\" -> 0, \"Reason\" -> \"Galilean Ward identity: f_lambda(0) = 1\"],\n  Protected[{h, hbar}, \"Frequency\" -> 0, \"Momentum\" -> 0, \"Reason\" -> \"h -> h + c (shift symmetry)\"],\n  NonRenormalised[lamk, \"Reason\" -> \"Galilean invariance\"]];"],
   input["kpzlamVertices = BuildVertices[kpzlamAnsatz];"],
   text["Regulators as in benchmarks/KPZLO: R_nu = nu q^2 rn(q^2/k^2) and R_D = -2 D rd(q^2/k^2), two independent shapes. The \
propagator now carries i w f_lambda(q^2): Delta(w) = c0 + c1 (i w) with c1 = f_lambda. The residues assume f_lambda > 0 at every q \
and Q (causal poles; stated on the ExternalMomentum card, not checked at run time)."],
   input["kpzlamRegulator = DefineRegulator[kpzlamAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], \"Shapes\" -> {rn, rd}];"],
   input["kpzlamPropagator = BuildPropagator[kpzlamVertices, kpzlamRegulator];"],
   input["kpzlamLoops = BuildLoopIntegrands[kpzlamVertices, kpzlamPropagator, kpzlamRegulator];"],

   section["3. LO equations, oracle, dimensionless form"],
   text["BuildLOEquations integrates each loop over the frequency by residues. The rows of f_D and f_nu are the w^0 coefficients of \
Gamma(0,2) and Gamma(1,1), as in KPZLO. The row of f_lambda is the i nu coefficient of d_s Gamma(1,1)(nu, p) at nu = 0, with d/d nu taken \
before the w integral: nu enters only G(Q), so the poles at Q become double poles. RunOracle compares every row with a numerical w \
integral (with an exact d/d nu) at a random point, in d = 2 and 3."],
   input["kpzlamEquations = BuildLOEquations[kpzlamLoops, kpzlamSymmetries];"],
   input["RunOracle[kpzlamEquations];"],
   text["AdimensionalEquations writes F(p^2) = f(p^2/k^2) at fixed p/k. f_D and f_nu get eta_X f_X + p f_X' as linear part; f_lambda, \
which has no eta, gets p f_lambda' only. The eta system has one row per eta (2 x 2), not one per function."],
   input["kpzlamDims = DeriveDimensions[kpzlamAnsatz, \"References\" -> {nuk -> \"nu\", Dk -> \"D\"}];\nkpzlamAdim = AdimensionalEquations[kpzlamEquations, kpzlamDims];"],
   input["{kpzlamAdim[\"Functions\"], kpzlamAdim[\"HasEta\"], kpzlamAdim[\"Etas\"], kpzlamAdim[\"P2Divided\"]}"],
   input["{kpzlamAdim[\"Linear\"], kpzlamAdim[\"BetaL\"], Dimensions[kpzlamAdim[\"EtaMatrix\"]]}"],

   section["4. The f_lambda row"],
   text["The new row, in the dimensionless variables (ap = p/k, aq = q/k, CosTheta = cos of the angle between p and q). It contains \
f_lambda at the external momentum (flam[ap^2]) as well as at q and Q. Its size compared with the two KPZ rows:"],
   input["LeafCount /@ AssociationThread[kpzlamAdim[\"Functions\"], kpzlamAdim[\"Loop\"]]"],
   input["Short[kpzlamAdim[\"Loop\"][[3]], 6]"],
   text["The Ward identity, read off the row itself: at p = 0 the f_lambda integrand vanishes identically (at every q, angle, d and \
every value of the functions), so f_lambda(0) = 1 is kept by the flow. The library checks this when it builds the equations; these \
are the checks it ran on this model:"],
   input["ExactlyZeroQ[kpzlamAdim[\"Loop\"][[3]] /. ap -> 0]"],
   input["Select[kpzlamEquations[\"Checks\"][[All, \"Check\"]], StringStartsQ[#, \"ExternalMomentum.\"] &]"],

   section["5. Reference 1: f_lambda = 1 gives KPZLO"],
   text["benchmarks/KPZLO derives the f_lambda = 1 model and checks it exactly against BMWflow's LO notes. Setting f_lambda -> 1, \
f_lambda' -> 0 here must give exactly the KPZLO rows of f_D and f_nu, their eta system and beta. ExactlyZeroQ substitutes random \
rationals for every symbol and every function atom, so the comparison is pointwise (before any integral) and holds at symbolic d."],
   input["Block[{Print}, Get[FileNameJoin[{ParentDirectory[kpzlamDir], \"KPZLO\", \"Pipeline.wl\"}]]];\nkpzlamOne = {flam[_] -> 1, Derivative[_][flam][_] -> 0};"],
   input["{ExactlyZeroQ[(kpzlamAdim[\"Loop\"][[1]] /. kpzlamOne) - kpzloAdim[\"Loop\"][[1]]], ExactlyZeroQ[(kpzlamAdim[\"Loop\"][[2]] /. kpzlamOne) - kpzloAdim[\"Loop\"][[2]]]}"],
   input["{ExactlyZeroQ[(kpzlamAdim[\"EtaMatrix\"] /. kpzlamOne) - kpzloAdim[\"EtaMatrix\"]], ExactlyZeroQ[(kpzlamAdim[\"EtaVector\"] /. kpzlamOne) - kpzloAdim[\"EtaVector\"]], ExactlyZeroQ[kpzlamAdim[\"BetaL\"] - kpzloAdim[\"BetaL\"]]}"],

   section["6. Reference 2: the Kloss et al. kernels at generic f_lambda"],
   text["BMWflow's flow/models/model_kpz.py keeps the Kloss et al. LO kernels with f_lambda slots at p, q and Q (set to 1 there). \
Reference.wl transcribes them in the variables of the notes (ph, qh, ch, one shape r, g = lam^2 D/nu^3 k^(d-2)). The library's f_D \
and f_nu rows must equal them at generic f_lambda; this tests f_lambda at the external momentum away from p = 0."],
   input["Get[FileNameJoin[{kpzlamDir, \"Reference.wl\"}]];\nkpzlamToRef = {ap -> ph, aq -> qh, CosTheta -> ch, rn -> r, rd -> r, alamk -> Sqrt[g]};\nColumn[{\"D\" -> Short[kpzlamRefLoop[\"D\"], 3], \"nu\" -> Short[kpzlamRefLoop[\"nu\"], 3]}, Spacings -> 2]"],
   input["{ExactlyZeroQ[(kpzlamAdim[\"Loop\"][[1]] /. kpzlamToRef) - kpzlamRefLoop[\"D\"]], ExactlyZeroQ[(kpzlamAdim[\"Loop\"][[2]] /. kpzlamToRef) - kpzlamRefLoop[\"nu\"]]}"],

   section["7. Reference 3: the derivative expansion at f_lambda = 1 + x p^2"],
   text["No source gives the f_lambda row. Its sign and normalisation are checked by a different route: the derivative expansion of \
the same model with f_lambda(p^2) = 1 + x p^2 and f_nu = f_D = 1, i.e. (hbar - x lap hbar)(dt h - lam/2 (grad h)^2) - nu hbar lap h - D hbar^2, \
which uses the DE projection, single-family residues and Assemble's normalisation. The flow of x must equal the p^2 coefficient of \
the f_lambda row (averaged over the direction of p), and the flows of nu and D the p^0 coefficients of the other two rows. Pipeline.wl \
builds the DE stage (kpzlamDEEquations); Checks.wl holds the Taylor-coefficient helpers."],
   input["Block[{Print}, Get[FileNameJoin[{kpzlamDir, \"Pipeline.wl\"}]]];\nGet[FileNameJoin[{kpzlamDir, \"Checks.wl\"}]];\nkpzlamFlamRow = kpzlamRow[kpzlamEquations, flam][\"Integrand\"];"],
   input["{ExactlyZeroQ[kpzlamPCoef[kpzlamFlamRow, 1] - kpzlamDE[xk]], ExactlyZeroQ[kpzlamPCoef[-kpzlamFlamRow, 1] - kpzlamDE[xk]]}"],
   text["True, False: the p^2 coefficient matches the DE flow of x, and the same row with its sign flipped does not."],

   section["8. Summary and the Python export"],
   text["The same checks as run.wl and the suite hook tests/acceptance/KPZLOLambda.wl. The next cell rebuilds the stage from \
Pipeline.wl (cards hidden) and checks that it gives the same loops, so that this notebook and the scripted benchmark cannot drift apart. \
The last cell writes the BMWflow module (contract v3) to a temporary directory; export.wl writes it to out/ for runBMWflow.py."],
   input["kpzlamResults = kpzlamChecks[];\nkpzlamSummary[kpzlamResults]"],
   input["If[AllTrue[kpzlamResults, #Pass &], \"all checks pass\", Abort[]]"],
   input["nbLoops = kpzlamAdim[\"Loop\"];\nBlock[{kpzlamAdim, Print}, Get[FileNameJoin[{kpzlamDir, \"Pipeline.wl\"}]]; If[kpzlamAdim[\"Loop\"] === nbLoops, \"same stage as Pipeline.wl\", Abort[]]]"],
   input["ExportPython[kpzlamAdim, \"G\" -> alamk^2, \"Model\" -> \"KPZLambda\", \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZLambda_LO.py\"}]]"],

   section["9. In BMWflow: exponents at the fixed point"],
   text["runBMWflow.py integrates both modules (f_lambda = 1 and flowing f_lambda) with BMWflow's ModelGenerated, Euler steps \
ds = 1e-3 down to s = -20, grids and regulator of BMWflow's KPZ_1D_twogrids example; z = 2 - eta_nu, chi = (2 - d + eta_D - eta_nu)/2. \
Literature: d = 1 exact; d = 2 Pagnani & Parisi, PRE 92, 010101(R) (2015); d = 3 Marinari, Pagnani & Parisi, J. Phys. A 33, 8181 (2000) \
(chi measured, z = 2 - chi). Results of 2026-11-03 (docs/handoff/BMWflow_LO.md, Status (v3)):"],
   input["Grid[{{\"d\", \"model\", \"z\", \"chi\"},\n  {1, \"f_lambda = 1\", 1.50000, 0.50000}, {1, \"flowing f_lambda\", 1.45004, 0.54997}, {1, \"exact\", 3/2, 1/2},\n  {2, \"f_lambda = 1\", 1.67116, 0.32884}, {2, \"flowing f_lambda\", 1.61277, 0.38723}, {2, \"simulations\", \"1.6131(4)\", \"0.3869(4)\"},\n  {3, \"f_lambda = 1\", 1.82762, 0.17237}, {3, \"flowing f_lambda\", 1.82058, 0.17941}, {3, \"simulations\", \"1.6865(15)\", \"0.3135(15)\"}},\n Frame -> All, Alignment -> Left]"],
   text["In d = 1 a flowing f_lambda breaks the time-reversal symmetry of KPZ (Int h'' f_lambda(-d^2) (h')^2 != 0 for f_lambda != 1), \
so its exponents are not protected there (design/LOFLAMBDA.md section 4). The plots below exist after running \
~/sources/BMWflow/.venv/bin/python benchmarks/KPZLOLambda/runBMWflow.py: eta(s) and the functions f_X(p) at s = -20."],
   input["With[{files = FileNames[\"kpz_flambda_d*.png\", FileNameJoin[{kpzlamDir, \"out\"}]]},\n If[files === {}, \"no plots yet: run runBMWflow.py\", Column[Import /@ files]]]"],
   text["Large p (runBMWflow.py --p-max 1e5 --n-p 175, three more decades at the same density; the *_pmax100000_np175 plots above). \
f_lambda falls from 1 around p ~ 1-10 and then stays constant through p = 1e5: f_lambda(infinity) ~ 0.541, 0.650, 0.786 in d = 1, 2, 3 \
(in d = 2, 3 after a shallow minimum near p ~ 5-10). f_D and f_nu fall as p^(-eta_D), p^(-eta_nu) of the same run; f_lambda has no eta \
and tends to a constant. The exponents move by at most 1e-3 from the p_max = 100 values (chi with f_lambda: d = 2 0.38794, d = 3 0.18008)."],
   input["Grid[Prepend[{{1, 1, 0.9066, 0.6156, 0.5445, 0.5415, 0.5414, 0.5415, 0.5424}, {2, 1, 0.9451, 0.7019, 0.6472, 0.6495, 0.6498, 0.6498, 0.6497},\n   {3, 1, 0.9791, 0.7992, 0.7761, 0.7859, 0.7862, 0.7862, 0.7862}},\n  {\"d \\\\ p\", 0, 1, 3, 10, 10^2, 10^3, 10^4, 10^5}], Frame -> All]"]};

Export[FileNameJoin[{dir, "KPZLOLambda.nb"}], Notebook[cells, StyleDefinitions -> "Default.nb"]];
