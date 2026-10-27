(* ::Package:: *)

(* makeNotebooks.wl -- writes examples/ShallowWaterDE4plus.nb, examples/KPZ.nb and examples/KPZLO.nb from the cell lists below.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/examples/makeNotebooks.wl
   Check that they run headless:  math -script ~/Physics/Sources/FRGscalarMSR/examples/runNotebook.wl <notebook.nb>
   Edit the cells here, not in the .nb files, then regenerate. *)

dir = DirectoryName[$InputFileName];
title[s_] := Cell[s, "Title"];
section[s_] := Cell[s, "Section"];
text[s_] := Cell[s, "Text"];
input[s_] := Cell[s, "Input"];
nb[cells_] := Notebook[cells, StyleDefinitions -> "Default.nb"];

load = input["Get[\"~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl\"]"];

shallowWater = {
   title["Shallow water, DE4plus ansatz"],
   text["This notebook derives the 9 flow equations of the shallow-water DE4plus ansatz (two scalars H, Psi with response fields \
Hbar, Psibar), the case DE4plus_fast.wl was written for. Every step returns an immutable stage and prints a summary card: one PASS line \
per check, the warnings and the assumptions. A failed check prints a red FAIL block and stops the evaluation. Full expressions are \
shown only on request, with the Show... commands. The modules are described in docs/MODULES.md, the conventions in docs/CONVENTIONS.md."],
   load,
   section["1. Fields and ansatz"],
   text["DefineTheory creates a response field <name>bar for each scalar. The ansatz is the action density, written with the \
package's own derivatives dt, grad, div, lap (not the System Dt, Grad: Dt[nk H] would add Dt[nk] terms). The parser prints the table \
of parsed monomials so that you can check it."],
   input["theory = DefineTheory[\"Scalars\" -> {H, Psi}, \"Dimension\" -> Dim]"],
   input["ansatz = DefineAnsatz[theory, \"Density\" ->\n    Hbar (dt[H] - nk Ak lap[H] + Hk lk lap[Psi] + lk div[H grad[Psi]])\n  + Psibar (dt[Psi] - nk ank lap[Psi] + gk H - g2k lap[H] + lk/2 grad[Psi].grad[Psi])\n  - Dhk aDhk grad[Hbar].grad[Hbar] - Dphik aDphik Psibar^2 - Dphik D2phik grad[Psibar].grad[Psibar] + Kk grad[Hbar].grad[Psibar],\n  \"Couplings\" -> {nk, Ak, Hk, lk, ank, gk, g2k, Dhk, aDhk, Dphik, aDphik, D2phik, Kk}]"],
   section["2. Vertices, regulator, propagator"],
   text["BuildVertices prints every Fourier vertex under the convention header (dt -> +I w, grad -> -I q, all momenta incoming, \
last leg = minus the sum of the others). The regulator entries are the amounts added to Gamma2_ab, including the 2 of a diagonal pair; \
the shapes stay abstract functions of y = q^2/k^2. The propagator is inverted in abstract atoms; Delta is the determinant whose roots \
are the poles, and its validity conditions say when they are causal (upper half w-plane)."],
   input["vertices = BuildVertices[ansatz];"],
   input["regulator = DefineRegulator[ansatz,\n  {Hbar, H} -> nk q^2 arnh[q^2/k^2], {Psibar, Psi} -> nk q^2 arnphi[q^2/k^2],\n  {Hbar, Hbar} -> -2 Dhk q^2 armh[q^2/k^2], {Psibar, Psibar} -> -2 Dphik armphi[q^2/k^2],\n  \"Shapes\" -> {arnh, arnphi, armh, armphi}];"],
   input["prop = BuildPropagator[vertices, regulator];"],
   section["3. Symmetries: what flows, what must not, what is dropped"],
   text["Every monomial (I w)^m (q^2)^n of every Gamma2 entry up to the truncation orders gets a class: declared in the ansatz \
(it gets an equation), Protected (a symmetry forbids it, so its flow is checked to vanish), Neglected (computed, reported, dropped), \
or zero by structure. An untagged monomial stops the evaluation. Ward identities are assumed through the Protected tags, not derived."],
   input["symmetries = DeclareSymmetries[ansatz,\n  Protected[{Hbar, _}, \"Momentum\" -> 0, \"Reason\" -> \"Hbar -> Hbar + c(t) (mass conservation)\"],\n  Protected[{Psi, _}, \"Momentum\" -> 0, \"Reason\" -> \"Psi -> Psi + c(t)\"],\n  Neglected[{H, Psibar}, \"Frequency\" -> {1, 2}, \"Momentum\" -> 0, \"Reason\" -> \"truncation, z = 1 power counting\"],\n  Neglected[{Psibar, Psibar}, \"Frequency\" -> 2, \"Momentum\" -> 0, \"Reason\" -> \"coloured forcing, dropped\"],\n  Neglected[_, \"Frequency\" -> {1, 2}, \"Momentum\" -> 1, \"Reason\" -> \"mixed i w q^2 / w^2 q^2 terms, dropped\"],\n  NonRenormalised[lk, \"Reason\" -> \"Galilean invariance\"],\n  \"TruncationOrders\" -> <|\"Frequency\" -> 2, \"Momentum\" -> 1|>];"],
   section["4. Loops, projection, frequency integrals, closure"],
   text["The one-loop integrand of d_s Gamma2 is projected on each monomial at zero external frequency and momentum, then \
integrated over the loop frequency exactly by residues (no square roots: the result is rational in the coefficients of Delta). \
CheckClosure asserts the vanishing flows and turns each declared monomial into d_s(combination) = RHS, the number in front of the \
coefficient (e.g. the -2 of -2 Dhk aDhk) being divided out."],
   input["loops = BuildLoopIntegrands[vertices, prop, regulator];"],
   input["projected = ProjectOnMonomials[loops, symmetries];"],
   input["integrated = IntegrateFrequencies[projected];"],
   input["flows = CheckClosure[integrated, symmetries];"],
   input["ShowClosureReport[flows]"],
   section["5. Equations"],
   text["AssembleEquations writes each equation compactly in scale-stripped atoms (e.g. ScaleAtom[H, Hbar, 0, 0] = Ak + arnh[q^2/k^2]) \
and explicitly: the integrand of Int d^dq/(2 Pi)^d, with dsCoupling[c] = k dc/dk entering linearly and the shapes and their derivatives \
abstract. It takes about 35 s."],
   input["equations = AssembleEquations[flows];"],
   input["ShowEquation[equations, gk]"],
   input["ShowEquation[equations, Dhk aDhk]"],
   section["6. Numerical oracle"],
   text["RunOracle recomputes every equation numerically at one random point, without the propagator algebra or the residues: matrix \
inverse of Gamma2 + R, numerical frequency integral, finite differences in the external frequency and momentum, in d = 2 and d = 3. \
It shares the vertices and the loop formula with the symbolic chain, so it does not test the loop normalisation; the KPZ notebook does."],
   input["oracle = RunOracle[equations];"],
   section["7. Saving"],
   input["file = SaveStage[equations, FileNameJoin[{$TemporaryDirectory, \"swDE4plus-equations.wdx\"}]]"],
   input["equations2 = LoadStage[file, \"Expect\" -> {theory, ansatz, regulator, symmetries}];\nequations2 === equations"],
   section["8. Dimensionless flows and the FlowRG file"],
   text["DeriveDimensions derives the dimension of every coupling from the ansatz (Gamma_k is dimensionless), in the basis of k and \
N+1 = 3 reference couplings, and checks them (independent references, consistent monomials, action invariant under the derived \
rescaling). CheckDimensions compares hand values exactly. Here every reference appears only times an amplitude (nk ank, Dhk aDhk, \
Dphik aDphik); \"Mode\" -> \"Running\" makes each eta follow from that equation (the amplitude is counted as 1), the convention \
of the old hand-built export. Without it the eta would be fixed at run time (MUST_FIX)."],
   input["dims = DeriveDimensions[ansatz, \"References\" -> {Dphik -> {\"phi\", \"Mode\" -> \"Running\"}, Dhk -> {\"h\", \"Mode\" -> \"Running\"},\n    nk -> {\"nu\", \"Mode\" -> \"Running\", \"Equation\" -> nk ank}}, \"Dimensionless\" -> {Ak, ank, aDhk, aDphik}];"],
   input["CheckDimensions[dims, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, lk -> k^(1 - Dim/2) Dphik^(-1/2) nk^(3/2), \"Frequency\" -> k^2 nk|>]"],
   text["AdimensionalEquations writes d_s chat = BetaL (not integrated) + Integral daq LoopMeasure BetaNL for each state coupling, and \
the linear system of the etas. Every beta is checked to be free of k and the references, structurally and by a random rescaling of the \
units. ExportJulia writes the FlowRG equations file for one d (this one uses r'', D2<shape>, which FlowRG must bind)."],
   input["adim = AdimensionalEquations[equations, dims, \"State\" -> {gk, g2k, lk, Hk, aDphik, D2phik, Kk, Ak, ank, aDhk}];\nShowBeta[adim, gk]"],
   input["ExportJulia[adim, \"Dimension\" -> 2, \"File\" -> FileNameJoin[{$TemporaryDirectory, \"swDE4plus\", \"Equations2D.jl\"}]]"]};

kpz = {
   title["KPZ, one loop, dimensionful"],
   text["The KPZ equation dt h = nu lap h + lambda/2 (grad h)^2 + noise, with noise variance 2 D, written as an MSR action. \
The ansatz has the form of the bare action. The Galilean symmetry protects the coefficient 1 of dt h (Ward identity) and the \
non-linearity lambda (not checkable in v1: that needs d_s Gamma3, so it is declared NonRenormalised)."],
   load,
   input["theory = DefineTheory[\"Scalars\" -> {h}, \"Dimension\" -> Dim];\nansatz = DefineAnsatz[theory, \"Density\" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2,\n  \"Couplings\" -> {nuk, lamk, Dk}];"],
   input["symmetries = DeclareSymmetries[ansatz,\n  Protected[{hbar, h}, \"Frequency\" -> 1, \"Momentum\" -> 0, \"Reason\" -> \"Galilean Ward identity\"],\n  Protected[{h, hbar}, \"Frequency\" -> 0, \"Momentum\" -> 0, \"Reason\" -> \"h -> h + c (shift symmetry)\"],\n  Neglected[_, \"Reason\" -> \"beyond the ansatz\"], NonRenormalised[lamk, \"Reason\" -> \"Galilean invariance\"],\n  \"TruncationOrders\" -> <|\"Frequency\" -> 2, \"Momentum\" -> 1|>];"],
   input["vertices = BuildVertices[ansatz];\nregulator = DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], \"Shapes\" -> {rn, rd}];\nprop = BuildPropagator[vertices, regulator];"],
   input["flows = CheckClosure[IntegrateFrequencies[ProjectOnMonomials[BuildLoopIntegrands[vertices, prop, regulator], symmetries]], symmetries];"],
   input["equations = AssembleEquations[flows];\nShowEquation[equations, nuk]\nShowEquation[equations, Dk]"],
   input["oracle = RunOracle[equations, \"ZeroFlows\" -> True];"],
   section["Check of the loop normalisation against Canet's flow"],
   text["L. Canet, arXiv:2509.01472, Eq. (20): with the Litim cutoff r(y) = (1/y - 1) theta(1 - y) in both sectors, the dimensionless \
coupling g = v_d k^(d-2) lambda^2 D/nu^3 flows as d_s g = g (d - 2 - eta_D + 3 eta_nu), with d_s nu = -eta_nu nu, d_s D = -eta_D D. \
We integrate our two equations over q with the Litim shape (Int d^dq/(2 Pi)^d f = v_d k^d/2 Int_0^1 y^(d/2-1) f dy), solve the \
linear system for eta_nu and eta_D, and compare. The difference is exactly 0 for every d: this checks the prefactor, sign and 2 Pi \
of the loop formula (not the Gamma4 tadpole, absent here)."],
   input["litim = {rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};\nlitimIntegral[f_] := Module[{g = Expand[Together[f /. q -> k Sqrt[y] /. litim] y^(d/2 - 1)], terms},\n  terms = If[Head[g] === Plus, List @@ g, {g}];\n  vd k^d/2 Total[(# /. y -> 1)/(Exponent[#, y] + 1) & /@ terms]];\nflowOf[c_] := litimIntegral[EquationOf[equations, c][\"Integrand\"] /. Dim -> d];"],
   input["etas = First[Solve[{-etanu nuk == (flowOf[nuk] /. {dsCoupling[nuk] -> -etanu nuk, dsCoupling[Dk] -> -etaD Dk}),\n     -etaD Dk == (flowOf[Dk] /. {dsCoupling[nuk] -> -etanu nuk, dsCoupling[Dk] -> -etaD Dk})}, {etanu, etaD}]];\nbracket = Together[(d - 2 - etaD + 3 etanu) /. etas /. lamk -> Sqrt[g nuk^3/(vd k^(d - 2) Dk)]]"],
   input["canet = (2 d^6 + 4 d^5 - d^4 (7 g + 8) - 2 d^3 (3 g + 8) + d^2 g (3 g + 20) - d (g - 8) g - 2 g^2)/\n   (2 d^5 + 8 d^4 + d^3 (8 - 3 g) - 2 d^2 g + d g (g + 8) - g^2);\nTogether[bracket - canet]"],
   section["Dimensionless flow and the FlowRG file"],
   text["The same, done by the library: dimensions from the ansatz in the basis (k, nuk, Dk), the dimensionless flow of lamk (linear part \
only: lambda is NonRenormalised) and the eta system, then the FlowRG file for d = 2. tests/acceptance/KPZAdimensional.wl checks that this \
output gives Canet's Eq. (20) exactly."],
   input["dims = DeriveDimensions[ansatz, \"References\" -> {nuk -> \"nu\", Dk -> \"D\"}];"],
   input["adim = AdimensionalEquations[equations, dims];\nShowBeta[adim, lamk]\nShowEtaSystem[adim]"],
   input["ExportJulia[adim, \"Dimension\" -> 2, \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZ\", \"Equations2D.jl\"}]]"]};

kpzlo = {
   title["KPZ at LO: couplings that are functions of p^2"],
   text["The second pipeline of the library (design/LOPIPELINE.md). The Gamma2 couplings of KPZ become functions of the momentum: \
Gamma^(1,1) = i w + nu p^2 f_nu(p^2), Gamma^(0,2) = -2 D f_D(p^2), with f_nu(0) = f_D(0) = 1, the vertex lambda constant. \
The library derives the flow of each function at generic p (one loop, two pole families: Delta at q and at Q = p + q), the \
anomalous dimensions from f(0) = 1, the dimensionless form at fixed p/k, and writes a Python module for BMWflow. \
The exact comparison with BMWflow's notes (LO_flow_equations.tex, kpz_notation.tex) is the benchmark benchmarks/KPZLO; this \
notebook shows the workflow and checks f = 1 against Canet's Eq. (20)."],
   load,
   section["1. Ansatz with functions"],
   text["fn[f, X] in the density is f(-Laplacian) acting on X: it multiplies the term by f(p^2), p the momentum of X's leg. \
\"Functions\" -> {f -> reference} states f(0) = 1 and names the k-dependent amplitude; the ansatz infers \"Pipeline\" -> \"LO\" from \
the fn. The symmetries take the same tags as in DE but no truncation orders: the Galilean Ward identity protects the 1 of i w, the \
shift symmetry h -> h + c makes the p^0 order of Gamma^(1,1) vanish."],
   input["theory = DefineTheory[\"Scalars\" -> {h}, \"Dimension\" -> Dim];\nansatz = DefineAnsatz[theory,\n  \"Density\" -> hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar],\n  \"Couplings\" -> {nuk, lamk, Dk}, \"Functions\" -> {fD -> Dk, fnu -> nuk}];"],
   input["symmetries = DeclareSymmetries[ansatz,\n  Protected[{hbar, h}, \"Frequency\" -> 1, \"Momentum\" -> 0, \"Reason\" -> \"Galilean Ward identity\"],\n  Protected[{h, hbar}, \"Frequency\" -> 0, \"Momentum\" -> 0, \"Reason\" -> \"h -> h + c (shift symmetry)\"],\n  NonRenormalised[lamk, \"Reason\" -> \"Galilean invariance\"]];"],
   section["2. Vertices, regulator, propagator, loops"],
   text["The same stages as in DE. The propagator atoms now contain fnu[s] and fD[s]; the regulator shapes are functions of q^2/k^2."],
   input["vertices = BuildVertices[ansatz];\nregulator = DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], \"Shapes\" -> {rn, rd}];\nprop = BuildPropagator[vertices, regulator];\nloops = BuildLoopIntegrands[vertices, prop, regulator];"],
   section["3. Flows at generic p"],
   text["BuildLOEquations replaces Projection ... Assemble. For each function it takes the loop at external momentum p and zero \
external frequency, integrates the loop frequency by residues at the poles of both families, and divides by the ansatz prefactor \
(nuk p^2 for fnu, -2 Dk for fD: P2Divided records the p^2). It checks that the entries outside the ansatz vanish at every p, the \
Ward identity, that the p^0 order of the fnu row vanishes, and builds the eta rows as the exact p -> 0 limit. The integrand is \
explicit in p, q, CosTheta = cos(p, q), f at q^2 and at Q^2 = p^2 + q^2 + 2 p q CosTheta, the shapes and dsCoupling."],
   input["equations = BuildLOEquations[loops, symmetries];\nShowLOEquation[equations, fD]"],
   input["oracle = RunOracle[equations, \"ZeroFlows\" -> True];"],
   section["4. Dimensionless flows"],
   text["F_k(p^2) = f(p^2/k^2), ap = p/k, aq = q/k. At fixed ap each function flows as d_s f = eta f + 2 ap^2 f'(ap^2) + Loop \
(integrated with Measure[\"Loop\"] over aq and theta). f(0) = 1 fixes the etas: (1 + <EtaMatrix>) eta = -<EtaVector>, with \
<.> = Integral daq Measure[\"Eta\"]. The constant coupling alamk = lamk Dk^(1/2) nuk^(-3/2) k^(d/2-1) has only a linear part."],
   input["dims = DeriveDimensions[ansatz, \"References\" -> {nuk -> \"nu\", Dk -> \"D\"}];\nadim = AdimensionalEquations[equations, dims];"],
   input["Thread[adim[\"Functions\"] -> adim[\"Linear\"]]"],
   input["{adim[\"Etas\"], adim[\"Couplings\"] -> adim[\"BetaL\"], adim[\"Measure\"][\"Eta\"]}"],
   input["Short[adim[\"Loop\"], 6]"],
   section["5. Check: f = 1 gives Canet's Eq. (20)"],
   text["At f_D = f_nu = 1 the eta system is the one of the DE KPZ flow (KPZ.nb). With the Litim cutoff r(y) = 1/y - 1 for y < 1 \
in both sectors, the aq integral is elementary (y = aq^2). Measure[\"Eta\"] = v_d aq^(d-1) and Canet's g = v_d alamk^2, so v_d \
drops out with alamk^2 = g/v_d: set the v_d of the measure to 1 and alamk^2 = g. The difference with Canet's d_s g/g is exactly 0 at every d."],
   input["litim = {fD -> Function[x, 1], fnu -> Function[x, 1], rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};\nlitimAvg[f_] := Module[{e = Expand[Together[(f /. litim /. aq -> Sqrt[y] /. alamk -> Sqrt[g] /. Dim -> d)] y^(d/2 - 1)/2], terms},\n  terms = If[Head[e] === Plus, List @@ e, {e}];\n  Total[(# /. y -> 1)/(Exponent[#, y] + 1) & /@ terms]];"],
   input["etas = Thread[adim[\"Etas\"] -> LinearSolve[IdentityMatrix[2] + Map[litimAvg, adim[\"EtaMatrix\"], {2}], -(litimAvg /@ adim[\"EtaVector\"])]];\nbracket = Together[2 First[adim[\"BetaL\"]]/alamk /. etas /. Dim -> d]"],
   input["canet = (2 d^6 + 4 d^5 - d^4 (7 g + 8) - 2 d^3 (3 g + 8) + d^2 g (3 g + 20) - d (g - 8) g - 2 g^2)/\n   (2 d^5 + 8 d^4 + d^3 (8 - 3 g) - 2 d^2 g + d g (g + 8) - g^2);\nTogether[bracket - canet]"],
   section["6. The BMWflow module"],
   text["ExportPython writes the module of the contract docs/handoff/BMWflow_LO.md: FUNCTIONS, ETAS, COUPLINGS, P2_DIVIDED and the \
functions linear, loop, eta_matrix, eta_vector, beta. The loop measure is not in the kernels (the adapter applies it); primes on f are \
d/d|momentum|. \"Model\" -> \"KPZ\" picks the row of PythonModelTable that maps alamk to BMWflow's g = alamk^2."],
   input["ExportPython[adim, \"Model\" -> \"KPZ\", \"File\" -> FileNameJoin[{$TemporaryDirectory, \"KPZLO\", \"KPZ_LO.py\"}]]"]};

Export[FileNameJoin[{dir, "ShallowWaterDE4plus.nb"}], nb[shallowWater]];
Export[FileNameJoin[{dir, "KPZ.nb"}], nb[kpz]];
Export[FileNameJoin[{dir, "KPZLO.nb"}], nb[kpzlo]];
Print["wrote ShallowWaterDE4plus.nb (", Count[shallowWater, Cell[_, "Input"]], " inputs), KPZ.nb (", Count[kpz, Cell[_, "Input"]], " inputs) and KPZLO.nb (", Count[kpzlo, Cell[_, "Input"]], " inputs)"];
