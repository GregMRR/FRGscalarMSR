(* ::Package:: *)

(* Checks of the KPZ LO benchmark: the library's dimensionless LO stage (Pipeline.wl) against BMWflow's LO_flow_equations.tex /
   kpz_notation.tex (Reference.wl). Every comparison is exact at symbolic d: the difference must vanish at every (p, q, c), tested
   with ExactlyZeroQ (random rationals for every symbol and for every atom f(x), f'(x), r(x), r'(x)).
   kpzloChecks[] gives one record per check: <|"Name", "Paper", "Pass", "Detail"|>. No test-harness dependency. *)

(* library -> reference variables: ap -> p, aq -> q, cos -> c, both shapes -> the single r of the notes, alamk = Sqrt[g] *)
kpzloToRef = {ap -> ph, aq -> qh, CosTheta -> ch, rn -> r, rd -> r, alamk -> Sqrt[g]};
kpzloLibrary[key_] := kpzloAdim[key] /. kpzloToRef;
kpzloLibraryLinear := kpzloLibrary["Linear"] /. {Derivative[1][fD][ph^2] -> fDprime/(2 ph), Derivative[1][fnu][ph^2] -> fnuprime/(2 ph)};
(* BMWflow's g = alamk^2: d_s g = 2 alamk d_s alamk *)
kpzloLibraryGFlow := Together[2 alamk First[kpzloAdim["BetaL"]] /. alamk -> Sqrt[g]];

kpzloRecord[name_, paper_, pass_, detail_] := <|"Name" -> name, "Paper" -> paper, "Pass" -> TrueQ[pass], "Detail" -> detail|>;

kpzloChecks[] := Module[{eqNames = kpzloEquations["Checks"][[All, "Check"]], loop = kpzloLibrary["Loop"], eta},
   {
    kpzloRecord["KPZLO.contract", "design/LOPIPELINE.md section 6",
     {kpzloAdim["Functions"], kpzloAdim["Etas"], kpzloAdim["Couplings"], kpzloAdim["P2Divided"]} === {{fD, fnu}, {etaD, etanu}, {alamk}, {0, 1}},
     "FUNCTIONS (fD, fnu), ETAS (etaD, etanu), COUPLINGS (alamk,), P2_DIVIDED (0, 1)"],
    kpzloRecord["KPZLO.loop.D", "I_D, LO_flow_equations.tex", ExactlyZeroQ[loop[[1]] - kpzloRefLoop["D"]],
     "loop integrand of d_s f_D = 2g (q^2 + pqc)^2 k_Q [A_a d_sR_D - 2 q^2 l_q k_q A_c d_sR_nu], at every p, q, c, d"],
    kpzloRecord["KPZLO.loop.nu", "I_nu, LO_flow_equations.tex", ExactlyZeroQ[loop[[2]] - kpzloRefLoop["nu"]],
     "loop integrand of d_s f_nu, already divided by p^2 (the shift symmetry makes its p^0 order vanish), at every p, q, c, d"],
    kpzloRecord["KPZLO.eta.matrix", "I_DD, I_Dnu, I_nuD, I_nunu, kpz_notation.tex", ExactlyZeroQ[kpzloLibrary["EtaMatrix"] - kpzloRefEta["Matrix"]],
     "eta matrix (row X = equation of f_X, column Y multiplies eta_Y): the p -> 0 limit of the loop, at every q and d"],
    kpzloRecord["KPZLO.eta.vector", "I_D0, I_nu0, kpz_notation.tex", ExactlyZeroQ[kpzloLibrary["EtaVector"] - kpzloRefEta["Vector"]],
     "eta vector, at every q and d"],
    kpzloRecord["KPZLO.linear", "LO_flow_equations.tex Eq. 1", ExactlyZeroQ[kpzloLibraryLinear - Values[kpzloRefLinear]],
     "linear part eta_X f_X + p f_X'(p) (flow at fixed p/kappa)"],
    kpzloRecord["KPZLO.coupling", "kpz_notation.tex, model_kpz.py g_rhs_calc", ExactlyZeroQ[kpzloLibraryGFlow - kpzloRefGFlow],
     "d_s g = 2 alamk beta_alamk = g (d - 2 - eta_D + 3 eta_nu), g = alamk^2"],
    kpzloRecord["KPZLO.symmetries", "Galilean and shift symmetries",
     SubsetQ[eqNames, {"ExternalMomentum.ward", "ExternalMomentum.dividedOrders", "ExternalMomentum.zeroEntries"}],
     "the 1 of i w flows to 0 at p = 0, the p^0 order of Gamma(1,1) vanishes (so f_nu's row is divided by p^2), Gamma(h,h) stays 0"]}];

kpzloSummary[checks_] := Grid[Prepend[{#["Name"], #["Paper"], If[#["Pass"], "PASS", "FAIL"], #["Detail"]} & /@ checks,
    {"check", "reference", "result", "what"}], Alignment -> Left, Frame -> All];
