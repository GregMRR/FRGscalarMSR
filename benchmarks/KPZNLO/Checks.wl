(* ::Package:: *)

(* Checks of the KPZ NLO benchmark (design/NLOPIPELINE.md). References: BMWflow's NLO kernels with f_lambda slots (Reference.wl), the
   LO benchmarks at w = 0 (KPZLO, KPZLOLambda), KCW Eq. 39 and Eq. 56, the numerical w integral (Oracle). Exact: ExactlyZeroQ at
   symbolic d. Needs kpzlo*, kpzlam* (LO pipelines), kpznlo*, kpzlamnlo* (Pipeline.wl) and Reference.wl.
   kpznloChecks[] gives <|"Name", "Paper", "Pass", "Detail"|> records. *)

kpznloToRef = {ap -> ph, aq -> qh, CosTheta -> ch, rn -> r, rd -> r, alamk -> Sqrt[g], aw -> w};
kpznloOne = {flam[_] -> 1, Derivative[_][flam][_] -> 0};
kpznloRecord[name_, paper_, pass_, detail_] := <|"Name" -> name, "Paper" -> paper, "Pass" -> TrueQ[pass], "Detail" -> detail|>;
kpznloCardOK[stage_, names_] := SubsetQ[stage["Checks"][[All, "Check"]], names] && AllTrue[stage["Checks"], #Status === "PASS" &];

kpznloChecks[] := Module[{oracle},
   oracle = Block[{Print}, RunOracle[kpzlamnloEquations, "ZeroFlows" -> True]];
   {
    kpznloRecord["KPZNLO.contract", "design/NLOPIPELINE.md N2, N6",
     {kpzlamnloAdim["Approximation"], kpzlamnloAdim["Functions"], kpzlamnloAdim["Etas"], kpzlamnloAdim["P2Divided"], kpznloAdim["Functions"]} ===
      {"NLO", {fD, fnu, flam}, {etaD, etanu}, {0, 1, 0}, {fD, fnu}},
     "NLO stages; FUNCTIONS (fD, fnu, flam) and (fD, fnu), ETAS and P2_DIVIDED as LO"],
    kpznloRecord["KPZNLO.reference.D", "Kloss et al. 2012 NLO kernels with f_lambda slots (BMWflow model_kpz.py)",
     ExactlyZeroQ[(kpzlamnloAdim["Loop"][[1]] /. kpznloToRef) - kpzlamnloRefLoop["D"]],
     "generic f_lambda, w: the loop integrand of fD at every p, q, c, w, d"],
    kpznloRecord["KPZNLO.reference.nu", "Kloss et al. 2012 NLO kernels with f_lambda slots (BMWflow model_kpz.py)",
     ExactlyZeroQ[(kpzlamnloAdim["Loop"][[2]] /. kpznloToRef) - kpzlamnloRefLoop["nu"]],
     "generic f_lambda, w: the loop integrand of fnu (divided by p^2)"],
    kpznloRecord["KPZNLO.flam1", "BMWflow NLO = KCW LO", ExactlyZeroQ[(Take[kpzlamnloAdim["Loop"], 2] /. kpznloOne) - kpznloAdim["Loop"]],
     "f_lambda = 1: the fD, fnu rows of the f_lambda model equal the KPZ NLO rows"],
    kpznloRecord["KPZNLO.slice", "the w = 0 slice is closed", ExactlyZeroQ[(kpzlamnloAdim["Loop"] /. aw -> 0) - kpzlamAdim["Loop"]] &&
       ExactlyZeroQ[(kpznloAdim["Loop"] /. aw -> 0) - kpzloAdim["Loop"]] &&
       ExactlyZeroQ[{kpzlamnloAdim["EtaMatrix"] - kpzlamAdim["EtaMatrix"], kpzlamnloAdim["EtaVector"] - kpzlamAdim["EtaVector"], kpzlamnloAdim["BetaL"] - kpzlamAdim["BetaL"]}],
     "at aw = 0 every loop row equals the LO row (KPZLO, KPZLOLambda; f_lambda row: O/(i w) -> d/d(i w)); eta system and beta unchanged"],
    kpznloRecord["KPZNLO.linear", "KCW Eq. 56 (eta_lambda = 0)",
     ExactlyZeroQ[kpzlamnloAdim["Linear"] - {etaD fD[ap^2, aw] + 2 ap^2 Derivative[1, 0][fD][ap^2, aw] + (2 - etanu) aw Derivative[0, 1][fD][ap^2, aw],
        etanu fnu[ap^2, aw] + 2 ap^2 Derivative[1, 0][fnu][ap^2, aw] + (2 - etanu) aw Derivative[0, 1][fnu][ap^2, aw],
        2 ap^2 Derivative[1, 0][flam][ap^2, aw] + (2 - etanu) aw Derivative[0, 1][flam][ap^2, aw]}],
     "d_s f_X = eta_X f_X + p d_p f_X + (2 - eta_nu) w d_w f_X, with 2 - eta_nu derived from the dimension of the frequency"],
    kpznloRecord["KPZNLO.flam.ward", "KCW Eq. 39: f_lambda(w, 0) = 1", ExactlyZeroQ[kpzlamnloAdim["Loop"][[3]] /. ap -> 0],
     "the f_lambda row vanishes at p = 0 for every w"],
    kpznloRecord["KPZNLO.flam.flows", "design/NLOPIPELINE.md N4", ! ExactlyZeroQ[D[kpzlamnloAdim["Loop"][[3]], aw]],
     "the f_lambda row depends on w (not the LO row copied to every w)"],
    kpznloRecord["KPZNLO.symmetries", "design/NLOPIPELINE.md N7",
     kpznloCardOK[kpzlamnloEquations, {"ExternalMomentum.nloLimit", "ExternalMomentum.nloNoiseEven", "ExternalMomentum.nloWardAllFrequencies",
        "ExternalMomentum.nloDividedOrders", "ExternalMomentum.nloZeroEntries"}] &&
      kpznloCardOK[kpznloEquations, {"ExternalMomentum.nloLimit", "ExternalMomentum.nloNoiseEven"}] &&
      MemberQ[kpznloEquations["Warnings"][[All, "Check"]], "ExternalMomentum.wardNeglectedNLO"],
     "NLO checks pass: rows at w = 0 = LO, Gamma(hbar, hbar) even in w, f_lambda row 0 at p = 0 and the p^0 order of fnu's row 0 at every w, Gamma(h, h) = 0; \
f_lambda = 1: the odd part of Gamma(hbar, h) dropped (INFO)"],
    kpznloRecord["KPZNLO.oracle", "numerical w integral (Oracle, NLO)", AllTrue[oracle["Checks"], #Status === "PASS" &] && Length[oracle["Results"]] == 24,
     "even / odd parts at w = 3/10, 7/5, 4 = NIntegrate of the numeric loop (relative 1e-8), d = 2, 3: worst " <>
      ToString[CForm[N[Max[oracle["Results"][[All, "Error"]]], 2]]]]}];

kpznloSummary[checks_] := Grid[Prepend[{#["Name"], #["Paper"], If[#["Pass"], "PASS", "FAIL"], #["Detail"]} & /@ checks,
    {"check", "reference", "result", "what"}], Alignment -> Left, Frame -> All];
