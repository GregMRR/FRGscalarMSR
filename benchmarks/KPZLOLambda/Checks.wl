(* ::Package:: *)

(* Checks of the KPZ LO + f_lambda benchmark (design/LOFLAMBDA.md). The reference is the KPZLO benchmark (f_lambda = 1, itself checked
   against BMWflow's LO notes): with f_lambda -> 1, f_lambda' -> 0 every equation of fD and fnu (loop integrands at generic p, eta rows,
   linear parts, coupling flow) equals KPZLO's, exactly (ExactlyZeroQ: random rationals for every symbol and function atom), at
   symbolic d, in the dimensionful Equations stage and in the dimensionless one. Needs kpzloEquations/kpzloAdim (KPZLO Pipeline.wl)
   and kpzlamEquations/kpzlamAdim (Pipeline.wl here). kpzlamChecks[] gives <|"Name", "Paper", "Pass", "Detail"|> records. *)

kpzlamOne = {flam[_] -> 1, Derivative[_][flam][_] -> 0};
kpzlamRow[stage_, f_] := SelectFirst[stage["Equations"], #Function === f &];
kpzlamRecord[name_, paper_, pass_, detail_] := <|"Name" -> name, "Paper" -> paper, "Pass" -> TrueQ[pass], "Detail" -> detail|>;
(* the dimensionless rows of fD, fnu (positions 1, 2 in both stages: FUNCTIONS order) *)
kpzlamFD := Take[kpzlamAdim["Functions"], 2] === kpzloAdim["Functions"];

(* generic f_lambda: the library in the reference variables of Reference.wl (as benchmarks/KPZLO/Checks.wl) *)
kpzlamToRef = {ap -> ph, aq -> qh, CosTheta -> ch, rn -> r, rd -> r, alamk -> Sqrt[g]};

(* DE route: the (p^2)^n coefficient of a dimensionful LO integrand f(p, q, c) at fD = fnu = 1, f_lambda = 1 + xk s, averaged over the
   direction of p in d dimensions (<c^k> = (k-1)!!/(d (d+2)...)), against the DE flow of the coupling c *)
kpzlamPoly = {fD -> (1 &), fnu -> (1 &), flam -> (1 + xk # &)};
kpzlamCAvg[k_] := If[OddQ[k], 0, (k - 1)!!/Product[Dim + 2 j, {j, 0, k/2 - 1}]];
kpzlamPCoef[e_, n_] := Expand[SeriesCoefficient[e /. kpzlamPoly, {p, 0, 2 n}]] /. CosTheta^kk_. :> kpzlamCAvg[kk];
kpzlamDE[c_] := SelectFirst[kpzlamDEEquations["Equations"], #Flowing === c &]["Integrand"] /. LoopMomentumSquared -> q^2;

kpzlamChecks[] := Module[{eqNames = kpzlamEquations["Checks"][[All, "Check"]], oracle},
   oracle = Block[{Print}, RunOracle[kpzlamEquations]];
   {
    kpzlamRecord["KPZLOLambda.contract", "design/LOFLAMBDA.md F2, F5",
     {kpzlamAdim["Functions"], kpzlamAdim["HasEta"], kpzlamAdim["Etas"], kpzlamAdim["Couplings"], kpzlamAdim["P2Divided"]} ===
      {{fD, fnu, flam}, {True, True, False}, {etaD, etanu}, {alamk}, {0, 1, 0}},
     "FUNCTIONS (fD, fnu, flam), ETAS (etaD, etanu, None), COUPLINGS (alamk,), P2_DIVIDED (0, 1, 0)"],
    kpzlamRecord["KPZLOLambda.flam1.equations", "KPZLO, Equations stage",
     AllTrue[{fD, fnu}, ExactlyZeroQ[(kpzlamRow[kpzlamEquations, #]["Integrand"] /. kpzlamOne) - kpzlamRow[kpzloEquations, #]["Integrand"]] &&
        ExactlyZeroQ[(kpzlamRow[kpzlamEquations, #]["EtaIntegrand"] /. kpzlamOne) - kpzlamRow[kpzloEquations, #]["EtaIntegrand"]] &],
     "f_lambda = 1: the dimensionful loop integrands of d_s(Dk fD), d_s(nuk fnu) and their p -> 0 limits equal KPZLO's, at every p, q, c, d"],
    kpzlamRecord["KPZLOLambda.flam1.loop", "KPZLO, dimensionless", kpzlamFD && ExactlyZeroQ[(Take[kpzlamAdim["Loop"], 2] /. kpzlamOne) - kpzloAdim["Loop"]],
     "f_lambda = 1: the dimensionless loop integrands of fD and fnu (fnu divided by p^2) equal KPZLO's"],
    kpzlamRecord["KPZLOLambda.flam1.eta", "KPZLO, dimensionless",
     ExactlyZeroQ[(kpzlamAdim["EtaMatrix"] /. kpzlamOne) - kpzloAdim["EtaMatrix"]] && ExactlyZeroQ[(kpzlamAdim["EtaVector"] /. kpzlamOne) - kpzloAdim["EtaVector"]],
     "f_lambda = 1: the eta matrix and vector (etaD, etanu) equal KPZLO's"],
    kpzlamRecord["KPZLOLambda.flam1.linear", "KPZLO", kpzlamFD && ExactlyZeroQ[Take[kpzlamAdim["Linear"], 2] - kpzloAdim["Linear"]],
     "linear parts eta_X f_X + p f_X' of fD, fnu unchanged"],
    kpzlamRecord["KPZLOLambda.coupling", "KPZLO; Galilean invariance", ExactlyZeroQ[kpzlamAdim["BetaL"] - kpzloAdim["BetaL"]],
     "beta of alamk unchanged: lambda is not renormalised and f_lambda(0) = 1, so d_s g = g (d - 2 - eta_D + 3 eta_nu)"],
    kpzlamRecord["KPZLOLambda.reference.D", "Kloss et al. 2012 kernels with f_lambda slots (BMWflow model_kpz.py)",
     ExactlyZeroQ[(kpzlamAdim["Loop"][[1]] /. kpzlamToRef) - kpzlamRefLoop["D"]],
     "generic f_lambda: the loop integrand of fD, with f_lambda at p (squared), q and Q, at every p, q, c, d"],
    kpzlamRecord["KPZLOLambda.reference.nu", "Kloss et al. 2012 kernels with f_lambda slots (BMWflow model_kpz.py)",
     ExactlyZeroQ[(kpzlamAdim["Loop"][[2]] /. kpzlamToRef) - kpzlamRefLoop["nu"]],
     "generic f_lambda: the loop integrand of fnu (divided by p^2), with f_lambda at p, q and Q"],
    kpzlamRecord["KPZLOLambda.de.flam", "derivative expansion, f_lambda = 1 + x p^2",
     ExactlyZeroQ[kpzlamPCoef[kpzlamRow[kpzlamEquations, flam]["Integrand"], 1] - kpzlamDE[xk]] &&
      ! ExactlyZeroQ[kpzlamPCoef[-kpzlamRow[kpzlamEquations, flam]["Integrand"], 1] - kpzlamDE[xk]],
     "p^2 coefficient (direction-averaged) of the f_lambda row = the DE flow of x (coefficient of i w p^2), and its sign flip is rejected: \
sign and normalisation of the d/d(i nu) row by an independent route"],
    kpzlamRecord["KPZLOLambda.de.nuD", "derivative expansion, f_lambda = 1 + x p^2",
     ExactlyZeroQ[kpzlamPCoef[kpzlamRow[kpzlamEquations, fnu]["Integrand"], 0] - kpzlamDE[nuk]] &&
      ExactlyZeroQ[kpzlamPCoef[kpzlamRow[kpzlamEquations, fD]["Integrand"], 0] - kpzlamDE[Dk]],
     "p^0 coefficients of the fnu row (divided by p^2) and of the fD row = the DE flows of nuk and Dk at f_lambda = 1 + x p^2"],
    kpzlamRecord["KPZLOLambda.flam.linear", "design/LOFLAMBDA.md F2", ExactlyZeroQ[kpzlamAdim["Linear"][[3]] - 2 ap^2 Derivative[1][flam][ap^2]],
     "f_lambda has no eta: d_s f_lambda = p f_lambda' + I_lambda"],
    kpzlamRecord["KPZLOLambda.flam.flows", "design/LOFLAMBDA.md F3", ! kpzlamRow[kpzlamEquations, flam]["IdenticallyZero"],
     "the d/d(i nu) row of f_lambda is not identically zero at p != 0"],
    kpzlamRecord["KPZLOLambda.symmetries", "Galilean and shift symmetries",
     SubsetQ[eqNames, {"ExternalMomentum.ward", "ExternalMomentum.unitReference", "ExternalMomentum.dividedOrders", "ExternalMomentum.zeroEntries",
       "ExternalMomentum.etaAtZero", "ExternalMomentum.firstDerivativesOnly"}],
     "I_lambda(p = 0) = 0 exactly (f_lambda(0) stays 1), the p^0 order of fnu's row vanishes, Gamma(h,h) stays 0, eta rows hold f_lambda(0) = 1 and f, f', r, r' at q only"],
    kpzlamRecord["KPZLOLambda.oracle", "numerical w integral (Oracle LO)", AllTrue[oracle["Checks"], #Status === "PASS" &],
     "the residue rows (f_lambda: double poles at Q) equal NIntegrate of the numeric loop (exact d/d nu) at a random point, d = 2, 3: " <>
      StringRiffle[(#Label <> " (d = " <> ToString[#Dimension] <> ") " <> ToString[CForm[N[#Error, 2]]]) & /@ oracle["Results"], ", "]]}];

kpzlamSummary[checks_] := Grid[Prepend[{#["Name"], #["Paper"], If[#["Pass"], "PASS", "FAIL"], #["Detail"]} & /@ checks,
    {"check", "reference", "result", "what"}], Alignment -> Left, Frame -> All];
