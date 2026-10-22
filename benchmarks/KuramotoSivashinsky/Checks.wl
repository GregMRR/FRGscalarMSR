(* ::Package:: *)

(* Checks of the KS benchmark: the library's flows (Pipeline.wl) against Gosteva, Wschebor, Canet, arXiv:2607.15784 (Reference.wl).
   Needs FRGscalarMSR, Pipeline.wl and Reference.wl loaded. Every comparison is exact: the difference must vanish at every loop
   momentum, tested with ExactlyZeroQ (random rationals for every symbol and every shape atom r(x), r'(x), ...).
   ksChecks[] gives one record per check: <|"Name", "Paper", "Pass", "Detail"|>. No test-harness dependency. *)

(* the library's d^dq/(2 Pi)^d integrand of d_s(comb), in the paper's variables (Eq. 25) at dimension d:
   q = k qh, nu = D Fh, tau = D Th/k^2, zeta = D Ch/k^2, lam^2 = Gh k D^2, d_s D = -eta D; rn -> rnu, rd -> rD.
   At d = 1, Int dq/(2 Pi) I(q) = Int dqh/(2 Pi) k I(k qh), so the paper's J = k^(1 + 2 j) I(k qh)/D, with j = 1 for the
   k^2 d_s zeta/D and k^2 d_s tau/D equations. *)
ksToPaper[d_] := {nuk -> Dk Fh, tauk -> Dk Th/k^2, zetak -> Dk Ch/k^2, lamk -> Sqrt[Gh k] Dk, dsCoupling[Dk] -> -eta Dk, Dim -> d};
ksLibraryFlow[comb_, kPower_, d_ : 1] := Together[k^kPower/Dk (EquationOf[ksEquations, comb]["Integrand"] /. ksToPaper[d] /.
       q -> k qh /. {rn -> rnu, rd -> rD})];
ksLibrary[c_String, d_ : 1] := Switch[c, "D", ksLibraryFlow[Dk, 1, d], "nu", ksLibraryFlow[nuk, 1, d],
   "zeta", ksLibraryFlow[zetak, 3, d], "tau", ksLibraryFlow[tauk, 3, d]];

(* the regulator choice (32) and the time-reversal point (19) *)
ksAlphaRegulator = {rnu -> (al rD[#] &)};
ksTimeReversal = {Fh -> al, Th -> al Ch};
ksThirdDerivatives[e_] := Union[Cases[e, Derivative[3][_][_], {0, Infinity}]];

ksRecord[name_, paper_, pass_, detail_] := <|"Name" -> name, "Paper" -> paper, "Pass" -> TrueQ[pass], "Detail" -> detail|>;

ksChecks[] := Module[{closureNames = ksClosure["Checks"][[All, "Check"]], eqs, tauAlpha, trs},
   eqs = Table[With[{diff = ksLibrary[c] - (ksRef[c] /. ksSRules)},
      ksRecord["KS.equation." <> c, Switch[c, "D", "Eq. (34)", "nu", "Eq. (35)", "zeta", "Eq. (36)", "tau", "Eq. (37)"],
       ExactlyZeroQ[diff], Switch[c, "D", "d_s D/D, Gamma(0,2) at p^0", "nu", "d_s nu/D, Gamma(1,1) at p^2",
        "zeta", "k^2 d_s zeta/D, Gamma(0,2) at p^2", "tau", "k^2 d_s tau/D, Gamma(1,1) at p^4 (the n = 2 projection)"] <>
        ": library - paper = 0 at every qh"]], {c, {"D", "nu", "zeta", "tau"}}];
   tauAlpha = ksLibrary["tau"] /. ksAlphaRegulator;
   trs[d_] := {(ksLibrary["nu", d] - al ksLibrary["D", d]) /. ksAlphaRegulator /. ksTimeReversal,
      (ksLibrary["tau", d] - al ksLibrary["zeta", d]) /. ksAlphaRegulator /. ksTimeReversal};
   Join[
    {ksRecord["KS.protected", "Sec. III A", MemberQ[closureNames, "Closure.protected"] && MemberQ[ksSymmetries["NonRenormalised"][[All, 1]], lamk],
      "the 1 of i w (Galilean) and the p^0 of Gamma(1,1) (shift) flow to exactly 0; lam declared non-renormalised"]},
    eqs,
    {ksRecord["KS.thirdDerivatives", "below Eq. (41)",
      ksThirdDerivatives[ksLibrary["tau"]] =!= {} && ksThirdDerivatives[tauAlpha] === {},
      "d_s tau contains r''' for general r_nu, r_D, and r''' cancels for r_nu = alpha r_D"],
     ksRecord["KS.timeReversal.nu", "Eqs. (19), (32)", ExactlyZeroQ[trs[1][[1]]],
      "d = 1, r_nu = alpha r_D, nu = alpha D, tau = alpha zeta: d_s nu - alpha d_s D = 0 at every qh"],
     ksRecord["KS.timeReversal.tau", "Eqs. (19), (32)", ExactlyZeroQ[trs[1][[2]]],
      "same point: d_s tau - alpha d_s zeta = 0 at every qh (p^4 projection against p^2 projection)"],
     ksRecord["KS.timeReversal.controlD2", "control", ! ExactlyZeroQ[trs[2][[1]]] && ! ExactlyZeroQ[trs[2][[2]]],
      "the same two combinations are nonzero at d = 2, where time reversal is absent"]}]];

ksSummary[checks_] := Grid[Prepend[{#["Name"], #["Paper"], If[#["Pass"], "PASS", "FAIL"], #["Detail"]} & /@ checks,
    {"check", "paper", "result", "what"}], Alignment -> Left, Frame -> All];
