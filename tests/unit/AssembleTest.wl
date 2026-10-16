(* ::Package:: *)

(* Assemble: scale stripping, KPZ equations (contract), oracle on KPZ including the flows that must vanish *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

strip = FRGscalarMSR`Private`stripScale;
test["scale stripping: nk s (Ak + r[s])", strip[nk s Ak + nk s r[s], s, {nk, Ak}, {r}], {nk, 1, Ak + r[s]}];
test["scale stripping: -2 D (a + D2 s + r[s])", strip[-2 D a - 2 D D2 s - 2 D r[s], s, {D, a, D2}, {r}], {-2 D, 0, a + D2 s + r[s]}];
test["scale stripping: single term H l s", strip[H l s, s, {H, l}, {r}], {H l, 1, 1}];

eqD = EquationOf[kpzEquations, Dk];
test["KPZ: 2 equations", Length[kpzEquations["Equations"]], 2];
test["KPZ: explicit integrand has no intermediate symbols",
  FreeQ[kpzEquations["Equations"][[All, "Integrand"]], GammaAtom | ScaleAtom | DeltaCoefficient | dsR | LoopMomentumSquared | \[CapitalDelta]], True];
test["KPZ: d_s Dk depends on dsCoupling[nuk] and dsCoupling[Dk]", eqD["DsCouplings"], {dsCoupling[Dk], dsCoupling[nuk]}];
test["KPZ: integrand linear in dsCoupling", ExactlyZeroQ[D[eqD["Integrand"], {dsCoupling[Dk], 2}]], True];
test["KPZ: compact form /. atoms /. d_s R rules = integrand",
  ExactlyZeroQ[(eqD["Compact"] /. Normal[kpzEquations["Atoms"]] /. kpzEquations["DsRRules"]) - eqD["Integrand"]], True];
test["KPZ: shapes kept abstract in q^2/k^2", Union[Cases[eqD["Integrand"], (rn | rd)[x_] | Derivative[_][rn | rd][x_] :> x, Infinity]], {q^2/k^2}];

kpzOracle = quiet[RunOracle[kpzEquations, "ZeroFlows" -> True]];
test["KPZ oracle (d = 2, 3), equations and every flow that must vanish", Length[kpzOracle["Results"]], 20];
test["KPZ oracle max error below 1e-8", Max[kpzOracle["Results"][[All, "Error"]]] < 10^-8, True];
