(* ::Package:: *)

(* Milestone 6 acceptance: KPZ, the only cross-check of the loop normalisation (E3/E4).
   L. Canet, arXiv:2509.01472, Eq. (20) (also LPA_KPZ_thermal_forcing/lpa/kpz.jl): for the ansatz of the bare action with the Litim
   cutoff r(y) = (1/y - 1) theta(1 - y) in both sectors, the dimensionless coupling g = v_d k^(d-2) lambda^2 D / nu^3 flows as
     d_s g = g (d - 2 - eta_D + 3 eta_nu) = g num(d, g)/den(d, g),   d_s nu = -eta_nu nu,  d_s D = -eta_D D.
   We integrate our dimensionful nu and D equations over q with the Litim shapes, solve for eta_nu, eta_D and compare the
   bracket exactly, for symbolic d. This tests the bubble prefactor, sign and 2 Pi, the frequency integrals and the momentum
   projection; it does not test the Gamma4 tadpole (absent in KPZ). *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

litim = {rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};
(* Int d^dq/(2 Pi)^d f(q^2) = v_d k^d/2 Int_0^1 y^(d/2 - 1) f(k^2 y) dy (the Litim integrand vanishes for y > 1, its r' has no delta term
   here because it multiplies r(1) = 0); every term is c y^a with a depending on d *)
litimIntegral[f_] := Module[{g = Expand[Together[f /. q -> k Sqrt[yy] /. litim] yy^(dd/2 - 1)], terms},
   terms = If[Head[g] === Plus, List @@ g, {g}];
   If[! AllTrue[terms, FreeQ[# /. yy -> 1, yy] && FreeQ[Together[#/yy^Exponent[#, yy]], yy] &], $Failed,
    vd k^dd/2 Total[(# /. yy -> 1)/(Exponent[#, yy] + 1) & /@ terms]]];
kpzFlow[c_] := litimIntegral[EquationOf[kpzEquations, c]["Integrand"] /. Dim -> dd];
etaRules = {dsCoupling[nuk] -> -etanu nuk, dsCoupling[Dk] -> -etaD Dk};
etaSol = Solve[{-etanu nuk == (kpzFlow[nuk] /. etaRules), -etaD Dk == (kpzFlow[Dk] /. etaRules)}, {etanu, etaD}];
bracket = Together[(dd - 2 - etaD + 3 etanu) /. First[etaSol] /. lamk -> Sqrt[gt nuk^3/(vd k^(dd - 2) Dk)]];
canet = (2 dd^6 + 4 dd^5 - dd^4 (7 gt + 8) - 2 dd^3 (3 gt + 8) + dd^2 gt (3 gt + 20) - dd (gt - 8) gt - 2 gt^2)/
   (2 dd^5 + 8 dd^4 + dd^3 (8 - 3 gt) - 2 dd^2 gt + dd gt (gt + 8) - gt^2);

test["KPZ Litim: q-integrals done exactly", FreeQ[{kpzFlow[nuk], kpzFlow[Dk]}, $Failed], True];
test["KPZ Litim: eta_nu, eta_D solved uniquely", Length[etaSol], 1];
test["KPZ Litim: bracket depends only on d and g (k, nu, D, v_d drop out)", Complement[Variables[bracket], {dd, gt}], {}];
test["KPZ Litim: d_s g / g = Canet Eq. (20), exactly, any d", Together[bracket - canet], 0];
test["KPZ Litim: d = 1 fixed point g* = 6/5 (al* = Sqrt[6/5] in kpz.jl)", Together[bracket /. {dd -> 1, gt -> 6/5}], 0];
test["KPZ: Galilean Ward identity (i w in Gamma_{hbar h}) flows to exactly 0", MemberQ[kpzClosure["Checks"][[All, "Check"]], "Closure.protected"], True];
