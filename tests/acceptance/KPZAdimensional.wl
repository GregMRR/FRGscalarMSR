(* ::Package:: *)

(* Adimensional acceptance on KPZ: the dimensionless output (betas, eta system, loop measure S_d aq^(d-1)/(2 Pi)^d) reproduces
   L. Canet, arXiv:2509.01472, Eq. (20), exactly for symbolic d. Canet's g = v_d k^(d-2) lambda^2 D/nu^3 with
   v_d = 1/(2^(d-1) Pi^(d/2) Gamma(d/2)) = S_d/(2 Pi)^d, so g = v_d alamk^2: this also fixes the (2 Pi)^-d of the measure. *)

Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];

kpzDims = quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]];
kpzAdim = quiet[AdimensionalEquations[kpzEquations, kpzDims]];

(* Integral_0^1 daq Measure f(aq^2) with the Litim shape r(y) = 1/y - 1 (y < 1): (v_d/2) Integral_0^1 y^(d/2-1) f(y) dy, termwise *)
adimLitim = {rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};
adimIntegral[f_] := Module[{g, terms},
   litimSmoothGuard[f /. Normal[kpzAdim["Atoms"]]];
   g = Expand[Together[(f /. Normal[kpzAdim["Atoms"]] /. aq -> Sqrt[yy] /. adimLitim) yy^(Dim/2 - 1)]];
   terms = If[Head[g] === Plus, List @@ g, {g}];
   If[! AllTrue[terms, FreeQ[Together[#/yy^Exponent[#, yy]], yy] &], $Failed,
    vd/2 Total[(# /. yy -> 1)/(Exponent[#, yy] + 1) & /@ terms]]];

test["KPZ adim: state and etas", {kpzAdim["Couplings"], kpzAdim["Etas"], kpzAdim["MustFix"]}, {{alamk}, {etanu, etaD}, {}}];
test["KPZ adim: measure = v_d aq^(d-1)", Together[kpzAdim["Measure"]/(aq^(Dim - 1)) - 1/(2^(Dim - 1) Pi^(Dim/2) Gamma[Dim/2])], 0];
test["KPZ adim: linear part of alamk", Expand[First[kpzAdim["BetaL"]]], Expand[alamk (Dim/2 - 1 + 3 etanu/2 - etaD/2)]];
kpzEtaSol = Module[{mnl = Map[adimIntegral, kpzAdim["EtaMNL"], {2}], bnl = adimIntegral /@ kpzAdim["EtaBNL"]},
   Solve[Thread[(kpzAdim["EtaML"] + mnl) . kpzAdim["Etas"] == kpzAdim["EtaBL"] + bnl], kpzAdim["Etas"]]];
test["KPZ adim: eta system solved uniquely", Length[kpzEtaSol], 1];
kpzAdimBracket = Together[2 First[kpzAdim["BetaL"]]/alamk /. First[kpzEtaSol] /. alamk -> Sqrt[gt/vd] /. Dim -> dd];
canetAdim = (2 dd^6 + 4 dd^5 - dd^4 (7 gt + 8) - 2 dd^3 (3 gt + 8) + dd^2 gt (3 gt + 20) - dd (gt - 8) gt - 2 gt^2)/
   (2 dd^5 + 8 dd^4 + dd^3 (8 - 3 gt) - 2 dd^2 gt + dd gt (gt + 8) - gt^2);
test["KPZ adim: d_s g / g = Canet Eq. (20), exactly, any d", Together[kpzAdimBracket - canetAdim], 0];
test["KPZ adim: d = 1 fixed point g* = 6/5", Together[kpzAdimBracket /. {dd -> 1, gt -> 6/5}], 0];

(* the exported file evaluated in Julia (and loaded by FlowRG) *)
Get[FileNameJoin[{testRoot, "reference", "JuliaRoundTrip.wl"}]];
Scan[juliaRoundTrip["KPZ" <> ToString[#] <> "D", kpzAdim, #] &, {1, 2, 3}];
