(* ::Package:: *)

(* Propagator.wl -- BuildPropagator: G = (Gamma2 + R)^-1 in abstract entry atoms, its common denominator Delta(w) Delta(-w),
   the coefficients of Delta and the conditions for its roots to lie in the upper half plane.
   Objective, assumptions and checks: docs/MODULES.md section "Propagator". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

BuildPropagator::usage = "BuildPropagator[vertices, regulator] inverts Gamma2 + R in abstract entry atoms and returns the \"Propagator\" stage.";
GammaAtom::usage = "GammaAtom[a, b, m][s] is the coefficient of (I w)^m in (Gamma2 + R)_ab(w) at q^2 = s, for a before b in field order \
(the (b, a) entry is the same polynomial at -w). Kept abstract until the end (primer rule 1).";
\[CapitalDelta]::usage = "\[CapitalDelta][w, s] = det of the (phi, phibar) block of Gamma2 + R at frequency w, q^2 = s: \
Delta(w) = c0 + c1 (I w) + c2 (I w)^2 with c_m = DeltaCoefficient[m][s]. Its roots are the poles of the retarded propagator.";
DeltaCoefficient::usage = "DeltaCoefficient[m][s] is the coefficient c_m of (I w)^m in Delta(w) at q^2 = s (abstract).";
PropagatorMatrix::usage = "PropagatorMatrix[prop, w, s] is G(w) at q^2 = s in field order, in GammaAtom atoms and \[CapitalDelta][+-w, s].";
FullMatrix::usage = "FullMatrix[prop, w, s] is (Gamma2 + R)(w) at q^2 = s in GammaAtom atoms.";
AtomRules::usage = "AtomRules[prop] are the rules GammaAtom[a, b, m] -> Function[s, explicit value] (couplings, shapes r[s/k^2]).";
DeltaRules::usage = "DeltaRules[prop] are the rules expressing \[CapitalDelta] and DeltaCoefficient in GammaAtom atoms.";
ValidityConditions::usage = "ValidityConditions[prop] are the inequalities (in DeltaCoefficient atoms) for all roots of Delta(w) to lie in the upper half plane.";
ShowPropagator::usage = "ShowPropagator[prop, {a, b}] prints G_ab(w, q) in atoms; ShowPropagator[prop] prints the atoms, Delta and the validity conditions.";

Begin["`Private`"];

SetAttributes[{GammaAtom, \[CapitalDelta], DeltaCoefficient}, Protected];

adjugate[m_] := If[Length[m] == 1, {{1}}, Transpose[Table[(-1)^(i + j) Det[Drop[m, {i}, {j}]], {i, Length[m]}, {j, Length[m]}]]];

BuildPropagator[vertices_, regulator_] := Block[{$FRGModule = "Propagator"},
  Module[{t0 = AbsoluteTime[], fields, theory, phis, bars, n, nf, pos, g2, rm, full, atoms = {}, atomVals = <||>, abs, iphi, ibar, B, Bm, C,
    detB, deg, cs, G, checks = {}, warns = {}, stage, valid, explanation, ex},
   validateStage[vertices, "Vertices"]; validateStage[regulator, "Regulator"];
   FRGAssert["Propagator.sameAnsatz", regulator["Upstream"] === vertices["Upstream"], {regulator["Upstream"], vertices["Upstream"]},
    "the regulator and the vertices were built from different ansatz stages"];
   theory = vertices["Ansatz"]["Theory"]; fields = vertices["Fields"]; phis = theory["Scalars"]; bars = theory["ResponseFields"];
   nf = Length[fields]; n = Length[phis]; pos[a_] := First[FirstPosition[fields, a]];

   (* explicit Gamma2 + R, polynomial in w *)
   g2 = Gamma2Matrix[vertices, ww, qv] /. scal[qv, qv] -> ss;
   rm = RegulatorMatrix[regulator, ss];
   full = Expand[g2 + rm];
   FRGAssert["Propagator.polynomial", AllTrue[Flatten[full], PolynomialQ[#, ww] &], full, "Gamma2 + R must be polynomial in w"];

   (* abstract atoms: coefficient of (I w)^m of entry (a, b), a <= b; numbers stay numbers *)
   abs = Table[Module[{i0 = Min[i, j], j0 = Max[i, j], e, a, b},
      a = fields[[i0]]; b = fields[[j0]]; e = full[[i0, j0]];
      Sum[Module[{c = Expand[Coefficient[e, ww, m] / I^m]},
        Which[c === 0, 0, NumericQ[c], c,
         True, (If[! KeyExistsQ[atomVals, GammaAtom[a, b, m]], atomVals[GammaAtom[a, b, m]] = Function[Evaluate[{ss}], Evaluate[c]]];
          GammaAtom[a, b, m][ss])] (I If[i <= j, 1, -1] ww)^m], {m, 0, Exponent[e, ww]}]], {i, nf}, {j, nf}];
   atoms = Keys[atomVals];
   AppendTo[checks, assertCheck["Propagator.atomsReal", FreeQ[(#[ss] & /@ Values[atomVals]), Complex], atomVals,
     "the coefficients of (I w)^m must be real (no explicit I in the ansatz or the regulator)", ToString[Length[atoms]] <> " real entry atoms"]];
   AppendTo[checks, assertCheck["Propagator.atomForm", ExactlyZeroQ[Flatten[abs /. (h : GammaAtom[__])[x_] :> atomVals[h][x]] - Flatten[full]] // TrueQ,
     abs, "internal: the atom form does not reproduce Gamma2 + R", "atom form = Gamma2 + R"]];

   (* block structure in the order (phi_1..phi_N | phibar_1..phibar_N) *)
   iphi = pos /@ phis; ibar = pos /@ bars;
   AppendTo[checks, assertCheck["Propagator.causalBlock", Flatten[abs[[iphi, iphi]]] === ConstantArray[0, n^2], abs[[iphi, iphi]],
     "the (phi, phi) block of Gamma2 + R must vanish (causality)", "(phi, phi) block = 0"]];
   B[w_] := abs[[iphi, ibar]] /. ww -> w;              (* rows phi at w, columns phibar *)
   C = abs[[ibar, ibar]];                              (* noise block *)
   detB = Expand[Det[B[ww]]];
   deg = Exponent[detB, ww];
   AppendTo[checks, assertCheck["Propagator.deltaDegree", 1 <= deg <= 2, deg,
     If[deg == 0, "Delta(w) does not depend on w: the frequency integrals diverge (is a dt term missing?)",
      "deg_w Delta = " <> ToString[deg] <> " > 2: not supported in v1 (future extension, docs/FUTURE.md)"],
     "deg_w Delta = " <> ToString[deg]]];
   cs = Table[Expand[Coefficient[detB, ww, m] / I^m], {m, 0, 2}];
   FRGAssert["Propagator.deltaCoefficients", Expand[Sum[cs[[m + 1]] (I ww)^m, {m, 0, 2}] - detB] === 0 && FreeQ[cs, ww], cs, "internal"];

   (* G = M^-1 with M = [[0, B(w)], [B(-w)^T, C]] in the (phi | phibar) order:
      G_{phibar phi} = B(w)^-1, G_{phi phibar} = (B(-w)^T)^-1, G_{phi phi} = -(B(-w)^T)^-1 C B(w)^-1, G_{phibar phibar} = 0 *)
   Module[{binv, dinv, gpp},
    binv = adjugate[B[ww]]/\[CapitalDelta][ww, ss];
    dinv = adjugate[Transpose[B[-ww]]]/\[CapitalDelta][-ww, ss];
    gpp = -dinv . C . binv;
    G = ConstantArray[0, {nf, nf}];
    G[[iphi, iphi]] = Map[Together, gpp, {2}]; G[[iphi, ibar]] = dinv; G[[ibar, iphi]] = binv];
   ex = {\[CapitalDelta][w_, ss] :> (detB /. ww -> w)};
   AppendTo[checks, assertCheck["Propagator.inverse",
     ExactlyZeroQ[(G /. ex) . abs - IdentityMatrix[nf]],
     G, "internal: G (Gamma2 + R) != 1", "G.(Gamma2 + R) = 1 exactly (random rationals for every atom, w, q^2)"]];
   AppendTo[checks, assertCheck["Propagator.inverseRight", ExactlyZeroQ[abs . (G /. ex) - IdentityMatrix[nf]], G, "internal", "(Gamma2 + R).G = 1"]];
   AppendTo[checks, assertCheck["Propagator.responseBlock", Flatten[G[[ibar, ibar]]] === ConstantArray[0, n^2], G, "internal", "G_{phibar phibar} = 0"]];

   (* validity: roots of Delta in w in the upper half plane <=> roots u = I w of c0 + c1 u + c2 u^2 have negative real part *)
   {valid, explanation} = If[deg == 2,
     {{DeltaCoefficient[1][ss]/DeltaCoefficient[2][ss] > 0, DeltaCoefficient[0][ss]/DeltaCoefficient[2][ss] > 0},
      "both roots of Delta(w) lie in the upper half w-plane (causal poles) iff c1/c2 > 0 and c0/c2 > 0: writing u = I w, a quadratic \
c2 u^2 + c1 u + c0 has both roots with Re u < 0 exactly when its three coefficients have the same sign"},
     {{DeltaCoefficient[0][ss]/DeltaCoefficient[1][ss] > 0}, "the root of Delta(w) = c0 + c1 (I w) is w = I c0/c1: upper half plane iff c0/c1 > 0"}];
   AppendTo[warns, infoCheck["ASSUMPTION", "Propagator.validity", displayString[valid /. ss -> Symbol["FRGscalarMSR`Propagator`Private`s"]] <> " at every loop momentum q (not checked by the library: exported as causal_poles by ExportJulia)"]];

   stage = stampStage[<|"Stage" -> "Propagator", "Pipeline" -> "DerivativeExpansion",
      "Upstream" -> {vertices["Fingerprint"], regulator["Fingerprint"]}, "Fields" -> fields, "Scale" -> regulator["Scale"],
      "Atoms" -> atoms, "AtomValues" -> atomVals,
      "M" -> Function[Evaluate[{ww, ss}], Evaluate[abs]], "G" -> Function[Evaluate[{ww, ss}], Evaluate[G]],
      "Delta" -> Function[Evaluate[{ww, ss}], Evaluate[detB]], "DeltaDegree" -> deg,
      "DeltaCoefficients" -> Function[Evaluate[{ss}], Evaluate[cs]],
      "ValidityConditions" -> Function[Evaluate[{ss}], Evaluate[valid]], "ValidityExplanation" -> explanation,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowPropagator[stage];
   PrintCard[stage];
   stage]];

PropagatorMatrix[p_Association, w_, s_] := p["G"][w, s];
FullMatrix[p_Association, w_, s_] := p["M"][w, s];
AtomRules[p_Association] := Normal[p["AtomValues"]];
DeltaRules[p_Association] := {\[CapitalDelta][w_, s_] :> p["Delta"][w, s], DeltaCoefficient[m_][s_] :> p["DeltaCoefficients"][s][[m + 1]]};
ValidityConditions[p_Association] := p["ValidityConditions"];

ShowPropagator[p_Association] := Module[{s = Symbol["FRGscalarMSR`Propagator`Private`s"], w = Symbol["FRGscalarMSR`Propagator`Private`\[Omega]"]},
   Print["Propagator atoms (s = q^2):"];
   KeyValueMap[Print["  ", displayString[#1[s]], " = ", displayString[#2[s]]] &, p["AtomValues"]];
   Print["  Delta(", displayString[w], ", s) = det[(phi, phibar) block] = ", displayString[Collect[p["Delta"][w, s], w]]];
   Print["  Delta = c0 + c1 (I w) + c2 (I w)^2 with {c0, c1, c2} = ", displayString[p["DeltaCoefficients"][s]]];
   Print["  validity: ", displayString[p["ValidityConditions"][s]], "   -- ", p["ValidityExplanation"]]];

ShowPropagator[p_Association, {a_, b_}] := Module[{s = Symbol["FRGscalarMSR`Propagator`Private`s"], w = Symbol["FRGscalarMSR`Propagator`Private`\[Omega]"],
   i = First[FirstPosition[p["Fields"], a]], j = First[FirstPosition[p["Fields"], b]]},
   Print["G[", a, ", ", b, "](", displayString[w], ", s) = ", displayString[p["G"][w, s][[i, j]]]]];

End[];
EndPackage[];
