(* ::Package:: *)

(* Adimensional.wl -- AdimensionalEquations: from the dimensionful equations d_s(comb) = Integral d^dq/(2 Pi)^d (integrand) and the
   derived dimensions, the dimensionless beta functions of the state couplings, split into a linear part (no momentum integral) and a
   loop part (integrand of Integral daq), and the linear system for the anomalous dimensions eta_j of the reference couplings.
   Objective, assumptions and checks: docs/MODULES.md section "Adimensional". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

AdimensionalEquations::usage = "AdimensionalEquations[equations, dims, \"State\" -> {gk, lk, ...}] returns the \"Adimensional\" stage: \
for each state coupling c (dimensionless name chat) d_s chat = BetaL + Integral daq Measure[aq] BetaNL, and the eta system \
(EtaML + Integral daq Measure EtaMNL) eta = EtaBL + Integral daq Measure EtaBNL. d_s ref_j = -eta_j ref_j.";
ShowBeta::usage = "ShowBeta[adim, c] prints the dimensionless flow of the state coupling c: linear part, loop integrand, atoms used.";
ShowEtaSystem::usage = "ShowEtaSystem[adim] prints the eta system row by row (running references: eta from their equation; fixed references: \
placeholder row eta = 0, to be fixed at run time with FlowRG fixed_etas).";
AdimAtom::usage = "AdimAtom[a, b, m, j] is a dimensionless propagator atom (the scale-stripped atom ScaleAtom[a, b, m, j] in dimensionless \
couplings and aq); its value is adim[\"Atoms\"][AdimAtom[...]].";
ShowCausalPoles::usage = "ShowCausalPoles[adim] prints the causal-pole conditions in minimal dimensionless form: each entry must be > 0 at \
every aq (exported as causal_poles; design/CAUSALITY.md).";
LoopMeasure::usage = "LoopMeasure[d, aq] = 2 Pi^(d/2)/Gamma[d/2] aq^(d-1)/(2 Pi)^d: Integral d^dq/(2 Pi)^d f(q^2) = k^d Integral daq LoopMeasure[d, aq] f(k^2 aq^2).";

Begin["`Private`"];

SetAttributes[AdimAtom, Protected];

LoopMeasure[d_, aq_] := 2 Pi^(d/2)/Gamma[d/2] aq^(d - 1)/(2 Pi)^d;

(* ---------- homogeneity in the units (generalises HomVec/HomClean of AdimensionalFlowEQGenericRegLargeScale.nb) ----------
   hom[e] = {v, rest}: e = (product of units with exponent vector v) x rest, rest free of units. Every term of a sum must carry the
   same v; units may not appear in exponents or in function arguments. Exponents may depend on the symbolic dimension. *)
homVec[e_, unitIndex_Association, nb_] := Module[{zero = ConstantArray[0, nb], hv},
   hv[s_Symbol] := If[KeyExistsQ[unitIndex, s], {UnitVector[nb, unitIndex[s]], 1}, {zero, s}];
   hv[x_?NumericQ] := {zero, x};
   hv[x_Plus] := Module[{h = hv /@ List @@ x, v},
     v = h[[1, 1]];
     If[! AllTrue[h[[All, 1]], AllTrue[Expand[# - v], # === 0 &] &], Throw[{"a sum mixes dimensions", DeleteDuplicates[h[[All, 1]]], x}, homTag]];
     {v, Plus @@ h[[All, 2]]}];
   hv[x_Times] := Module[{h = hv /@ List @@ x}, {Expand[Total[h[[All, 1]]]], Times @@ h[[All, 2]]}];
   hv[Power[b_, n_]] := Module[{hb = hv[b], hn = hv[n]},
     If[hn[[1]] =!= zero, Throw[{"a unit in an exponent", Power[b, n]}, homTag]];
     {Expand[n hb[[1]]], hb[[2]]^n}];
   hv[f_[args___]] := Module[{h = hv /@ {args}},
     If[! AllTrue[h[[All, 1]], # === zero &], Throw[{"a unit inside a function argument", f[args]}, homTag]];
     {zero, f @@ h[[All, 2]]}];
   hv[e]];

(* dimension vector of e (units substituted); aborts with the reason if e is not homogeneous *)
dimensionVector[check_, e_, unitIndex_, nb_, what_] := Module[{r = Catch[homVec[e, unitIndex, nb], homTag]},
   FRGAssert[check, MatchQ[r, {_List, _}] && Length[r[[1]]] == nb, {what, r}, "not homogeneous in the units k, references: " <> ToString[First[r]]];
   r];
(* rest of e after asserting that its dimension is exactly v *)
homClean[check_, e_, v_, unitIndex_, nb_, what_] := Module[{r = dimensionVector[check, e, unitIndex, nb, what]},
   FRGAssert[check, AllTrue[Expand[r[[1]] - v], # === 0 &], {what, "residual dimension" -> Expand[r[[1]] - v]},
    "the expression still depends on the units k, references (dimension mismatch)"];
   r[[2]]];

(* independent route (review E8): e must be unchanged when every unit is multiplied by a random factor. Units are random perfect
   powers so that every fractional power stays rational; d is set to each test dimension. *)
unitInvariantQ[e_, units_, dimSym_, ds_, seed_] := AllTrue[ds, Function[d, Module[{ed = e /. dimSym -> d, den, pt, lam, memo, atom, val, syms},
     den = LCM @@ Denominator[Cases[ed, Power[_, x_?NumericQ] :> x, {0, Infinity}] /. {} -> {1}];
     SeedRandom[seed + d];
     syms = Complement[Union[Cases[ed, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> False]], units];
     pt = Thread[syms -> RandomInteger[{2, 40}, Length[syms]]/RandomInteger[{41, 97}, Length[syms]]];
     lam = Thread[units -> (RandomInteger[{2, 7}, Length[units]]/RandomInteger[{2, 7}, Length[units]])^den];
     memo[x_] := memo[x] = RandomInteger[{2, 60}]/RandomInteger[{61, 97}];
     (* not Core's whole-atom rule (design/MOMENTUMPROJECTION.md step 0): enough here, AdimAtom[..] takes no argument and shapes are symbols *)
     atom = x : (Derivative[__][_Symbol][___] | (h_Symbol /; Context[h] =!= "System`")[___]) :> memo[x];
     val[u_] := Together[(ed /. Thread[units -> (u units)] /. pt /. Thread[units -> (Last /@ lam)] /. atom)];
     Together[val[1] - val[Last /@ lam]] === 0]]];

(* ---------- causal-pole conditions (design/CAUSALITY.md) ---------- *)

(* numeric coefficient of one term of an expanded polynomial *)
causalPolesCoefficient[t_] := Which[NumberQ[t], t, Head[t] === Times, Times @@ Select[List @@ t, NumberQ], True, 1];
(* exponent of aq in one term, outside function arguments: rn[aq^2] does not count *)
causalPolesAqExponent[t_, aq_] := Exponent[t /. (h_Symbol /; Context[h] =!= "System`")[___] :> 1, aq];
(* {content, minimal}: rest = content minimal, content = positive numeric content x aq^(smallest exponent) *)
causalPolesStrip[rest_, aq_] := Module[{e = Expand[rest], terms, num, m},
   If[e === 0, Return[{1, 0}]];
   terms = If[Head[e] === Plus, List @@ e, {e}];
   num = Abs[GCD @@ (causalPolesCoefficient /@ terms)];
   m = Min[causalPolesAqExponent[#, aq] & /@ terms];
   {num aq^m, Expand[e/(num aq^m)]}];
(* a product of a positive number and powers of the given symbols *)
causalPolesPositiveQ[s_, syms_] := AllTrue[If[Head[s] === Times, List @@ s, {s}],
   (NumberQ[#] && # > 0) || MatchQ[#, Alternatives @@ syms] || MatchQ[#, Power[Alternatives @@ syms, _]] &];

(* ---------- AdimensionalEquations ---------- *)

Options[AdimensionalEquations] = {"State" -> Automatic};

AdimensionalEquations[equations_, dims_, opts : OptionsPattern[]] := Block[{$FRGModule = "Adimensional"},
  Module[{t0 = AbsoluteTime[], checks = {}, warns = {}, ansatz, theory, fields, cps, refs, etas, names, vecs, basis, nb, unitIndex, kk, qU, aq,
    dimSym, nonren, regC, reg, eqs, combs, toUnits, form, gamma2Dim, kind, defEq, running, fixed, state, dsC, lhs, unknowns, uEqs, mat, rhsVec,
    sol, loopSym, dsOf, betaL, betaNL, betaNLExplicit, atomVals, atomVecs, adimAtoms, compactPiece, explicitPiece, etaRows, mL, mNL, bL, bNL,
    mustFix, stage, testDims, cand, bad, compactAtoms, etaRowsExplicit, mNLx, bNLx, causalRaw, causalOrig, causalStripped, causalMin,
    causalExpl, causalKeep},
   validateStage[equations, "Equations"]; validateStage[dims, "Dimensions"];
   FRGAssert["Adimensional.options", MatchQ[{opts}, {(Rule | RuleDelayed)["State", _] ...}], {opts}, "only the option \"State\" is accepted"];
   ansatz = AnsatzOf[equations];
   AppendTo[checks, assertCheck["Adimensional.sameAnsatz", ansatz["Fingerprint"] === dims["Ansatz"]["Fingerprint"],
     {ansatz["Fingerprint"], dims["Ansatz"]["Fingerprint"]}, "the dimensions were derived from another ansatz: DeriveDimensions[AnsatzOf[equations], ...]",
     "dimensions and equations from the same ansatz"]];
   theory = ansatz["Theory"]; fields = theory["Fields"]; cps = ansatz["Couplings"];
   refs = dims["References"]; etas = dims["Etas"]; names = dims["Names"]; vecs = dims["Vectors"]; basis = dims["Basis"]; nb = Length[basis];
   kk = dims["Scale"]; dimSym = dims["Dimension"]; qU = equations["LoopMomentum"]; aq = Symbol["Global`aq"];
   FRGAssert["Adimensional.scale", equations["Scale"] === kk, {equations["Scale"], kk}, "the regulator scale and the \"Scale\" of DeriveDimensions differ"];
   unitIndex = AssociationThread[basis -> Range[nb]];
   nonren = First /@ EmbeddedStage[equations, "Symmetries"]["NonRenormalised"];
   testDims = If[dimSym === 2, {2}, {2, 3}];

   AppendTo[checks, assertCheck["Adimensional.referencesFlow", Intersection[refs, nonren] === {}, Intersection[refs, nonren],
     "a NonRenormalised coupling cannot be a reference: choose another one", "no reference is NonRenormalised"]];

   (* D4: d_s R may contain d_s of reference couplings only *)
   regC = Union[Cases[equations["DsRRules"][[All, 2]], dsCoupling[c_] :> c, {0, Infinity}]];
   AppendTo[checks, assertCheck["Adimensional.regulatorReferences", SubsetQ[refs, regC], Complement[regC, refs],
     "regulator prefactors may contain reference couplings only (e.g. nk q^2 arnh, not nk (1 + Ak) q^2 arnh); lifting this is in docs/FUTURE.md",
     "regulator prefactors contain references only: " <> StringRiffle[ToString /@ regC, ", "]]];
   (* units: c -> chat [c] for non-references, q -> k aq, d_s ref -> -eta ref *)
   form[c_] := basisForm[vecs[c], basis];
   toUnits = Join[(# -> names[#] form[#]) & /@ Complement[cps, refs], {qU -> kk aq}, (dsCoupling[#] -> -etas[#] #) & /@ refs];
   gamma2Dim[{a_, b_}] := Expand[dimSym UnitVector[nb, 1] - dims["Time"] - vecs[a] - vecs[b]];

   (* check 3: each regulator entry has the dimension of the Gamma2 entry it adds to *)
   reg = EmbeddedStage[equations, "Regulator"];
   bad = Select[Normal[reg["Entries"]], ! AllTrue[Expand[dimensionVector["Adimensional.regulatorDimension",
          (#[[2]] /. Symbol["FRGscalarMSR`Private`sq"] -> qU^2) /. toUnits, unitIndex, nb, #[[1]]][[1]] - gamma2Dim[#[[1]]]], # === 0 &] &];
   AppendTo[checks, assertCheck["Adimensional.regulatorDimension", bad === {}, bad,
     "a regulator entry R_ab does not have the dimension of Gamma2_ab (Dim e_k - [t] - [a] - [b]): its prefactor is wrong",
     "every R_ab has the dimension of Gamma2_ab"]];

   (* check 5: the Gamma_k coefficient of each equation has the dimension of its monomial (iw)^m (p^2)^n in Gamma2_ab *)
   eqs = equations["Equations"]; combs = eqs[[All, "Flowing"]];
   bad = Select[eqs, ! AllTrue[Expand[dimensionVector["Adimensional.lhsDimension", #GammaK /. toUnits, unitIndex, nb, #Flowing][[1]] -
          (gamma2Dim[#Pair] + #Frequency dims["Time"] - 2 #Momentum UnitVector[nb, 1])], # === 0 &] &];
   AppendTo[checks, assertCheck["Adimensional.lhsDimension", bad === {}, bad[[All, "Flowing"]],
     "the coefficient of a declared monomial does not have the dimension of its Gamma2 entry", "each Gamma_k coefficient has the dimension of its monomial"]];

   (* references: running (eta from the equation of ref or ref x amplitude) or fixed (eta given at run time) *)
   cand[r_] := Select[Range[Length[eqs]], With[{q = Together[combs[[#]]/r]},
        PolynomialQ[q, cps] && FreeQ[q /. Thread[dims["Dimensionless"] -> 1], Alternatives @@ cps] && Exponent[combs[[#]], r] == 1] &];
   kind = <||>; defEq = <||>;
   Do[Module[{c = cand[r], alone, mode = dims["ReferenceModes"][r], eqSpec = dims["ReferenceEquations"][r]},
      alone = Select[c, Expand[combs[[#]] - r] === 0 &];
      Which[
       alone =!= {},
       FRGAssert["Adimensional.referenceMode", mode =!= "Fixed", r, "this reference appears alone in an equation, so its eta is determined (Running); drop \"Mode\" -> \"Fixed\""];
       kind[r] = "Running"; defEq[r] = First[alone],
       c === {},
       FRGAssert["Adimensional.referenceEquation", False, r, "no equation has this reference alone or times dimensionless amplitudes: it cannot set an eta"],
       mode === "Running",
       If[eqSpec =!= Automatic, c = Select[c, Expand[combs[[#]] - eqSpec] === 0 &]];
       FRGAssert["Adimensional.referenceEquation", Length[c] == 1, {r, combs[[cand[r]]]},
        "several (or no) equations can define this eta: choose one with \"Equation\" -> comb in DeriveDimensions \"References\""];
       kind[r] = "Running"; defEq[r] = First[c],
       True, kind[r] = "Fixed"]], {r, refs}];
   running = Select[refs, kind[#] === "Running" &]; fixed = Select[refs, kind[#] === "Fixed" &];
   AppendTo[warns, infoCheck["INFO", "Adimensional.references", StringRiffle[(ToString[#] <> ": " <> kind[#] <>
        If[kind[#] === "Running", " (eta from d_s(" <> inputForm[combs[[defEq[#]]]] <> "))", " (eta" <> dims["EtaNames"][#] <> " fixed at run time)"]) & /@ refs, ", "]]];

   (* Leibniz: d_s(comb) = sum_c (d comb/d c) d_s c, d_s ref = -eta ref, d_s c = 0 for NonRenormalised c; one loop symbol per equation *)
   loopSym = Array[loopOf, Length[eqs]];
   lhs = Function[comb, Expand[Total[(D[comb, #] dsC[#]) & /@ Variables[comb]] /. Join[(dsC[#] -> -etas[#] #) & /@ refs, (dsC[#] -> 0) & /@ nonren]]] /@ combs;
   unknowns = Union[Cases[lhs, _dsC, {0, Infinity}]];
   uEqs = Select[Range[Length[eqs]], ! FreeQ[lhs[[#]], dsC] &];
   bad = Complement[Range[Length[eqs]], uEqs, Values[defEq]];
   FRGAssert["Adimensional.unusedEquation", bad === {}, combs[[bad]], "these equations determine no coupling and no eta (check the references and NonRenormalised tags)"];
   FRGAssert["Adimensional.leibnizSquare", Length[uEqs] == Length[unknowns], {combs[[uEqs]], unknowns},
    "the equations do not determine each flowing coupling exactly once (Leibniz system not square)"];
   sol = If[unknowns === {}, <||>,
     mat = Table[Coefficient[lhs[[e]], u], {e, uEqs}, {u, unknowns}];
     FRGAssert["Adimensional.leibnizSolvable", ! ExactlyZeroQ[Det[mat]], {combs[[uEqs]], unknowns}, "the Leibniz system is singular"];
     rhsVec = loopSym[[uEqs]] - (lhs[[uEqs]] /. _dsC -> 0);
     AssociationThread[(First /@ unknowns) -> Expand[LinearSolve[mat, rhsVec]]]];
   AppendTo[checks, assertCheck["Adimensional.leibniz", True, Null, "", ToString[Length[unknowns]] <> " flowing couplings solved from " <>
       ToString[Length[uEqs]] <> " equations: " <> StringRiffle[ToString /@ Keys[sol], ", "]]];

   (* state *)
   state = Replace[OptionValue["State"], Automatic -> Select[cps, ! MemberQ[refs, #] &]];
   FRGAssert["Adimensional.state", ListQ[state] && DuplicateFreeQ[state] && SubsetQ[Complement[cps, refs], state], state,
    "\"State\" -> {couplings of the ansatz that are not references}"];
   bad = Select[state, ! KeyExistsQ[sol, #] && ! MemberQ[nonren, #] &];
   FRGAssert["Adimensional.noBeta", bad === {}, bad, "these state couplings get no flow: no equation determines them and they are not NonRenormalised"];
   bad = Complement[Keys[sol], state];
   FRGAssert["Adimensional.stateMissing", bad === {}, bad, "these couplings have a flow equation but are missing from \"State\""];
   dsOf[c_] := Lookup[sol, c, 0];

   (* dimensionless atoms from the scale-stripped atoms of Assemble *)
   atomVals = equations["Atoms"];
   atomVecs = Association[KeyValueMap[#1 -> dimensionVector["Adimensional.atomDimension", #2 /. toUnits, unitIndex, nb, #1] &, atomVals]];
   adimAtoms = Association[KeyValueMap[(AdimAtom @@ #1) -> #2[[2]] &, atomVecs]];
   compactAtoms = KeyValueMap[#1 -> basisForm[#2[[1]], basis] (AdimAtom @@ #1) &, atomVecs];

   (* causal-pole conditions (design/CAUSALITY.md): clear denominators, units out (references > 0), positive content stripped *)
   causalRaw = equations["ValidityConditions"]["Atoms"] /. True -> Nothing;
   FRGAssert["Adimensional.causalPoles", FreeQ[causalRaw, False] && MatchQ[causalRaw, {Greater[_, 0] ...}], causalRaw,
    "the causal-pole conditions must be a list of e > 0; False means Delta has an acausal root at every state"];
   {causalOrig, causalStripped, causalMin} = If[causalRaw === {}, {{}, {}, {}}, Transpose[Function[cond, Module[{nd, p, x, vr, cont, minimal},
         nd = NumeratorDenominator[Together[First[cond]]];
         p = If[NumberQ[nd[[2]]], Sign[nd[[2]]] nd[[1]], nd[[1]] nd[[2]]];
         x = (p /. compactAtoms) /. toUnits;
         vr = dimensionVector["Adimensional.causalPolesDimension", x, unitIndex, nb, cond];
         {cont, minimal} = causalPolesStrip[vr[[2]], aq];
         {x, basisForm[vr[[1]], basis] cont, minimal}]] /@ causalRaw]];
   bad = Pick[causalRaw, MapThread[! (causalPolesPositiveQ[#2, Append[basis, aq]] && ExactlyZeroQ[#1 - #2 #3]) &, {causalOrig, causalStripped, causalMin}]];
   AppendTo[checks, assertCheck["Adimensional.causalPolesMinimal", bad === {}, bad,
     "internal: a minimal causal-pole condition is not the original one divided by a positive factor (number, references, k, aq)",
     "minimal causal-pole conditions = original / (positive number x powers of k, references, aq), exact"]];
   causalExpl = causalMin /. Normal[adimAtoms];
   bad = Pick[causalRaw, (NumberQ[#] && # <= 0) & /@ (Together /@ causalExpl)];
   FRGAssert["Adimensional.causalPoles", bad === {}, bad,
    "a causal-pole condition reduces to a constant <= 0: Delta has an acausal root at every state (check the ansatz signs)"];
   causalKeep = ! NumberQ[Together[#]] & /@ causalExpl;
   {causalStripped, causalMin, causalExpl} = Pick[#, causalKeep] & /@ {causalStripped, causalMin, causalExpl};
   AppendTo[warns, infoCheck["ASSUMPTION", "Adimensional.referencesPositive", "references " <> StringRiffle[ToString /@ refs, ", "] <>
      ", k and aq are > 0 (their powers are stripped from the causal-pole conditions)"]];
   AppendTo[warns, infoCheck["ASSUMPTION", "Adimensional.causalPoles", If[causalExpl === {}, "causal poles: trivially satisfied",
      StringRiffle[(displayString[#] <> " > 0") & /@ causalExpl, ", "] <> " at every aq (causal poles; exported as causal_poles, not checked here)"]]];

   (* a loop piece  alpha x Integral d^dq/(2 Pi)^d I(q) / [c]  as the integrand of Integral daq LoopMeasure: alpha k^d I(k aq)/[c] *)
   explicitPiece[alpha_, e_, v_, what_] := Module[{x = (alpha kk^dimSym eqs[[e]]["Integrand"]) /. toUnits},
     FRGAssert["Adimensional.unitInvariance", unitInvariantQ[x, basis, dimSym, testDims, e], what,
      "the loop integrand changes when the units are rescaled: it is not dimensionless (independent route)"];
     homClean["Adimensional.homogeneous", x, v, unitIndex, nb, what]];
   compactPiece[alpha_, e_, v_, what_] := Module[{x = (alpha kk^dimSym (eqs[[e]]["Compact"] /. equations["DsRRules"] /. compactAtoms)) /. toUnits},
     homClean["Adimensional.homogeneousCompact", x, v, unitIndex, nb, what]];

   (* beta functions: d_s chat = chat (-x + Sum y eta) + (d_s c)/[c], local part not integrated, loop part per Integral daq *)
   betaL = <||>; betaNL = <||>; betaNLExplicit = <||>;
   Do[Module[{ds = dsOf[c], v = vecs[c], hat = names[c], local, alphas},
      alphas = Coefficient[ds, #] & /@ loopSym;
      local = Expand[ds /. Thread[loopSym -> 0]];
      betaL[c] = Expand[hat linearPart[v, Values[etas]] + homClean["Adimensional.homogeneousLocal", (local/form[c]) /. toUnits, ConstantArray[0, nb], unitIndex, nb, c]];
      betaNL[c] = Total[MapThread[If[#1 === 0 || eqs[[#2]]["IdenticallyZero"], 0, compactPiece[#1/form[c], #2, ConstantArray[0, nb], c]] &, {alphas, Range[Length[eqs]]}]];
      betaNLExplicit[c] = Total[MapThread[If[#1 === 0 || eqs[[#2]]["IdenticallyZero"], 0, explicitPiece[#1/form[c], #2, ConstantArray[0, nb], c]] &, {alphas, Range[Length[eqs]]}]]],
    {c, state}];
   AppendTo[checks, assertCheck["Adimensional.homogeneous", True, Null, "",
     "every beta free of k and the references: homogeneity of every sum (structural) and invariance under random unit rescaling (numerical, d = " <>
      StringRiffle[ToString /@ testDims, ", "] <> ")"]];
   AppendTo[checks, assertCheck["Adimensional.linearNoMomentum", FreeQ[Values[betaL], aq | qU | kk], Select[betaL, ! FreeQ[#, aq | qU | kk] &],
     "internal: a linear part depends on the loop momentum", "linear parts free of aq (not integrated)"]];

   (* eta system: running ref r: eta_r = -(loop of its equation)/r (amplitude counted as 1); fixed ref: placeholder row eta = 0 *)
   etaRows = Table[If[kind[r] === "Running", compactPiece[-1/r, defEq[r], ConstantArray[0, nb], "eta" <> dims["EtaNames"][r]], None], {r, refs}];
   etaRowsExplicit = Table[If[kind[r] === "Running", explicitPiece[-1/r, defEq[r], ConstantArray[0, nb], "eta" <> dims["EtaNames"][r]], 0], {r, refs}];
   {bNL, mNL} = Module[{rows = Replace[etaRows, None -> 0, {1}], ca},
     ca = Normal[CoefficientArrays[rows, Values[etas]]];
     FRGAssert["Adimensional.etaLinear", Length[ca] == 2, rows, "the eta system is not linear in the etas"];
     {ca[[1]], -ca[[2]]}];
   {bNLx, mNLx} = Module[{ca = Normal[CoefficientArrays[etaRowsExplicit, Values[etas]]]}, {ca[[1]], -ca[[2]]}];
   mL = IdentityMatrix[Length[refs]]; bL = ConstantArray[0, Length[refs]];
   mustFix = dims["EtaNames"] /@ fixed;
   AppendTo[checks, assertCheck["Adimensional.etaSystem", True, Null, "",
     "eta system: " <> ToString[Length[running]] <> " running row(s)" <> If[fixed === {}, "", ", placeholder row(s) eta = 0 for " <>
       StringRiffle[mustFix, ", "] <> " (MUST_FIX: give them in FlowRG fixed_etas)"]]];

   (* compact = explicit, at random points (atoms substituted) *)
   bad = Join[Select[state, ! ExactlyZeroQ[(betaNL[#] /. Normal[adimAtoms]) - betaNLExplicit[#], "Seeds" -> {1, 2}] &],
     Pick[refs, MapThread[! ExactlyZeroQ[(Replace[#1, None -> 0] /. Normal[adimAtoms]) - #2, "Seeds" -> {1, 2}] &, {etaRows, etaRowsExplicit}]]];
   AppendTo[checks, assertCheck["Adimensional.compactExact", bad === {}, bad, "internal: the compact loop integrand differs from the explicit one",
     "compact (dimensionless atoms) = explicit loop integrands (betas and eta rows), exact at random points"]];

   stage = stampStage[<|"Stage" -> "Adimensional", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {equations["Fingerprint"], dims["Fingerprint"]},
      "Equations" -> equations, "Dimensions" -> dims, "State" -> state, "Couplings" -> (names /@ state),
      "References" -> refs, "ReferenceKinds" -> kind, "Etas" -> Values[etas], "EtaNames" -> (dims["EtaNames"] /@ refs), "MustFix" -> mustFix,
      "BetaL" -> (betaL /@ state), "BetaNL" -> (betaNL /@ state), "BetaNLExplicit" -> (betaNLExplicit /@ state),
      "EtaML" -> mL, "EtaMNL" -> mNL, "EtaBL" -> bL, "EtaBNL" -> bNL, "EtaMNLExplicit" -> mNLx, "EtaBNLExplicit" -> bNLx, "Atoms" -> adimAtoms,
      "CausalPoles" -> <|"Compact" -> causalMin, "Explicit" -> causalExpl, "Stripped" -> causalStripped|>,
      "LoopMomentum" -> aq, "Measure" -> LoopMeasure[dimSym, aq], "Dimension" -> dimSym, "Shapes" -> equations["Shapes"],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

(* ---------- inspection ---------- *)

ShowBeta[adim_Association, c_] := Block[{$FRGModule = "Adimensional"}, Module[{i, used},
   validateStage[adim, "Adimensional"];
   FRGAssert["Adimensional.notInState", MemberQ[adim["State"], c], c, "not a state coupling; state: " <> ToString[adim["State"]]];
   i = First[FirstPosition[adim["State"], c]];
   Print["d_s ", adim["Couplings"][[i]], " = ", displayString[adim["BetaL"][[i]]], "   (linear part, not integrated)"];
   Print["    + Integral daq ", displayString[adim["Measure"]], " x"];
   Print["      ", displayString[adim["BetaNL"][[i]]]];
   used = Union[Cases[adim["BetaNL"][[i]], _AdimAtom, {0, Infinity}]];
   If[used =!= {}, Print["atoms:"]; Scan[Print["  ", displayString[#], " = ", displayString[adim["Atoms"][#]]] &, used]]]];

ShowEtaSystem[adim_Association] := Block[{$FRGModule = "Adimensional"}, Module[{refs, n},
   validateStage[adim, "Adimensional"];
   refs = adim["References"]; n = Length[refs];
   Print["(EtaML + Integral daq Measure EtaMNL) eta = EtaBL + Integral daq Measure EtaBNL,  eta = ", adim["Etas"]];
   Do[Print["row ", adim["Etas"][[j]], " (", adim["ReferenceKinds"][refs[[j]]], If[adim["ReferenceKinds"][refs[[j]]] === "Fixed", ": placeholder eta = 0, MUST_FIX", ""], "):"];
    Print["   ML  = ", adim["EtaML"][[j]], "   BL = ", adim["EtaBL"][[j]]];
    Print["   MNL = ", inputForm[adim["EtaMNL"][[j]], 6]];
    Print["   BNL = ", inputForm[adim["EtaBNL"][[j]], 6]], {j, n}]]];

ShowCausalPoles[adim_Association] := Block[{$FRGModule = "Adimensional"}, Module[{cp, used},
   validateStage[adim, "Adimensional"];
   cp = adim["CausalPoles"];
   If[cp["Compact"] === {}, Print["causal poles: trivially satisfied (no condition left)"]; Return[Null]];
   Print["causal poles: each entry > 0 at every aq (references, k, aq > 0 assumed; exported as causal_poles):"];
   MapThread[Print["  ", displayString[#1], " > 0", If[#1 === #2, "", "   (compact: " <> displayString[#2] <> ")"]] &, {cp["Explicit"], cp["Compact"]}];
   used = Union[Cases[cp["Compact"], _AdimAtom, {0, Infinity}]];
   If[used =!= {}, Print["atoms:"]; Scan[Print["  ", displayString[#], " = ", displayString[adim["Atoms"][#]]] &, used]]]];

End[];
EndPackage[];
