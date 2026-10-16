(* ::Package:: *)

(* Regulator.wl -- DefineRegulator: the regulator R_k as Gamma2-level entries per field pair, and its derivative d_s R
   (s = ln k) by the chain rule, with abstract shape functions r(y), y = q^2/k^2.
   Objective, assumptions and checks: docs/MODULES.md section "Regulator". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

DefineRegulator::usage = "DefineRegulator[ansatz, {a, b} -> R_ab, ..., \"Shapes\" -> {r1, r2, ...}] declares the regulator per field pair. \
R_ab is the Gamma2-level entry (what adds to Gamma2_ab, including the factor 2 of a diagonal pair), written with the loop momentum q and \
the RG scale k, shapes as r[q^2/k^2], e.g. {Hbar, H} -> nk q^2 arnh[q^2/k^2]. Returns the \"Regulator\" stage.";
dsCoupling::usage = "dsCoupling[c] = k dc/dk, the scale derivative of the coupling c (becomes -eta c downstream).";
dsR::usage = "dsR[a, b] is the abstract entry (a, b) of d_s R (fields in field order); it is replaced by the derived d_s R at assembly. \
Map to DE4plus_fast.wl: dsNh = dsR[H, Hbar], dsMh = dsR[Hbar, Hbar], dsNphi = dsR[Psi, Psibar], dsMphi = dsR[Psibar, Psibar].";
RegulatorMatrix::usage = "RegulatorMatrix[regulator, s] is the matrix R_ab at q^2 = s, in field order.";
DsRegulatorMatrix::usage = "DsRegulatorMatrix[regulator, s] is the derived matrix d_s R_ab at q^2 = s (dsCoupling[c] symbols, shape derivatives).";
DsRAbstractMatrix::usage = "DsRAbstractMatrix[regulator] is the matrix of abstract symbols dsR[a, b] (0 where R_ab = 0).";
DsRRules::usage = "DsRRules[regulator, s] are the rules dsR[a, b] -> derived d_s R_ab at q^2 = s.";
ShowRegulator::usage = "ShowRegulator[regulator] prints R_ab and d_s R_ab.";

Begin["`Private`"];

SetAttributes[{dsCoupling, dsR}, Protected];

symbolsNamed[e_, name_] := Union[Cases[e, x_Symbol /; SymbolName[x] === name, {0, Infinity}, Heads -> True]];

(* q^n -> sq^(n/2), n even: R depends on q only through q^2 *)
toS[e_, qsym_] := e /. Power[qsym, n_Integer?EvenQ] :> sq^(n/2);

Options[DefineRegulator] = {"Shapes" -> {}};

DefineRegulator[ansatz_, args___] := Block[{$FRGModule = "Regulator"},
  Module[{t0 = AbsoluteTime[], theory, fields, resp, cps, pairs, opts, shapes, qs, ks, kk, entries, checks = {}, warns = {}, allowed, bad,
    rmat, dsmat, dsAbs, pos, stage, used},
   validateStage[ansatz, "Ansatz"];
   theory = ansatz["Theory"]; fields = theory["Fields"]; resp = theory["ResponseFields"]; cps = ansatz["Couplings"];
   pairs = Cases[{args}, ({_, _} -> _) | ({_, _} :> _)];
   opts = Cases[{args}, (Rule | RuleDelayed)[_String, _]];
   FRGAssert["Regulator.arguments", Length[pairs] + Length[opts] == Length[{args}] && MatchQ[opts, {("Shapes" -> _) ...}],
    Complement[{args}, pairs, opts], "arguments are {a, b} -> R_ab rules and \"Shapes\" -> {r1, ...}"];
   shapes = Lookup[opts, "Shapes", {}];
   FRGAssert["Regulator.shapeList", ListQ[shapes] && shapes =!= {}, shapes, "\"Shapes\" -> {r1, r2, ...} is required: the shape functions stay abstract"];
   assertFreshSymbols["Regulator.shapes", shapes, "shape functions"];
   FRGAssert["Regulator.shapesDistinct", Intersection[shapes, Join[fields, cps, {theory["Dimension"]}]] === {}, Intersection[shapes, Join[fields, cps]],
    "a shape function has the name of a field, a coupling or the dimension"];
   FRGAssert["Regulator.pairs", pairs =!= {}, {args}, "give at least one pair {a, b} -> R_ab"];

   (* pairs: fields of the theory, at least one response field, each unordered pair once *)
   AppendTo[checks, assertCheck["Regulator.pairFields", AllTrue[pairs[[All, 1]], SubsetQ[fields, #] &], Select[pairs[[All, 1]], ! SubsetQ[fields, #] &],
     "regulator pairs must be fields of the theory", "pairs are fields of the theory"]];
   AppendTo[checks, assertCheck["Regulator.responsePair", AllTrue[pairs[[All, 1]], IntersectingQ[#, resp] &], Select[pairs[[All, 1]], ! IntersectingQ[#, resp] &],
     "a regulator entry between two physical fields would break causality (Gamma2_{phi phi} = 0): only phibar phi and phibar phibar pairs",
     "only phibar-phi / phibar-phibar pairs"]];
   pos[a_] := First[FirstPosition[fields, a]];
   pairs = (SortBy[#[[1]], pos] -> #[[2]]) & /@ pairs;
   AppendTo[checks, assertCheck["Regulator.pairOnce", DuplicateFreeQ[pairs[[All, 1]]], pairs[[All, 1]],
     "a pair is given twice ({a, b} and {b, a} are the same entry: R is symmetric because it is even in q and w-independent)", "each pair given once"]];

   (* the user writes q (loop momentum) and k (scale): any context, matched by name *)
   qs = symbolsNamed[pairs[[All, 2]], "q"]; ks = symbolsNamed[pairs[[All, 2]], "k"];
   FRGAssert["Regulator.momentumSymbol", Length[qs] <= 1 && Length[ks] == 1, {qs, ks},
    "write the entries with one loop-momentum symbol q and one scale symbol k, shapes as r[q^2/k^2]"];
   kk = First[ks];
   FRGAssert["Regulator.scaleFree", AllTrue[Join[ks, qs], OwnValues[#] === {} && DownValues[#] === {} &], Join[ks, qs], "k and q must have no values"];
   entries = Association[(#[[1]] -> If[qs === {}, #[[2]], toS[#[[2]], First[qs]]]) & /@ pairs];
   AppendTo[checks, assertCheck["Regulator.evenInQ", FreeQ[Values[entries], Alternatives @@ qs], Select[entries, ! FreeQ[#, Alternatives @@ qs] &],
     "R_ab must depend on q through q^2 only (even powers of q)", "R depends on q^2 only"]];

   allowed = Join[cps, {kk, sq}];
   bad = Complement[Cases[Values[entries], x_Symbol /; Context[x] =!= "System`", {0, Infinity}, Heads -> False], allowed];
   AppendTo[checks, assertCheck["Regulator.knownSymbols", bad === {}, bad,
     "R_ab may contain the couplings of the ansatz, q, k and the declared shapes only (in particular no frequency: R is w-independent)",
     "only couplings, q, k and shapes; w-independent"]];
   bad = Complement[Cases[Values[entries], (h_Symbol /; Context[h] =!= "System`")[___] :> h, {0, Infinity}], shapes];
   AppendTo[checks, assertCheck["Regulator.onlyDeclaredShapes", bad === {}, bad, "every function in R_ab must be a declared shape", "the declared shapes are the only free functions"]];
   bad = Cases[Values[entries], (r : Alternatives @@ shapes)[x___] /; {x} =!= {sq/kk^2}, {0, Infinity}];
   AppendTo[checks, assertCheck["Regulator.shapeArgument", bad === {} && FreeQ[Values[entries], Derivative[__][Alternatives @@ shapes]] && FreeQ[Values[entries] /. (Alternatives @@ shapes)[sq/kk^2] -> 1, Alternatives @@ shapes],
     bad, "shapes appear only as r[q^2/k^2] (no derivatives, no other argument)", "shapes appear as r[q^2/k^2]"]];
   used = Select[shapes, ! FreeQ[Values[entries], #] &];
   If[used =!= shapes, AppendTo[warns, infoCheck["WARN", "Regulator.unusedShapes", "declared but unused: " <> ToString[Complement[shapes, used]]]]];
   AppendTo[checks, assertCheck["Regulator.nonzero", AllTrue[Values[entries], # =!= 0 &], Select[entries, # === 0 &], "drop pairs with R_ab = 0", "no zero entry"]];

   (* d_s R by the chain rule: k d/dk at fixed q; couplings c(k) -> dsCoupling[c] *)
   rmat = Table[Lookup[entries, Key[SortBy[{a, b}, pos]], 0], {a, fields}, {b, fields}];
   dsmat = Map[Expand[Sum[D[#, c] dsCoupling[c], {c, cps}] + kk D[#, kk]] &, rmat, {2}];
   dsAbs = Table[If[rmat[[i, j]] === 0, 0, dsR @@ SortBy[{fields[[i]], fields[[j]]}, pos]], {i, Length[fields]}, {j, Length[fields]}];
   (* the same derivative by a second route: couplings as functions c(k), total k d/dk, c'(k) -> dsCoupling[c]/k *)
   AppendTo[checks, assertCheck["Regulator.dsRChainRule",
     AllTrue[Flatten[dsmat - Map[(kk D[# /. Thread[cps -> (cf[#][kk] & /@ cps)], kk] /. Derivative[1][cf[c_]][kk] :> dsCoupling[c]/kk /. cf[c_][kk] :> c) &, rmat, {2}]], ExactlyZeroQ],
     dsmat, "internal: the two computations of d_s R differ",
     "d_s R = sum_c dR/dc dsCoupling[c] + k dR/dk  (k d/dk r(q^2/k^2) = -2 y r'(y)); checked by a second route"]];

   stage = stampStage[<|"Stage" -> "Regulator", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {ansatz["Fingerprint"]},
      "Fields" -> fields, "Scale" -> kk, "Shapes" -> shapes, "Couplings" -> cps,
      "Entries" -> entries, "RMatrix" -> Function[Evaluate[{sq}], Evaluate[rmat]], "DsRMatrix" -> Function[Evaluate[{sq}], Evaluate[dsmat]],
      "DsRAbstract" -> dsAbs, "DsCouplings" -> Union[Cases[dsmat, _dsCoupling, Infinity]],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowRegulator[stage];
   PrintCard[stage];
   stage]];

RegulatorMatrix[r_Association, s_] := r["RMatrix"][s];
DsRegulatorMatrix[r_Association, s_] := r["DsRMatrix"][s];
DsRAbstractMatrix[r_Association] := r["DsRAbstract"];
DsRRules[r_Association, s_] := DeleteDuplicates[DeleteCases[Flatten[MapThread[Rule, {r["DsRAbstract"], r["DsRMatrix"][s]}, 2]], 0 -> _]];

ShowRegulator[r_Association] := Module[{qq = Symbol["FRGscalarMSR`Regulator`Private`q"]},
   Print["Regulator (Gamma2-level entries; R is symmetric, w-independent; shapes abstract, y = q^2/k^2):"];
   KeyValueMap[Print["  R[", #1[[1]], ", ", #1[[2]], "] = ", displayString[#2 /. sq -> qq^2]] &, r["Entries"]];
   Print["d_s R = k dR/dk at fixed q (dsCoupling[c] = k dc/dk):"];
   Scan[Print["  ", displayString[#[[1]]], " = ", displayString[#[[2]] /. sq -> qq^2]] &, DsRRules[r, sq]]];

End[];
EndPackage[];
