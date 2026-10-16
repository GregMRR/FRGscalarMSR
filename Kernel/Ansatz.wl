(* ::Package:: *)

(* Ansatz.wl -- DefineAnsatz: parse an MSR action density written with the package operators dt, grad, div, lap into a
   table of monomials  coefficient x product of (derivatives of single fields).
   Objective, assumptions and checks: docs/MODULES.md section "Ansatz". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

dt::usage = "dt[e]: time derivative in the ansatz density (package-owned and inert; not System`Dt).";
grad::usage = "grad[e]: spatial gradient (a vector) in the ansatz density; contract two of them with a dot: grad[a].grad[b].";
div::usage = "div[v]: divergence of a vector expression in the ansatz density, e.g. div[H grad[Psi]].";
lap::usage = "lap[e]: Laplacian in the ansatz density.";
DefineAnsatz::usage = "DefineAnsatz[theory, \"Density\" -> density, \"Couplings\" -> {c1, c2, ...}] parses the density \
(Gamma_k = Integral dt d^dx density) into monomials, prints the parsed table and returns the \"Ansatz\" stage.";
FieldFactor::usage = "FieldFactor[phi, nt, {i, ...}] is the factor dt^nt d_i ... phi of a parsed monomial; a repeated index \
SpaceIndex[1] is summed over (one contraction at most, since momentum order <= 2).";
SpaceIndex::usage = "SpaceIndex[n] is a summed spatial index in a parsed monomial.";
ShowMonomials::usage = "ShowMonomials[ansatz] prints the table of parsed monomials (coefficient | fields | operators).";
MonomialString::usage = "MonomialString[monomial] writes the field part of a parsed monomial back in the DSL.";

Begin["`Private`"];

SetAttributes[{dt, grad, div, lap, FieldFactor, SpaceIndex}, Protected];

opHeads = dt | grad | div | lap;
systemOps = Dt | D | Grad | Div | Laplacian | Curl | Derivative | Integrate | Inactive | Inactivate;

(* ---------- vector rank of a raw density expression (0 scalar, 1 vector); fails on ill-formed input ---------- *)

rank[e_, fld_] := Which[
   FreeQ[e, fld | opHeads], 0,
   MatchQ[e, fld], 0,
   True, rankOf[e, fld]];
rankOf[grad[a_], fld_] := (needRank[a, 0, fld, "grad[...] acts on a scalar"]; 1);
rankOf[dt[a_], fld_] := rank[a, fld];
rankOf[div[a_], fld_] := (needRank[a, 1, fld, "div[...] acts on a vector, e.g. div[H grad[Psi]]"]; 0);
rankOf[lap[a_], fld_] := (needRank[a, 0, fld, "lap[...] acts on a scalar"]; 0);
rankOf[Dot[a_, b_], fld_] := (needRank[a, 1, fld, "a dot contracts two vectors"]; needRank[b, 1, fld, "a dot contracts two vectors"]; 0);
rankOf[e_Plus, fld_] := Module[{r = rank[#, fld] & /@ List @@ e},
   FRGAssert["Ansatz.vectorSum", Equal @@ r, e, "a sum mixes scalar and vector terms"]; First[r]];
rankOf[e_Times, fld_] := Module[{r = Total[rank[#, fld] & /@ List @@ e]},
   FRGAssert["Ansatz.vectorProduct", r <= 1, e, "a product of two vectors must be written as a dot: grad[a].grad[b]"]; r];
rankOf[Power[a_, n_], fld_] := (FRGAssert["Ansatz.polynomial", IntegerQ[n] && n > 0, Power[a, n],
    "the density must be polynomial in the fields (positive integer powers only)"]; needRank[a, 0, fld, "only scalars can be raised to a power"]; 0);
rankOf[e_, fld_] := FRGAssert["Ansatz.unknownStructure", False, e,
   "unknown structure: use products, sums, positive powers and the operators dt, grad, div, lap, and . between two grads"];
needRank[a_, r_, fld_, hint_] := FRGAssert["Ansatz.vectorRank", rank[a, fld] === r, a, hint];

(* ---------- conversion to products of FieldFactor[phi, nt, {indices}] ---------- *)

(* derivations acting on a sum of products of FieldFactors *)
derive[op_][e_Plus] := derive[op] /@ e;
derive[op_][e_Times] := Sum[MapAt[derive[op], e, i], {i, Length[e]}];
derive[op_][Power[f_, n_Integer]] := n f^(n - 1) derive[op][f];
derive[op_][FieldFactor[f_, nt_, ix_]] := Replace[op, {"t" -> FieldFactor[f, nt + 1, ix], i_ :> FieldFactor[f, nt, Append[ix, i]]}];
derive[op_][c_] /; FreeQ[c, FieldFactor] := 0;

toFactors[e_, fld_, new_] := Module[{S, V},
   S[x_] /; FreeQ[x, fld | opHeads] := x;
   S[x : fld] := FieldFactor[x, 0, {}];
   S[x_Plus] := S /@ x;
   S[x_Times] := S /@ x;
   S[Power[x_, n_]] := Expand[S[x]^n];
   S[dt[x_]] := Expand[derive["t"][S[x]]];
   S[lap[x_]] := With[{j = new[]}, Expand[derive[j][derive[j][S[x]]]]];
   S[div[v_]] := With[{j = new[]}, Expand[derive[j][V[v, j]]]];
   S[Dot[u_, v_]] := With[{j = new[]}, Expand[V[u, j] V[v, j]]];
   V[grad[x_], i_] := Expand[derive[i][S[x]]];
   V[dt[v_], i_] := Expand[derive["t"][V[v, i]]];
   V[x_Plus, i_] := V[#, i] & /@ x;
   V[x_Times, i_] := Module[{fs = List @@ x, r},
     r = rank[#, fld] & /@ fs;
     Expand[S[Times @@ Pick[fs, r, 0]] V[First[Pick[fs, r, 1]], i]]];
   V[x_, i_] := FRGAssert["Ansatz.vectorRank", False, x, "a vector was expected here (grad[...] or a sum/product containing one)"];
   S[x_] := FRGAssert["Ansatz.vectorRank", False, x, "a scalar was expected here; a vector must be contracted (grad[a].grad[b]) or divided (div[...])"];
   S[e]];

(* ---------- monomials ---------- *)

(* rename the (at most one) summed index to SpaceIndex[1] so that equal monomials combine *)
canonical[factors_List] := Module[{ix = Union[Flatten[factors[[All, 3]]]]},
   Sort[factors /. Thread[ix -> Table[SpaceIndex[n], {n, Length[ix]}]]]];

splitTerm[t_] := Module[{fs = If[Head[t] === Times, List @@ t, {t}]},
   FRGAssert["Ansatz.polynomial", AllTrue[fs, MatchQ[#, _FieldFactor | Power[_FieldFactor, _Integer?Positive] | _?(FreeQ[#, FieldFactor] &)] &], t,
    "the density must be polynomial in the fields"];
   {Times @@ Select[fs, FreeQ[#, FieldFactor] &],
    Flatten[Replace[Select[fs, ! FreeQ[#, FieldFactor] &], {Power[f_, n_] :> Table[f, n], f_ :> {f}}, {1}]]}];

fieldString[FieldFactor[f_, nt_, ix_]] := Nest["dt[" <> # <> "]" &,
   Which[ix === {}, ToString[f], Length[ix] == 2 && ix[[1]] === ix[[2]], "lap[" <> ToString[f] <> "]",
    Length[ix] == 1, "grad[" <> ToString[f] <> "]", True, "d" <> ToString[ix] <> "[" <> ToString[f] <> "]"], nt];
(* factors carrying one index are the two ends of the contraction: printed as grad[a].grad[b] *)
MonomialString[m_Association] := Module[{fs = m["Factors"], pair},
   pair = Select[fs, Length[#[[3]]] == 1 &];
   StringRiffle[Join[fieldString /@ Select[fs, Length[#[[3]]] =!= 1 &], If[pair === {}, {}, {StringRiffle[fieldString /@ pair, "."]}]], " "]];

(* factors stored response fields first, then in field order (display order; identity was fixed by canonical) *)
monomialRecord[{coef_, fac0_}, fieldOrder_, responses_] := Module[{nIx = Length[Flatten[fac0[[All, 3]]]], fac},
   fac = SortBy[fac0, {! MemberQ[responses, #[[1]]], First[FirstPosition[fieldOrder, #[[1]]]], #[[2]], Length[#[[3]]]} &];
   <|"Coefficient" -> coef, "Factors" -> fac, "Fields" -> fac[[All, 1]],
    "FieldCount" -> Length[fac], "ResponseCount" -> Count[fac[[All, 1]], Alternatives @@ responses],
    "TimeOrder" -> Total[fac[[All, 2]]], "MomentumOrder" -> nIx|>];

(* ---------- DefineAnsatz ---------- *)

Options[DefineAnsatz] = {"Density" -> None, "Couplings" -> {}, "Pipeline" -> "DerivativeExpansion"};

DefineAnsatz[theory_, opts : OptionsPattern[]] := Block[{$FRGModule = "Ansatz"},
  Module[{t0 = AbsoluteTime[], dens, cps, fields, fld, responses, checks = {}, warns = {}, counter = 0, new, ops, expanded, terms, mons,
    badSyms, unused, stage},
   validateStage[theory, "Theory"];
   FRGAssert["Ansatz.options", MatchQ[{opts}, {(Rule | RuleDelayed)["Density" | "Couplings" | "Pipeline", _] ...}], {opts},
    "only the options \"Density\", \"Couplings\", \"Pipeline\" are accepted"];
   dens = OptionValue["Density"]; cps = OptionValue["Couplings"];
   fields = theory["Fields"]; responses = theory["ResponseFields"]; fld = Alternatives @@ fields;
   AppendTo[checks, assertCheck["Ansatz.pipeline", OptionValue["Pipeline"] === "DerivativeExpansion", OptionValue["Pipeline"],
     "v1 implements only \"DerivativeExpansion\" (couplings depend on k only); see docs/FUTURE.md", "DerivativeExpansion"]];
   FRGAssert["Ansatz.density", dens =!= None && dens =!= 0, dens, "\"Density\" -> expression is required"];
   FRGAssert["Ansatz.couplingList", ListQ[cps], cps, "\"Couplings\" -> {c1, c2, ...}"];
   assertFreshSymbols["Ansatz.couplings", cps, "couplings"];
   FRGAssert["Ansatz.couplingsNotFields", Intersection[cps, Append[fields, theory["Dimension"]]] === {}, Intersection[cps, Append[fields, theory["Dimension"]]],
    "a coupling has the name of a field or of the dimension"];

   (* E7: System derivative heads would silently differentiate the couplings (Dt[nk H] = H Dt[nk] + nk Dt[H]) *)
   AppendTo[checks, assertCheck["Ansatz.noSystemDerivatives", FreeQ[dens, systemOps], Cases[dens, systemOps[___] | Derivative[__][_], {0, Infinity}, Heads -> True],
     "use the package operators dt, grad, div, lap (System Dt, Grad, Div, Laplacian, D are not allowed)", "no Dt/D/Grad/Div/Laplacian"]];

   badSyms = Complement[Cases[dens, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> True], fields, cps, {theory["Dimension"], dt, grad, div, lap}];
   AppendTo[checks, assertCheck["Ansatz.knownSymbols", badSyms === {}, badSyms,
     If[MemberQ[SymbolName /@ badSyms, "t" | "x"], "couplings depend on k only: no explicit t or x dependence is allowed",
      "every symbol must be a declared coupling (\"Couplings\") or a field of the theory"], "every symbol is a field, a coupling or Dim"]];

   (* operators act on field expressions only: couplings depend on k only and must stand outside (no dt[nk H], no lap[nk]) *)
   ops = Cases[dens, opHeads[___], {0, Infinity}];
   FRGAssert["Ansatz.operatorArity", AllTrue[ops, Length[#] == 1 &], Select[ops, Length[#] =!= 1 &], "dt, grad, div, lap take exactly one argument"];
   AppendTo[checks, assertCheck["Ansatz.couplingsOutsideOperators", AllTrue[ops, FreeQ[First[#], Alternatives @@ cps] &],
     Select[ops, ! FreeQ[First[#], Alternatives @@ cps] &],
     "couplings depend on k only: write nk dt[H], not dt[nk H] (and never differentiate a coupling)", "no coupling under dt/grad/div/lap"]];
   FRGAssert["Ansatz.operatorOnField", AllTrue[ops, ! FreeQ[First[#], fld] &], Select[ops, FreeQ[First[#], fld] &],
    "an operator acts on an expression without fields"];

   FRGAssert["Ansatz.scalarDensity", rank[dens, fld] === 0, dens, "the density must be a scalar (contract gradients with a dot)"];
   new[] := SpaceIndex[{"dummy", ++counter}];
   expanded = Expand[toFactors[dens, fld, new]];
   FRGAssert["Ansatz.fieldsInFactors", FreeQ[expanded /. _FieldFactor -> 1, fld | opHeads | Dot], expanded,
    "some field or operator was not converted (internal error, or a non-polynomial structure such as Exp[H])"];
   terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
   mons = splitTerm /@ terms;
   mons = KeyValueMap[{Together[Total[#2]], #1} &, GroupBy[mons, canonical[#[[2]]] & -> First]];
   If[MemberQ[mons[[All, 1]], 0],
    AppendTo[warns, infoCheck["INFO", "Ansatz.cancelled", ToString[Count[mons[[All, 1]], 0]] <> " monomial(s) cancel and are dropped"]]];
   mons = monomialRecord[#, fields, responses] & /@ DeleteCases[mons, {0, _}];
   mons = SortBy[mons, {#FieldCount, Sort[First /@ (FirstPosition[fields, #] & /@ #Fields)], #TimeOrder, #MomentumOrder} &];

   AppendTo[checks, assertCheck["Ansatz.MSR", AllTrue[mons, #ResponseCount >= 1 &], MonomialString /@ Select[mons, #ResponseCount == 0 &],
     "every monomial must contain at least one response field (MSR action)", "every monomial has >= 1 response field"]];
   AppendTo[checks, assertCheck["Ansatz.noLinear", AllTrue[mons, #FieldCount >= 2 &], MonomialString /@ Select[mons, #FieldCount < 2 &],
     "terms linear in the fields are not supported (expand around a stationary point)", "no linear term"]];
   AppendTo[checks, assertCheck["Ansatz.quartic", AllTrue[mons, #FieldCount <= 4 &], MonomialString /@ Select[mons, #FieldCount > 4 &],
     "v1 supports at most quartic monomials", "at most quartic (max " <> ToString[Max[#FieldCount & /@ mons]] <> ")"]];
   AppendTo[checks, assertCheck["Ansatz.momentumOrder", AllTrue[mons, #MomentumOrder <= 2 &], MonomialString /@ Select[mons, #MomentumOrder > 2 &],
     "momentum order p^4 or higher: the p-Laplacian projection is valid for p^2 only (E2, see docs/FUTURE.md)", "momentum order <= p^2"]];
   FRGAssert["Ansatz.indexPairs", AllTrue[mons, AllTrue[Tally[Flatten[#Factors[[All, 3]]]][[All, 2]], EqualTo[2]] &],
    Select[mons, ! AllTrue[Tally[Flatten[#Factors[[All, 3]]]][[All, 2]], EqualTo[2]] &], "internal: every spatial index must appear exactly twice"];
   AppendTo[checks, assertCheck["Ansatz.couplingsInCoefficients", AllTrue[mons, FreeQ[#Coefficient, fld | opHeads | SpaceIndex] &], mons,
     "internal: a coefficient contains a field", "coefficients are field-free"]];
   unused = Select[cps, FreeQ[mons[[All, "Coefficient"]], #] &];
   If[unused =!= {}, AppendTo[warns, infoCheck["WARN", "Ansatz.unusedCouplings", "declared but absent from the density: " <> ToString[unused]]]];

   stage = stampStage[<|"Stage" -> "Ansatz", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {theory["Fingerprint"]},
      "Theory" -> theory, "Density" -> dens, "Couplings" -> cps, "Monomials" -> mons,
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowMonomials[stage];
   PrintCard[stage];
   stage]];

ShowMonomials[ansatz_Association] := Module[{mons = validateStage[ansatz, "Ansatz"]["Monomials"], fo = ansatz["Theory"]["Fields"], rows},
   rows = {inputForm[#Coefficient], ToString[#Fields], MonomialString[#],
       "dt^" <> ToString[#TimeOrder] <> " p^" <> ToString[#MomentumOrder]} & /@ mons;
   Print["Parsed density: ", Length[mons], " monomials   (Gamma_k = Integral dt d^dx [ sum of coefficient x operators ])"];
   Print[If[$FrontEnd === Null,
     StringRiffle[StringRiffle[MapThread[StringPadRight, {#, {22, 26, 34, 10}}], " | "] & /@ Prepend[rows, {"coefficient", "fields", "operators", "order"}], "\n"],
     Grid[Prepend[rows, {"coefficient", "fields", "operators", "order"}], Frame -> All, Alignment -> Left]]];];

End[];
EndPackage[];
