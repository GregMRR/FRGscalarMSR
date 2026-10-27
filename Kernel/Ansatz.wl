(* ::Package:: *)

(* Ansatz.wl -- DefineAnsatz: parse an MSR action density written with the package operators dt, grad, div, lap into a
   table of monomials  coefficient x product of (derivatives of single fields).
   Objective, assumptions and checks: docs/MODULES.md section "Ansatz". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

dt::usage = "dt[e]: time derivative in the ansatz density (package-owned and inert; not System`Dt).";
grad::usage = "grad[e]: spatial gradient (a vector) in the ansatz density; contract two of them with a dot: grad[a].grad[b].";
div::usage = "div[v]: divergence of a vector expression in the ansatz density, e.g. div[H grad[Psi]].";
lap::usage = "lap[e]: Laplacian in the ansatz density.";
jac::usage = "jac[a, b] = eps_ij d_i a d_j b (eps_12 = +1) in the ansatz density, a and b scalars. It exists only in d = 2: declare \
\"Dimension\" -> FixedDimension[2] in DefineTheory. At most one jac per monomial.";
DefineAnsatz::usage = "DefineAnsatz[theory, \"Density\" -> density, \"Couplings\" -> {c1, c2, ...}] parses the density \
(Gamma_k = Integral dt d^dx density) into monomials, prints the parsed table and returns the \"Ansatz\" stage.";
FieldFactor::usage = "FieldFactor[phi, nt, {i, ...}] is the factor dt^nt d_i ... phi of a parsed monomial; every repeated \
SpaceIndex[n] is summed over (one index per contraction, numbered canonically).";
SpaceIndex::usage = "SpaceIndex[n] is a summed spatial index in a parsed monomial.";
EpsIndex::usage = "EpsIndex[s] (s = 1, 2) is the s-th index of eps_ij in a parsed monomial: jac[a, b] gives EpsIndex[1] on the factor of a \
and EpsIndex[2] on the factor of b.";
ShowMonomials::usage = "ShowMonomials[ansatz] prints the table of parsed monomials (coefficient | fields | operators).";
MonomialString::usage = "MonomialString[monomial] writes the field part of a parsed monomial back in the DSL.";
fn::usage = "fn[f, X] in the ansatz density is f(-Laplacian) acting on X: X is one field, optionally under grad/lap, and f a function \
declared in DefineAnsatz \"Functions\" -> {f -> reference}. In Fourier space f takes the momentum of X: in a bilinear term \
nuk fn[fnu, lap[h]] gives nuk p^2 fnu(p^2). An ansatz with fn uses the LO pipeline (design/LOPIPELINE.md).";
MonomialFunction::usage = "MonomialFunction[monomial] is the function f of the fn factor of a parsed monomial, or None.";

Begin["`Private`"];

SetAttributes[{dt, grad, div, lap, jac, fn, FieldFactor, SpaceIndex, EpsIndex}, Protected];

opHeads = dt | grad | div | lap | jac;
systemOps = Dt | D | Grad | Div | Laplacian | Curl | Derivative | Integrate | Inactive | Inactivate;

(* ---------- vector rank of a raw density expression (0 scalar, 1 vector); fails on ill-formed input ---------- *)

rank[e_, fld_] := Which[
   FreeQ[e, fld | opHeads], 0,
   MatchQ[e, fld], 0,
   True, rankOf[e, fld]];
rankOf[grad[a_], fld_] := (needRank[a, 0, fld, "grad[...] acts on a scalar"]; 1);
rankOf[dt[a_], fld_] := rank[a, fld];
rankOf[fn[f_, a_], fld_] := rank[a, fld];
rankOf[div[a_], fld_] := (needRank[a, 1, fld, "div[...] acts on a vector, e.g. div[H grad[Psi]]"]; 0);
rankOf[lap[a_], fld_] := (needRank[a, 0, fld, "lap[...] acts on a scalar"]; 0);
rankOf[jac[a_, b_], fld_] := (needRank[a, 0, fld, "jac[a, b] acts on two scalars"]; needRank[b, 0, fld, "jac[a, b] acts on two scalars"]; 0);
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

(* ---------- conversion to products of FieldFactor[phi, nt, {indices}] (LO: FieldFactor[phi, nt, {indices}, f]) ---------- *)

fnTag[e_, f_] := e /. FieldFactor[a_, nt_, ix_] :> FieldFactor[a, nt, ix, f];
MonomialFunction[m_Association] := FirstCase[m["Factors"], FieldFactor[_, _, _, f_] :> f, None];

(* the argument X of fn[f, X]: one field, optionally under grad/lap *)
fnArgumentQ[x_, fld_] := MatchQ[x, fld] || MatchQ[x, (grad | lap)[y_] /; fnArgumentQ[y, fld]];

(* derivations acting on a sum of products of FieldFactors *)
derive[op_][e_Plus] := derive[op] /@ e;
derive[op_][e_Times] := Sum[MapAt[derive[op], e, i], {i, Length[e]}];
derive[op_][Power[f_, n_Integer]] := n f^(n - 1) derive[op][f];
derive[op_][FieldFactor[f_, nt_, ix_, fnf___]] := Replace[op, {"t" -> FieldFactor[f, nt + 1, ix, fnf], i_ :> FieldFactor[f, nt, Append[ix, i], fnf]}];
derive[op_][c_] /; FreeQ[c, FieldFactor] := 0;

toFactors[e_, fld_, new_] := Module[{S, V},
   S[x_] /; FreeQ[x, fld | opHeads] := x;
   S[x : fld] := FieldFactor[x, 0, {}];
   S[x_Plus] := S /@ x;
   S[x_Times] := S /@ x;
   (* each factor of a power gets its own dummies: lap[h]^2, (grad h.grad h)^2 (design/MOMENTUMPROJECTION.md, review E14) *)
   S[Power[x_, n_Integer?Positive]] := Expand[Times @@ Table[S[x], {n}]];
   S[Power[x_, n_]] := Expand[S[x]^n];
   S[dt[x_]] := Expand[derive["t"][S[x]]];
   (* fn[f, x]: the factors of x (one field) carry f in a 4th slot (LO pipeline, design/LOPIPELINE.md L3) *)
   S[fn[f_, x_]] := fnTag[S[x], f];
   V[fn[f_, x_], i_] := fnTag[V[x, i], f];
   S[lap[x_]] := With[{j = new[]}, Expand[derive[j][derive[j][S[x]]]]];
   S[div[v_]] := With[{j = new[]}, Expand[derive[j][V[v, j]]]];
   S[Dot[u_, v_]] := With[{j = new[]}, Expand[V[u, j] V[v, j]]];
   S[jac[a_, b_]] := With[{j = new[]}, Expand[derive[EpsIndex[j, 1]][S[a]] derive[EpsIndex[j, 2]][S[b]]]];
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

(* canonical form of a monomial, so that equal monomials combine; returns {sign, factors}. For every numbering of the summed
   indices as SpaceIndex[1..]: sort each factor's index list (derivatives commute), sort the factor list; keep the smallest
   (a class invariant for any number of contractions: lap[lap[A]] = div[grad[lap[A]]], review E9).
   With eps (one pair EpsIndex[1], EpsIndex[2]): eps_ij is antisymmetric, so the factor list with the slots swapped is the same
   monomial with sign -1; the smaller of the two minima is kept, and equal minima mean zero (sign 0). *)
canonMinForm[factors_List, ix_List] := First[Sort[Function[perm,
     Sort[(factors /. Thread[ix -> (SpaceIndex /@ perm)]) /. FieldFactor[f_, nt_, i_, fnf___] :> FieldFactor[f, nt, Sort[i], fnf]]] /@ Permutations[Range[Length[ix]]]]];
canonical[factors0_List] := Module[{factors = factors0 /. EpsIndex[_, s_] :> EpsIndex[s], ix, l, ls},
   ix = Union[Cases[Flatten[factors[[All, 3]]], Except[_EpsIndex]]];
   l = canonMinForm[factors, ix];
   If[FreeQ[l, EpsIndex], {1, l},
    ls = canonMinForm[factors /. {EpsIndex[1] -> EpsIndex[2], EpsIndex[2] -> EpsIndex[1]}, ix];
    Which[ls === l, {0, l}, OrderedQ[{l, ls}], {1, l}, True, {-1, ls}]]];

arity[_jac] = 2;
arity[_] = 1;

splitTerm[t_] := Module[{fs = If[Head[t] === Times, List @@ t, {t}]},
   FRGAssert["Ansatz.polynomial", AllTrue[fs, MatchQ[#, _FieldFactor | Power[_FieldFactor, _Integer?Positive] | _?(FreeQ[#, FieldFactor] &)] &], t,
    "the density must be polynomial in the fields"];
   {Times @@ Select[fs, FreeQ[#, FieldFactor] &],
    Flatten[Replace[Select[fs, ! FreeQ[#, FieldFactor] &], {Power[f_, n_] :> Table[f, n], f_ :> {f}}, {1}]]}];

fieldString[FieldFactor[f_, nt_, ix_, fnf___]] := Nest["dt[" <> # <> "]" &,
   If[{fnf} === {}, #, "fn[" <> ToString[fnf] <> ", " <> # <> "]"] &[
    Which[ix === {}, ToString[f], EvenQ[Length[ix]] && AllTrue[Tally[ix][[All, 2]], EvenQ], Nest["lap[" <> # <> "]" &, ToString[f], Length[ix]/2],
     Length[ix] == 1, "grad[" <> ToString[f] <> "]", True, "d" <> ToString[ix] <> "[" <> ToString[f] <> "]"]], nt];
(* factors carrying one index are the two ends of the contraction: printed as grad[a].grad[b]; the factors carrying EpsIndex are
   printed as jac[slot 1, slot 2] *)
MonomialString[m_Association] := Module[{fs0 = m["Factors"], fs, pair, eps},
   eps = SortBy[Select[fs0, ! FreeQ[#[[3]], EpsIndex] &], First[Cases[#[[3]], EpsIndex[s_] :> s]] &];
   fs = Select[fs0, FreeQ[#[[3]], EpsIndex] &];
   pair = Select[fs, Length[#[[3]]] == 1 &];
   StringRiffle[Join[fieldString /@ Select[fs, Length[#[[3]]] =!= 1 &], If[pair === {}, {}, {StringRiffle[fieldString /@ pair, "."]}],
     If[eps === {}, {}, {"jac[" <> StringRiffle[fieldString[MapAt[DeleteCases[#, _EpsIndex] &, #, 3]] & /@ eps, ", "] <> "]"}]], " "]];

(* every SpaceIndex appears twice; EpsIndex[1] and EpsIndex[2] once each, or neither *)
indexPairsQ[m_] := Module[{ix = Flatten[m["Factors"][[All, 3]]]},
   AllTrue[Tally[Cases[ix, _SpaceIndex]][[All, 2]], EqualTo[2]] && MemberQ[{{0, 0}, {1, 1}}, {Count[ix, EpsIndex[1]], Count[ix, EpsIndex[2]]}]];

(* factors stored response fields first, then in field order (display order; identity was fixed by canonical) *)
monomialRecord[{coef_, fac0_}, fieldOrder_, responses_] := Module[{nIx = Length[Flatten[fac0[[All, 3]]]], fac},
   fac = SortBy[fac0, {! MemberQ[responses, #[[1]]], First[FirstPosition[fieldOrder, #[[1]]]], #[[2]], Length[#[[3]]]} &];
   <|"Coefficient" -> coef, "Factors" -> fac, "Fields" -> fac[[All, 1]],
    "FieldCount" -> Length[fac], "ResponseCount" -> Count[fac[[All, 1]], Alternatives @@ responses],
    "TimeOrder" -> Total[fac[[All, 2]]], "MomentumOrder" -> nIx|>];

(* ---------- DefineAnsatz ---------- *)

Options[DefineAnsatz] = {"Density" -> None, "Couplings" -> {}, "Pipeline" -> Automatic, "Functions" -> {}};

DefineAnsatz[theory_, opts : OptionsPattern[]] := Block[{$FRGModule = "Ansatz"},
  Module[{t0 = AbsoluteTime[], dens, cps, fields, fld, responses, checks = {}, warns = {}, counter = 0, new, ops, expanded, terms, mons,
    badSyms, unused, stage, jacQ, epsQ, loQ, pipeline, funcs, fsyms, frefs, fnOps, bad, kinds},
   validateStage[theory, "Theory"];
   FRGAssert["Ansatz.options", MatchQ[{opts}, {(Rule | RuleDelayed)["Density" | "Couplings" | "Pipeline" | "Functions", _] ...}], {opts},
    "only the options \"Density\", \"Couplings\", \"Pipeline\", \"Functions\" are accepted"];
   dens = OptionValue["Density"]; cps = OptionValue["Couplings"];
   fields = theory["Fields"]; responses = theory["ResponseFields"]; fld = Alternatives @@ fields;
   (* the pipeline is inferred: fn[...] in the density means couplings that are functions of p^2 (LO, design/LOPIPELINE.md L4) *)
   loQ = ! FreeQ[dens, fn];
   pipeline = If[loQ, "LO", "DerivativeExpansion"];
   AppendTo[checks, assertCheck["Ansatz.pipeline", MemberQ[{Automatic, pipeline}, OptionValue["Pipeline"]], {OptionValue["Pipeline"], "inferred" -> pipeline},
     "\"Pipeline\" is inferred from the density (fn[...] present: \"LO\", otherwise \"DerivativeExpansion\"); an explicit value must match it", pipeline]];
   funcs = OptionValue["Functions"];
   FRGAssert["Ansatz.functionReference", MatchQ[funcs, {(_Symbol -> _Symbol) ...}], funcs, "\"Functions\" -> {f1 -> reference1, ...}"];
   fsyms = funcs[[All, 1]]; frefs = funcs[[All, 2]];
   If[! loQ, FRGAssert["Ansatz.functionReference", funcs === {}, funcs, "\"Functions\" are declared but the density has no fn[f, X]"]];
   FRGAssert["Ansatz.density", dens =!= None && dens =!= 0, dens, "\"Density\" -> expression is required"];
   FRGAssert["Ansatz.couplingList", ListQ[cps], cps, "\"Couplings\" -> {c1, c2, ...}"];
   assertFreshSymbols["Ansatz.couplings", cps, "couplings"];
   FRGAssert["Ansatz.couplingsNotFields", Intersection[cps, Append[fields, theory["Dimension"]]] === {}, Intersection[cps, Append[fields, theory["Dimension"]]],
    "a coupling has the name of a field or of the dimension"];

   (* E7: System derivative heads would silently differentiate the couplings (Dt[nk H] = H Dt[nk] + nk Dt[H]) *)
   AppendTo[checks, assertCheck["Ansatz.noSystemDerivatives", FreeQ[dens, systemOps], Cases[dens, systemOps[___] | Derivative[__][_], {0, Infinity}, Heads -> True],
     "use the package operators dt, grad, div, lap (System Dt, Grad, Div, Laplacian, D are not allowed)", "no Dt/D/Grad/Div/Laplacian"]];

   (* functions (LO): fresh symbols tied to distinct references; every fn[f, X] has a declared f and X = one field under grad/lap *)
   If[loQ,
    assertFreshSymbols["Ansatz.functions", fsyms, "functions"];
    fnOps = Cases[dens, _fn, {0, Infinity}];
    AppendTo[checks, assertCheck["Ansatz.functionReference",
      SubsetQ[cps, frefs] && DuplicateFreeQ[frefs] && Intersection[fsyms, Join[cps, fields, {theory["Dimension"]}]] === {} &&
       SubsetQ[fsyms, Cases[fnOps, fn[f_, ___] :> f]],
      <|"functions" -> funcs, "undeclared" -> Complement[Cases[fnOps, fn[f_, ___] :> f], fsyms]|>,
      "each fn[f, X] needs f declared in \"Functions\" -> {f -> reference}, the reference a coupling of the ansatz, one reference per function \
(f(0) = 1 fixes the eta of its reference)", "functions " <> StringRiffle[(ToString[#[[1]]] <> " -> " <> ToString[#[[2]]]) & /@ funcs, ", "]]];
    AppendTo[checks, assertCheck["Ansatz.fnOneField", AllTrue[fnOps, MatchQ[#, fn[_Symbol, x_] /; fnArgumentQ[x, fld]] &],
      Select[fnOps, ! MatchQ[#, fn[_Symbol, x_] /; fnArgumentQ[x, fld]] &],
      "fn[f, X]: X is one field, optionally under grad/lap (e.g. fn[fnu, lap[h]]); no coupling, product or fn inside", "fn[f, X] acts on one field"]]];

   badSyms = Complement[Cases[dens, s_Symbol /; Context[s] =!= "System`", {0, Infinity}, Heads -> True], fields, cps, fsyms, {theory["Dimension"], dt, grad, div, lap, jac, fn}];
   AppendTo[checks, assertCheck["Ansatz.knownSymbols", badSyms === {}, badSyms,
     If[MemberQ[SymbolName /@ badSyms, "t" | "x"], "couplings depend on k only: no explicit t or x dependence is allowed",
      "every symbol must be a declared coupling (\"Couplings\") or a field of the theory"], "every symbol is a field, a coupling or Dim"]];

   (* operators act on field expressions only: couplings depend on k only and must stand outside (no dt[nk H], no lap[nk]) *)
   ops = Cases[dens, opHeads[___], {0, Infinity}];
   FRGAssert["Ansatz.operatorArity", AllTrue[ops, Length[#] == arity[#] &], Select[ops, Length[#] =!= arity[#] &],
    "dt, grad, div, lap take exactly one argument; jac takes two"];
   AppendTo[checks, assertCheck["Ansatz.couplingsOutsideOperators", AllTrue[ops, FreeQ[List @@ #, Alternatives @@ cps] &],
     Select[ops, ! FreeQ[List @@ #, Alternatives @@ cps] &],
     "couplings depend on k only: write nk dt[H], not dt[nk H] (and never differentiate a coupling)", "no coupling under dt/grad/div/lap"]];
   FRGAssert["Ansatz.operatorOnField", AllTrue[ops, AllTrue[List @@ #, ! FreeQ[#, fld] &] &], Select[ops, ! AllTrue[List @@ #, ! FreeQ[#, fld] &] &],
    "an operator acts on an expression without fields"];
   jacQ = ! FreeQ[dens, jac];
   If[jacQ, AppendTo[checks, assertCheck["Ansatz.epsilonDimension", theory["Dimension"] === 2, theory["Dimension"],
      "jac (eps_ij) exists only in d = 2: declare \"Dimension\" -> FixedDimension[2]", "jac with d = 2 (declared)"]]];

   FRGAssert["Ansatz.scalarDensity", rank[dens, fld] === 0, dens, "the density must be a scalar (contract gradients with a dot)"];
   new[] := SpaceIndex[{"dummy", ++counter}];
   expanded = Expand[toFactors[dens, fld, new]];
   FRGAssert["Ansatz.fieldsInFactors", FreeQ[expanded /. _FieldFactor -> 1, fld | opHeads | Dot], expanded,
    "some field or operator was not converted (internal error, or a non-polynomial structure such as Exp[H])"];
   terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
   mons = splitTerm /@ terms;
   If[jacQ, AppendTo[checks, assertCheck["Ansatz.singleEpsilon", AllTrue[mons, Count[#[[2, All, 3]], EpsIndex[_, 1], Infinity] <= 1 &],
      Select[mons, Count[#[[2, All, 3]], EpsIndex[_, 1], Infinity] > 1 &], "at most one jac per monomial (the eps reduction handles one pair of eps indices)", "at most one jac per monomial"]]];
   mons = Function[m, With[{c = canonical[m[[2]]]}, {c[[1]] m[[1]], c[[2]]}]] /@ mons;
   mons = KeyValueMap[{Together[Total[#2]], #1} &, GroupBy[mons, Last -> First]];
   If[MemberQ[mons[[All, 1]], 0],
    AppendTo[warns, infoCheck["INFO", "Ansatz.cancelled", ToString[Count[mons[[All, 1]], 0]] <> " monomial(s) cancel and are dropped"]]];
   mons = monomialRecord[#, fields, responses] & /@ DeleteCases[mons, {0, _}];
   mons = SortBy[mons, {#FieldCount, Sort[First /@ (FirstPosition[fields, #] & /@ #Fields)], #TimeOrder, #MomentumOrder} &];
   epsQ = ! FreeQ[mons[[All, "Factors"]], EpsIndex];
   If[epsQ, AppendTo[checks, assertCheck["Ansatz.jacQuadratic", NoneTrue[mons, #FieldCount == 2 && ! FreeQ[#Factors, EpsIndex] &],
      MonomialString /@ Select[mons, #FieldCount == 2 && ! FreeQ[#Factors, EpsIndex] &],
      "a quadratic jac term is a total derivative and does not contribute: remove it", "no quadratic jac term (total derivative)"]]];

   AppendTo[checks, assertCheck["Ansatz.MSR", AllTrue[mons, #ResponseCount >= 1 &], MonomialString /@ Select[mons, #ResponseCount == 0 &],
     "every monomial must contain at least one response field (MSR action)", "every monomial has >= 1 response field"]];
   AppendTo[checks, assertCheck["Ansatz.noLinear", AllTrue[mons, #FieldCount >= 2 &], MonomialString /@ Select[mons, #FieldCount < 2 &],
     "terms linear in the fields are not supported (expand around a stationary point)", "no linear term"]];
   AppendTo[checks, assertCheck["Ansatz.quartic", AllTrue[mons, #FieldCount <= 4 &], MonomialString /@ Select[mons, #FieldCount > 4 &],
     "v1 supports at most quartic monomials", "at most quartic (max " <> ToString[Max[#FieldCount & /@ mons]] <> ")"]];
   FRGAssert["Ansatz.indexPairs", AllTrue[mons, indexPairsQ], Select[mons, ! indexPairsQ[#] &],
    "internal: every spatial index must appear exactly twice, and EpsIndex[1], EpsIndex[2] once each or not at all"];
   If[loQ, checks = Join[checks, loMonomialChecks[mons, funcs, cps, fld]]];
   AppendTo[checks, assertCheck["Ansatz.couplingsInCoefficients", AllTrue[mons, FreeQ[#Coefficient, fld | opHeads | SpaceIndex | EpsIndex] &], mons,
     "internal: a coefficient contains a field", "coefficients are field-free"]];
   unused = Select[cps, FreeQ[mons[[All, "Coefficient"]], #] &];
   If[unused =!= {}, AppendTo[warns, infoCheck["WARN", "Ansatz.unusedCouplings", "declared but absent from the density: " <> ToString[unused]]]];

   (* LO-only keys: a DE stage keeps exactly its former keys (design/LOPIPELINE.md step 0) *)
   kinds = AssociationThread[cps -> (If[MemberQ[frefs, #], "Function", "Constant"] & /@ cps)];
   stage = stampStage[<|"Stage" -> "Ansatz", "Pipeline" -> pipeline, "Upstream" -> {theory["Fingerprint"]},
      "Theory" -> theory, "Density" -> dens, "Couplings" -> cps, "Monomials" -> mons, Sequence @@ If[epsQ, {"Epsilon" -> True}, {}],
      Sequence @@ If[loQ, {"Functions" -> funcs, "CouplingKinds" -> kinds}, {}],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowMonomials[stage];
   PrintCard[stage];
   stage]];

(* LO checks on the parsed monomials (design/LOPIPELINE.md L1-L4) *)
loMonomialChecks[mons_, funcs_, cps_, fld_] := Module[{withF = Select[mons, MonomialFunction[#] =!= None &], bil, bad, out = {}, refOf = Association[funcs]},
   AppendTo[out, assertCheck["Ansatz.fnBilinear", AllTrue[withF, #FieldCount == 2 &], MonomialString /@ Select[withF, #FieldCount =!= 2 &],
     "functions of p^2 only in bilinear (Gamma2) terms: vertex couplings are k-only constants (L1)", "functions in bilinear terms only"]];
   bad = Select[withF, Count[#Factors, FieldFactor[_, _, _, _]] > 1 &];
   AppendTo[out, assertCheck["Ansatz.fnOneField", bad === {}, MonomialString /@ bad, "at most one fn[...] per monomial", "one fn per monomial"]];
   AppendTo[out, assertCheck["Ansatz.functionTimesOmega", AllTrue[withF, #TimeOrder == 0 &], MonomialString /@ Select[withF, #TimeOrder > 0 &],
     "a function times i w (dt) is not in LO: only the w^0 coefficients are functions (L2; the d/d(i nu) extension is in docs/FUTURE.md)",
     "no function times i w"]];
   bad = Select[withF, ! (KeyExistsQ[refOf, MonomialFunction[#]] && NumericQ[Together[#Coefficient/refOf[MonomialFunction[#]]]]) &];
   AppendTo[out, assertCheck["Ansatz.functionReference", bad === {} && DuplicateFreeQ[MonomialFunction /@ withF] && Sort[MonomialFunction /@ withF] === Sort[Keys[refOf]],
     <|"coefficient not number x reference" -> (MonomialString /@ bad), "functions used" -> (MonomialFunction /@ withF)|>,
     "each function appears in exactly one monomial, with coefficient (number) x its reference (\"Functions\" -> {f -> reference})",
     "each function in one monomial, coefficient number x reference"]];
   bil = Select[mons, #FieldCount == 2 &];
   bad = Select[bil, (#TimeOrder == 0 && MonomialFunction[#] === None) || (#TimeOrder > 0 && ! NumericQ[#Coefficient]) &];
   AppendTo[out, assertCheck["Ansatz.mixedNotYet", bad === {}, MonomialString /@ bad,
     "an LO ansatz writes every w^0 bilinear term with a function and every i w term with a numeric coefficient; mixed ansatze \
(DE monomials + functions) are not implemented yet (docs/FUTURE.md)", "bilinear terms: functions at w^0, numbers at i w"]];
   out];

ShowMonomials[ansatz_Association] := Module[{mons = validateStage[ansatz, "Ansatz"]["Monomials"], fo = ansatz["Theory"]["Fields"], rows},
   rows = {inputForm[#Coefficient], ToString[#Fields], MonomialString[#],
       "dt^" <> ToString[#TimeOrder] <> " p^" <> ToString[#MomentumOrder]} & /@ mons;
   Print["Parsed density: ", Length[mons], " monomials   (Gamma_k = Integral dt d^dx [ sum of coefficient x operators ])"];
   Print[If[$FrontEnd === Null,
     StringRiffle[StringRiffle[MapThread[StringPadRight, {#, {22, 26, 34, 10}}], " | "] & /@ Prepend[rows, {"coefficient", "fields", "operators", "order"}], "\n"],
     Grid[Prepend[rows, {"coefficient", "fields", "operators", "order"}], Frame -> All, Alignment -> Left]]];];

End[];
EndPackage[];
