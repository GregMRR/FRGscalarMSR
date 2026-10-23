(* ::Package:: *)

(* Vertices.wl -- BuildVertices: Fourier transform of the parsed density, vertices Gamma^(n) as symmetrised functional
   derivatives, the Gamma2 matrix and the Gamma3 matrices in the leg convention of docs/CONVENTIONS.md.
   Objective, assumptions and checks: docs/MODULES.md section "Vertices". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

BuildVertices::usage = "BuildVertices[ansatz] Fourier-transforms every monomial (dt -> +I w, grad -> -I q, all momenta incoming), \
builds every nonzero vertex Gamma^(n), n = 2..4, checks it against an independent plane-wave substitution into the density, \
prints the vertex table with its convention header and returns the \"Vertices\" stage.";
VertexValue::usage = "VertexValue[vertices, {a1, ..., an}, {{w1, q1}, ..., {w(n-1), q(n-1)}}] is Gamma^(n)_{a1...an} with leg j at \
(wj, qj) and the last leg at -(sum of the others); vectors are symbols, scalar products scal[a, b].";
Gamma2Matrix::usage = "Gamma2Matrix[vertices, w, q] is the matrix Gamma2_ab(w, q) (leg a at (w, q), leg b at (-w, -q)) in field order, \
without regulator.";
Gamma3Matrix::usage = "Gamma3Matrix[vertices, i, {nu, p}, {w, q}] is the matrix with entries Gamma^(3)_{i a b}: leg i at (nu, p), \
leg a at (w, q), leg b at -(nu + w, p + q); i is a field or its position in the field order.";
ShowVertex::usage = "ShowVertex[vertices, {a1, ..., an}] prints Gamma^(n)_{a1...an}(w1, q1; ...; w(n-1), q(n-1)), last leg = -(sum).";
ShowGamma2::usage = "ShowGamma2[vertices] prints the nonzero entries Gamma2_ab(w, q) (without regulator).";
ShowVertexTable::usage = "ShowVertexTable[vertices] prints the convention header and every nonzero vertex (done by BuildVertices).";

Begin["`Private`"];

(* ---------- Fourier symbol of a monomial at a given assignment of legs ---------- *)

(* leg l carries frequency wv[l] and momentum qv[l]; a factor dt^nt d_i.. phi at leg l gives (I wv[l])^nt (-I qv[l]_i)... *)
factorSymbol[FieldFactor[_, nt_, ix_], l_] := {(I wv[l])^nt (-I)^Length[ix], Thread[ix -> qv[l]]};

(* a SpaceIndex pair gives scal[qa, qb]; the eps pair EpsIndex[1], EpsIndex[2] gives wedge[q(slot 1), q(slot 2)] *)
monomialAt[m_Association, perm_List] := Module[{parts = MapThread[factorSymbol, {m["Factors"], perm}], g},
   g = GroupBy[Flatten[parts[[All, 2]]], First -> Last];
   m["Coefficient"] Times @@ parts[[All, 1]] Times @@ (scal @@ # & /@ Values[KeySelect[g, FreeQ[#, EpsIndex] &]]) *
    If[KeyExistsQ[g, EpsIndex[1]], wedge[First[g[EpsIndex[1]]], First[g[EpsIndex[2]]]], 1]];

(* Gamma^(n) for the field multiset ms (sorted in field order), all n momenta independent (before delta stripping):
   sum over every assignment of the monomial's factors to the legs carrying the same field = full symmetrisation *)
vertexExpr[ms_List, mons_List, sortF_] := Expand[Total[Function[m,
      Total[monomialAt[m, #] & /@ Select[Permutations[Range[Length[ms]]], m["Fields"] === ms[[#]] &]]] /@
     Select[mons, #FieldCount == Length[ms] && sortF[#Fields] === ms &]]];

(* ---------- evaluation in the leg convention ---------- *)

fieldPos[v_, a_] := First[FirstPosition[v["Fields"], a, {0}]];

VertexValue[v_Association, tuple_List, legs_List] := Module[{full, ord, ms, e},
   FRGAssert["Vertices.legCount", Length[legs] == Length[tuple] - 1 && AllTrue[legs, Length[#] == 2 &], {tuple, legs},
    "give n - 1 pairs {w, q} for an n-point vertex (the last leg carries minus their sum)"];
   FRGAssert["Vertices.knownFields", AllTrue[tuple, fieldPos[v, #] > 0 &], tuple, "unknown field"];
   full = Append[legs, -Total[legs]];
   ord = Ordering[fieldPos[v, #] & /@ tuple];   (* canonical leg j is original leg ord[[j]]; ties are identical fields *)
   ms = tuple[[ord]];
   e = Lookup[v["VertexTable"], Key[ms], 0];
   Expand[e /. {wv[j_] :> full[[ord[[j]], 1]], qv[j_] :> full[[ord[[j]], 2]]}]];

Gamma2Matrix[v_Association, w_, q_] := Table[VertexValue[v, {a, b}, {{w, q}}], {a, v["Fields"]}, {b, v["Fields"]}];

Gamma3Matrix[v_Association, i0_, {nu_, p_}, {w_, q_}] := Module[{i = If[IntegerQ[i0], v["Fields"][[i0]], i0]},
   Table[VertexValue[v, {i, a, b}, {{nu, p}, {w, q}}], {a, v["Fields"]}, {b, v["Fields"]}]];

(* ---------- independent check: plane waves substituted into the raw density ---------- *)

(* phi_a(t, x) = sum over the legs l of field a of eps[l] Exp[I (ww[l] t - kq[l].X)]  (nc = 3 components, 2 for a declared d = 2);
   the coefficient of eps[1]...eps[n] Exp[I sum(ww t - kq.X)] in the density is Gamma^(n)(all legs) *)
planeWaveCoefficient[dens_, fields_, ms_List, nc_] := Module[{n = Length[ms], X = Array[xx, nc], kq, phase, rules, e},
   kq[l_] := Table[kc[l, c], {c, nc}];
   phase = I Sum[ww[l] tt - kq[l] . X, {l, n}];
   rules = Thread[fields -> (Sum[If[ms[[l]] === #, eps[l] Exp[I (ww[l] tt - kq[l] . X)], 0], {l, n}] & /@ fields)];
   e = (dens /. rules) /. {dt -> (D[#, tt] &), grad -> (Grad[#, X] &), div -> (If[# === 0, 0, Div[#, X]] &), lap -> (Laplacian[#, X] &),
      jac -> (D[#1, xx[1]] D[#2, xx[2]] - D[#1, xx[2]] D[#2, xx[1]] &)};
   e = D[e, Sequence @@ Table[eps[l], {l, n}]] /. eps[_] -> 0;
   e = Expand[e] /. Power[E, z_] :> Exp[Expand[z - phase]];
   FRGAssert["Vertices.planeWavePhase", FreeQ[e, tt | xx], e, "internal: the plane-wave coefficient still depends on t or x"];
   Expand[e]];

symbolAsPlaneWave[e_, nc_] := Module[{vec = # /. qv[l_] :> Table[kc[l, c], {c, nc}] &},
   Expand[e /. {scal[a_, b_] :> vec[a] . vec[b], wedge[a_, b_] :> vec[a][[1]] vec[b][[2]] - vec[a][[2]] vec[b][[1]]} /. wv[l_] :> ww[l]]];

(* ---------- display ---------- *)

legSymbols[n_] := Table[{Symbol["FRGscalarMSR`Vertices`Private`\[Omega]" <> ToString[j]], Symbol["FRGscalarMSR`Vertices`Private`q" <> ToString[j]]}, {j, n}];
vertexName[tuple_] := "Gamma(" <> ToString[Length[tuple]] <> ")_{" <> StringRiffle[ToString /@ tuple, " "] <> "}";
legString[legs_] := "(" <> StringRiffle[displayString[#[[1]]] <> "," <> displayString[#[[2]]] & /@ legs, "; "] <> ")";

conventionHeader = {
   "Fourier convention: phi(t,x) = Int dw d^dq/(2Pi)^(d+1) Exp[-I (q.x - w t)] phi(w,q)   =>   dt -> +I w,  grad -> -I q",
   "Functional derivative: d phi(w,q)/d phi(w',q') = (2Pi)^(d+1) delta(w-w') delta(q-q')",
   "d^n Gamma/d phi_a1(w1,q1)...d phi_an(wn,qn) = (2Pi)^(d+1) delta(sum w) delta(sum q) Gamma(n)_{a1..an}(w1,q1; ...; w(n-1),q(n-1))",
   "  all momenta incoming, delta stripped, LAST LEG = -(sum of the others), symmetrised over identical fields",
   "Gamma2_ab(w,q): leg a at (w,q), leg b at (-w,-q).  Vectors are symbols: q1.q2 = scal[q1,q2], q1^2 = scal[q1,q1]"};

ShowVertex[v_Association, tuple_List] := Module[{legs = legSymbols[Length[tuple] - 1]},
   Print[vertexName[tuple], legString[legs], " = ", displayString[VertexValue[v, tuple, legs]]]];

ShowGamma2[v_Association] := Module[{w = Symbol["FRGscalarMSR`Vertices`Private`\[Omega]"], q = Symbol["FRGscalarMSR`Vertices`Private`q"], m},
   m = Gamma2Matrix[v, w, q];
   Print["Gamma2_ab(", displayString[w], ",q), leg a at (", displayString[w], ",q), leg b at (-", displayString[w], ",-q), regulator not included:"];
   Do[If[m[[a, b]] =!= 0, Print["  Gamma2[", v["Fields"][[a]], ", ", v["Fields"][[b]], "] = ", displayString[Collect[m[[a, b]], w]]]],
    {a, Length[m]}, {b, Length[m]}]];

ShowVertexTable[v_Association] := (
   Print["---- Fourier vertices (", Length[v["VertexTable"]], " nonzero field multisets) ----"];
   Scan[Print, conventionHeader];
   ShowGamma2[v];
   Scan[If[Length[#] > 2, ShowVertex[v, #]] &, Keys[v["VertexTable"]]];
   Print["---- end of vertex table ----"]);

(* ---------- BuildVertices ---------- *)

BuildVertices[ansatz_] := Block[{$FRGModule = "Vertices"},
  Module[{t0 = AbsoluteTime[], theory, fields, phis, bars, sortF, mons, nmax, msets, table, checks = {}, warns = {}, stage, w, q, g2, bad,
    pw, conv, nc, fns, fnPat, toConst},
   validateStage[ansatz, "Ansatz"];
   theory = ansatz["Theory"]; fields = theory["Fields"]; phis = theory["Scalars"]; bars = theory["ResponseFields"];
   fns = Lookup[ansatz, "MomentumFunctions", {}]; fnPat = Alternatives @@ fns;
   (* a momentum function F of a quadratic monomial is F(p^2) with p the momentum of a leg (both legs carry +-p): F[scal[qv[1], qv[1]]] *)
   mons = If[#FieldCount == 2, MapAt[# /. f : fnPat :> f[scal[qv[1], qv[1]]] &, #, "Coefficient"], #] & /@ ansatz["Monomials"];
   sortF[l_] := SortBy[l, First[FirstPosition[fields, #]] &];
   nmax = Max[4, Max[#FieldCount & /@ mons]];
   msets = Join @@ Table[Union[sortF /@ Tuples[fields, n]], {n, 2, nmax}];
   table = Association[DeleteCases[(# -> vertexExpr[#, mons, sortF]) & /@ msets, _ -> 0]];

   (* (a) plane waves into the raw density, every multiset (also those with no vertex: nothing may be lost). A momentum function is
      replaced by a constant on both sides: the check covers everything but its argument, which is a definition (F(-Laplacian)) *)
   nc = If[theory["Dimension"] === 2, 2, 3];
   toConst = Thread[fns -> (pwConst /@ fns)];
   bad = Select[msets, Expand[planeWaveCoefficient[ansatz["Density"] /. toConst, fields, #, nc] -
         symbolAsPlaneWave[Lookup[table, Key[#], 0] /. f : fnPat[_] :> pwConst[Head[f]], nc]] =!= 0 &];
   AppendTo[checks, assertCheck["Vertices.planeWave", bad === {}, bad,
     "the Fourier vertex differs from the plane-wave substitution into the density: parser or Fourier-transform error",
     "plane waves in the raw density = derived vertex, " <> ToString[Length[msets]] <> " field multisets (n = 2.." <> ToString[nmax] <> ", d = " <>
      ToString[nc] <> " components" <> If[fns === {}, "", "; momentum functions as constants, their argument p^2 by definition"] <> ")"]];

   stage = <|"Fields" -> fields, "VertexTable" -> table|>;
   g2 = Gamma2Matrix[stage, w, q];
   (* implementation checks (E11): near-tautological by construction, kept as guards *)
   AppendTo[checks, assertCheck["Vertices.gamma2Symmetry", Expand[g2 - Transpose[Gamma2Matrix[stage, -w, -q]]] === ConstantArray[0, Dimensions[g2]],
     g2, "internal: Gamma2_ab(w,q) must equal Gamma2_ba(-w,-q)", "(implementation check) Gamma2_ab(w,q) = Gamma2_ba(-w,-q)"]];
   AppendTo[checks, assertCheck["Vertices.causality", AllTrue[Tuples[phis, 2], VertexValue[stage, #, {{w, q}}] === 0 &],
     Select[Tuples[phis, 2], VertexValue[stage, #, {{w, q}}] =!= 0 &], "a term quadratic in physical fields without response field (MSR)",
     "(implementation check) Gamma2_{phi phi} = 0 for physical fields"]];
   AppendTo[checks, assertCheck["Vertices.gamma2Polynomial", AllTrue[Flatten[g2], PolynomialQ[# /. scal[q, q] -> ss /. fnPat[ss] -> fnAtom, {w, ss}] && FreeQ[# /. scal[q, q] -> ss, q] &],
     g2, "Gamma2 entries must be polynomials in I w and q^2", "(implementation check) Gamma2 entries polynomial in I w and q^2" <>
      If[fns === {}, "", " (momentum functions F(q^2) as coefficients)"]]];
   (* design/MOMENTUMFUNCTIONS.md A3: in every Gamma2 slot (entry, power of I w) a function enters as one term
      number x (at most one coupling) x F(q^2) x (q^2)^n0 *)
   If[fns =!= {}, Module[{slots, badT},
     slots = Flatten[Table[{fields[[i]], fields[[j]], m, Expand[Coefficient[g2[[i, j]] /. scal[q, q] -> ss, w, m]/I^m]},
        {i, Length[fields]}, {j, i, Length[fields]}, {m, 0, Max[0, Exponent[g2[[i, j]], w]]}], 2];
     slots = Select[slots, ! FreeQ[#[[4]], fnPat] &];
     badT = Select[slots, Function[sl, Module[{fa = Union[Cases[sl[[4]], fnPat[_], {0, Infinity}]]},
          ! AllTrue[fa, Function[f, Module[{t = D[sl[[4]] /. f -> fz, fz]},
               FreeQ[t, fnPat | fz] && Head[Expand[t]] =!= Plus && With[{tt = Expand[t] /. ss -> 1},
                 NumericQ[tt] || MatchQ[tt, _?NumericQ _Symbol | _Symbol]]]]]]]];
     AppendTo[checks, assertCheck["Vertices.functionalTerm", badT === {}, {#[[1 ;; 2]], "(I w)^" <> ToString[#[[3]]], displayString[#[[4]]]} & /@ badT,
       "one term ref x Fnu x (p^2)^n0 per function: absorb any extra factor (another coupling, a power, an inverse, a second power of p^2) into the function",
       "momentum functions enter each Gamma2 slot as number x (at most one coupling) x F(q^2) x (q^2)^n0"]]]];

   (* convention check (E1): the free response G_{phibar phi}(w) = 1 / Gamma2_{phi phibar}(w) has its pole in the upper half plane
      (with the couplings at random positive values, q^2 = 1/2) *)
   conv = Table[Module[{e = VertexValue[stage, {phis[[a]], bars[[a]]}, {{w, q}}], cs, roots},
       cs = Complement[Cases[e, s_Symbol /; Context[s] =!= "System`", Infinity], {w}];
       e = e /. scal[q, q] -> 1/2 /. fnPat[_] -> 1 /. RandomRationalPoint[cs, 1];
       roots = If[Exponent[e, w] >= 1, w /. Solve[e == 0, w], {}];
       Which[
        roots === {}, infoCheck["WARN", "Vertices.retardedPole[" <> ToString[phis[[a]]] <> "]", "no dt term in Gamma2_{phi phibar}: convention check skipped"],
        AnyTrue[roots, Im[#] < 0 &], assertCheck["Vertices.retardedPole[" <> ToString[phis[[a]]] <> "]", False, roots,
         "a free response pole is in the lower half plane: sign of the damping term, or of the Fourier convention"],
        AnyTrue[roots, Im[#] == 0 &], infoCheck["WARN", "Vertices.retardedPole[" <> ToString[phis[[a]]] <> "]",
         "pole on the real axis at " <> inputForm[roots] <> " (no damping): check inconclusive"],
        True, assertCheck["Vertices.retardedPole[" <> ToString[phis[[a]]] <> "]", True, roots, "",
         "pole of G_{" <> ToString[bars[[a]]] <> " " <> ToString[phis[[a]]] <> "} at " <> inputForm[roots] <> " (upper half plane)"]]],
     {a, Length[phis]}];
   checks = Join[checks, Select[conv, #Status === "PASS" &]]; warns = Join[warns, Select[conv, #Status =!= "PASS" &]];
   If[AnyTrue[Keys[table], Length[#] == 4 &], AppendTo[warns, infoCheck["WARN", "Vertices.gamma4", "Gamma(4) != 0: the tadpole term of the flow is hand-derived, untested"]]];

   stage = stampStage[<|"Stage" -> "Vertices", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {ansatz["Fingerprint"]},
      "Ansatz" -> ansatz, "Fields" -> fields, "VertexTable" -> table, "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   ShowVertexTable[stage];
   PrintCard[stage];
   stage]];

End[];
EndPackage[];
