(* ::Package:: *)

(* Projection.wl -- ProjectOnMonomials: for every classified monomial (I nu)^m (p^2)^n of d_s Gamma2_ab(nu, p), apply
     (d/d(I nu))^m / m!  and, for n = 1, the p-Laplacian / (2 d), at nu = 0, p = 0 (pipeline 1: external variables set to zero),
   to the loop integrand; the same operator applied to Gamma_k (regulator off) is asserted to give its coefficient (LHS check).
   Objective, assumptions and checks: docs/MODULES.md section "Projection". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

ProjectOnMonomials::usage = "ProjectOnMonomials[loops, symmetries] returns the \"Projection\" stage: one w-integrand per classified \
monomial (all classes except \"Structure\"), at nu = p = 0, as a function of w and q^2 = LoopMomentumSquared.";
LoopMomentumSquared::usage = "LoopMomentumSquared is the symbol for q^2 in projected integrands and in the flow equations.";
ShowIntegrand::usage = "ShowIntegrand[projected, coefficient] or ShowIntegrand[projected, {a, b}, {m, n}] prints the w-integrand of a monomial.";

Begin["`Private`"];

(* p-Laplacian with vectors as bare symbols (DE4plus_fast.wl section 3): D[D[f, p], p] /. dimFactors at p = 0 is the Laplacian;
   a bare vector q left by a first derivative of scal[p, q] is marked by lam, so that q.q = lam^2 -> q^2 *)
dimFactors[dim_, q_, lam_] := {Derivative[0, 1][scal][a_, b_] :> (a /. q -> lam), Derivative[1, 0][scal][a_, b_] :> (b /. q -> lam),
   Derivative[0, 2][scal][a_, b_] -> 0, Derivative[2, 0][scal][a_, b_] -> 0, Derivative[1, 1][scal][a_, b_] -> dim};

(* derivatives of Delta(w, s) = sum_m c_m(s) (I w)^m *)
deltaCompletion[deg_] := Derivative[i_, j_][\[CapitalDelta]][w_, s_] :>
   (D[Sum[DeltaCoefficient[mm][xs] (I yw)^mm, {mm, 0, deg}], {yw, i}, {xs, j}] /. {yw -> w, xs -> s});

(* with eps (d = 2): products of wedges are reduced by F2; every term is then wedge[p, q]^e f(p^2, p.q, q^2, nu) with e = 0 or 1, and
   the e = 1 part has a vanishing (I nu)^m (p^2)^n projection (n <= 1, F4 of design/LEVICIVITA.md): it is dropped, after checking
   that its projection vanishes on explicit 2-vectors. $wedgeDropped counts the integrands with a nonzero dropped part. *)
dropWedge[e0_, m_, n_, {nu_, p_}, {w_, q_}] := Module[{e = FRGscalarMSR`Internal`reduceWedge2D[e0], ew, odd, ex, pc, qc, vec},
   ew = e /. {wedge[p, q] -> wz, wedge[q, p] -> -wz};
   FRGAssert["Projection.wedgeOdd", FreeQ[ew, wedge] && Exponent[ew, wz] <= 1 && FreeQ[ew, Power[wz, _?Negative]],
    Cases[ew, _wedge | Power[wz, _], {0, Infinity}], "internal: a wedge other than +-wedge[p, q] to the power <= 1 is left after the 2D reduction"];
   odd = Coefficient[ew, wz, 1];
   If[odd =!= 0,
    vec = {p -> {pc[1], pc[2]}, q -> {qc[1], qc[2]}};
    ex = (pc[1] qc[2] - pc[2] qc[1]) odd /. scal[a_, b_] :> Expand[(a /. vec) . (b /. vec)];
    If[m > 0, ex = D[ex, {nu, m}]];
    If[n == 1, ex = D[ex, {pc[1], 2}] + D[ex, {pc[2], 2}]];
    ex = ex /. {pc[_] -> 0, nu -> 0};
    FRGAssert["Projection.wedgeOdd", ExactlyZeroQ[ex, "Seeds" -> {1}], Short[ex],
     "internal: the part odd in wedge[p, q] has a nonzero projection (explicit 2-vectors)"];
    $wedgeDropped++];
   ew /. wz -> 0];

project[e0_, m_, n_, {nu_, p_}, {w_, q_}, dim_, deg_, eps_ : False] := Module[{e = e0, lin},
   If[eps, e = dropWedge[e, m, n, {nu, p}, {w, q}]];
   If[m > 0, e = D[e, {nu, m}]];
   If[n == 1, e = D[e, p, p] /. dimFactors[dim, q, lam]];
   e = e /. p -> 0 /. nu -> 0;
   If[n == 1,
    lin = D[e, lam] /. lam -> 0;
    FRGAssert["Projection.evenInQ", ExactlyZeroQ[lin, "Seeds" -> {1}], Short[lin], "internal: a term odd in the loop momentum survived the p-Laplacian"];
    e = (e /. lam -> 0) + (D[e, {lam, 2}] /. lam -> 0)/2 scal[q, q]];
   e = e /. scal[q, q] -> LoopMomentumSquared /. deltaCompletion[deg];
   FRGAssert["Projection.complete", FreeQ[e, q | p | nu | lam | wedge | Derivative[__][scal] | Derivative[__][\[CapitalDelta]]], e,
    "internal: external variables, vectors or Delta derivatives left after the projection"];
   e/(I^m m!)/If[n == 1, 2 dim, 1]];

ProjectOnMonomials[loops_, symmetries_] := Block[{$FRGModule = "Projection"},
  Module[{t0 = AbsoluteTime[], v, prop, dim, deg, targets, lhsBad, byPair, res, checks = {}, warns = {}, stage, wq, eps},
   validateStage[loops, "Loops"]; validateStage[symmetries, "Symmetries"];
   v = loops["Vertices"]; prop = loops["Propagator"];
   FRGAssert["Projection.sameAnsatz", v["Upstream"] === symmetries["Upstream"], {v["Upstream"], symmetries["Upstream"]},
    "the loops and the symmetries were built from different ansatz stages"];
   dim = symmetries["Dimension"]; deg = prop["DeltaDegree"]; wq = LoopVariables[loops];
   eps = Lookup[v["Ansatz"], "Epsilon", False];
   If[eps, AppendTo[checks, assertCheck["Projection.epsilonDimension", dim === 2, dim, "internal: eps (jac) needs the declared d = 2", "eps with d = 2"]]];
   $wedgeDropped = 0;
   targets = Select[symmetries["Classification"], #Class =!= "Structure" &];

   (* LHS check: the projection operator applied to Gamma_k (regulator off) reproduces the classified coefficient *)
   lhsBad = Select[targets, Expand[project[VertexValue[v, #Pair, {{nu, p}}], #Frequency, #Momentum, {nu, p}, {zz, zq}, dim, deg, eps] - #GammaK] =!= 0 &];
   AppendTo[checks, assertCheck["Projection.LHS", lhsBad === {}, lhsBad,
     "internal: the projection operator does not reproduce the Gamma_k coefficient (normalisation of the frequency or momentum projection)",
     "projection of Gamma_k = its coefficient, " <> ToString[Length[targets]] <> " monomials ((d/d(I nu))^m/m!, Laplacian_p/(2 " <> ToString[dim] <> "))"]];

   byPair = GroupBy[targets, #Pair &];
   res = Join @@ KeyValueMap[Function[{pair, ms}, Module[{loop = LoopIntegrand[loops, pair, nu, p]},
        Append[#, "Integrand" -> project[loop, #Frequency, #Momentum, {nu, p}, wq, dim, deg, eps]] & /@ ms]], byPair];
   AppendTo[checks, assertCheck["Projection.integrandsBuilt", Length[res] == Length[targets], Length[res], "internal",
     ToString[Length[res]] <> " integrands at nu = p = 0 (largest LeafCount " <> ToString[Max[LeafCount[#Integrand] & /@ res]] <> ")"]];
   If[eps, AppendTo[checks, assertCheck["Projection.wedgeOdd", True, Null, "",
      "2D reduction (F2): one wedge at most; the part odd in wedge[p, q] dropped (F4) in " <> ToString[$wedgeDropped] <>
       " integrand(s), its projection checked to vanish on explicit 2-vectors"]]];

   stage = stampStage[<|"Stage" -> "Projection", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {loops["Fingerprint"], symmetries["Fingerprint"]},
      "Loops" -> loops, "Symmetries" -> symmetries, "Integrands" -> res, "LoopFrequency" -> First[wq],
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   PrintCard[stage];
   stage]];

ShowIntegrand[pr_Association, {a_, b_}, {m_, n_}] := Module[{r = SelectFirst[pr["Integrands"], #Pair === {a, b} && #Frequency == m && #Momentum == n &]},
   If[MissingQ[r], Print["no such monomial (class Structure, or outside the truncation)"],
    Print["w-integrand of (I nu)^", m, " (p^2)^", n, " in d_s Gamma2[", a, ", ", b, "] (", r["Class"], ", LeafCount ", LeafCount[r["Integrand"]], "):"];
    Print[displayString[r["Integrand"]]]]];
ShowIntegrand[pr_Association, coef_] := Module[{r = SelectFirst[pr["Integrands"], Expand[#GammaK - coef] === 0 &]},
   If[MissingQ[r], Print["no monomial with Gamma_k coefficient ", coef], ShowIntegrand[pr, r["Pair"], {r["Frequency"], r["Momentum"]}]]];

End[];
EndPackage[];
