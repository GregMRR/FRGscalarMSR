(* ::Package:: *)

(* Projection.wl -- ProjectOnMonomials: for every classified monomial (I nu)^m (p^2)^n of d_s Gamma2_ab(nu, p), apply
     (d/d(I nu))^m / m!  and the (p^2)^n coefficient (direction average of the s^(2n) Taylor coefficient of p = s e, any n;
   design/MOMENTUMPROJECTION.md M1), at nu = 0, p = 0 (pipeline 1: external variables set to zero), to the loop integrand; the same operator applied to Gamma_k (regulator off) is asserted to give its coefficient (LHS check).
   Objective, assumptions and checks: docs/MODULES.md section "Projection". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

ProjectOnMonomials::usage = "ProjectOnMonomials[loops, symmetries] returns the \"Projection\" stage: one w-integrand per classified \
monomial (all classes except \"Structure\"), at nu = p = 0, as a function of w and q^2 = LoopMomentumSquared.";
LoopMomentumSquared::usage = "LoopMomentumSquared is the symbol for q^2 in projected integrands and in the flow equations.";
ExternalMomentumSquared::usage = "ExternalMomentumSquared is the symbol for p^2 in the integrands of the flows at generic p (functional slots, \
\"AllMomenta\" tags); MomentumProduct is p.q there (design/MOMENTUMFUNCTIONS.md C).";
MomentumProduct::usage = "MomentumProduct is the symbol for p.q (midpoint loop momentum q) in the integrands of the flows at generic p.";
ShowIntegrand::usage = "ShowIntegrand[projected, coefficient] or ShowIntegrand[projected, {a, b}, {m, n}] prints the w-integrand of a monomial.";

Begin["`Private`"];

(* derivatives of Delta(w, s) = sum_m c_m(s) (I w)^m *)
deltaCompletion[deg_] := Derivative[i_, j_][\[CapitalDelta]][w_, s_] :>
   (D[Sum[DeltaCoefficient[mm][xs] (I yw)^mm, {mm, 0, deg}], {yw, i}, {xs, j}] /. {yw -> w, xs -> s});

(* <c^k> |q|^k, c = cosine between a random unit direction e and q in d dimensions (design/MOMENTUMPROJECTION.md M1) *)
projCAverage[k_, dim_, qq_] := If[OddQ[k], 0, qq^(k/2) (k - 1)!!/Product[dim + 2 j, {j, 0, k/2 - 1}]];

(* M1: the (p^2)^n coefficient of f(p) is the direction average of the s^(2n) Taylor coefficient of f(p = s e). The loop holds p only
   in scal[p, p] = s^2 and scal[p, q] = s c (scal is bilinear: (p + q)^2 is already expanded); c^k averages to projCAverage *)
directionalAverage[e0_, n_, p_, q_, dim_] := Module[{e = e0 /. {scal[p, p] -> projS^2, scal[p, q] -> projS projC}, odd},
   If[n == 0, Return[e /. projS -> 0]];
   e = (D[e, {projS, 2 n}] /. projS -> 0)/(2 n)!;
   odd = Table[D[e, {projC, k}] /. projC -> 0, {k, 1, 2 n - 1, 2}];
   FRGAssert["Projection.noOddC", ExactlyZeroQ[odd, "Seeds" -> {1}], Short[odd],
    "internal invariant: a term odd in p.q survived at order s^(2n) (impossible by degree counting: c-degree is 2n - 2a)"];
   Sum[(D[e, {projC, k}] /. projC -> 0)/k! projCAverage[k, dim, scal[q, q]], {k, 0, 2 n, 2}]];

(* N_n = (Laplacian_p)^n (p^2)^n = 4^n n! (d/2)_n  (M2) *)
lapNorm[n_, d_] := 4^n n! Pochhammer[d/2, n];

(* the same projection on explicit d-vectors: (Laplacian_p)^n at p = 0 / N_n; returns {coefficient, q.q}. Independent of
   directionalAverage and projCAverage (used by Projection.explicitVectors and the tests) *)
explicitProjection[e0_, m_, n_, {nu_, p_}, {w_, q_}, d_Integer, deg_] := Module[{pc = Array[projPc, d], qc = Array[projQc, d], vec, e},
   vec = {p -> pc, q -> qc};
   e = e0 /. {scal[a_, b_] :> Expand[(a /. vec) . (b /. vec)], wedge[a_, b_] :> With[{x = a /. vec, y = b /. vec}, x[[1]] y[[2]] - x[[2]] y[[1]]]};
   If[m > 0, e = D[e, {nu, m}]];
   e = Nest[Function[g, Sum[D[g, {x, 2}], {x, pc}]], e, n] /. Thread[pc -> 0] /. nu -> 0;
   {(e /. deltaCompletion[deg])/(I^m m! lapNorm[n, d]), qc . qc}];

(* Projection.explicitVectors (n >= 2): the symbolic projection equals the explicit one at the declared d, or at d = 1..n+1 for a
   symbolic d (the p^(2n) result spans n+1 independent d-structures, review E4) *)
checkExplicitVectors[loop_, integrand_, m_, n_, {nu_, p_}, wq_, dim_, deg_] := Module[{ds = If[IntegerQ[dim], {dim}, Range[n + 1]], bad},
   bad = Select[ds, Function[d, With[{x = explicitProjection[loop, m, n, {nu, p}, wq, d, deg]},
        ! ExactlyZeroQ[x[[1]] - (integrand /. LoopMomentumSquared -> x[[2]] /. If[IntegerQ[dim], {}, dim -> d]), "Seeds" -> {1}]]]];
   FRGAssert["Projection.explicitVectors", bad === {}, bad,
    "internal: the directional average differs from (Laplacian_p)^n / N_n on explicit d-vectors in these dimensions"]];

(* Projection.normalisation: the full projection of (p.q)^(2k) (p^2)^(n-k) against <c^(2k)> q^(2k) = q^(2k) Gamma(k+1/2) Gamma(d/2) /
   (Sqrt[Pi] Gamma(k+d/2)), written independently of projCAverage (review E17) *)
checkNormalisation[nmax_, dim_] := Module[{bad},
   bad = Select[Flatten[Table[{n, k}, {n, 1, nmax}, {k, 0, n}], 1], Function[nk, With[{n = nk[[1]], k = nk[[2]]},
        Together[FunctionExpand[project[scal[nP, nQ]^(2 k) scal[nP, nP]^(n - k), 0, n, {nNu, nP}, {nW, nQ}, dim, 0] -
            LoopMomentumSquared^k Gamma[k + 1/2] Gamma[dim/2]/(Sqrt[Pi] Gamma[k + dim/2])]] =!= 0]]];
   FRGAssert["Projection.normalisation", bad === {}, bad,
    "internal: the projection of (p.q)^(2k) (p^2)^(n-k) differs from <c^(2k)> q^(2k) (these {n, k})"]];

(* with eps (d = 2): products of wedges are reduced by F2; every term is then wedge[p, q]^e f(p^2, p.q, q^2, nu) with e = 0 or 1, and
   the e = 1 part has a vanishing (I nu)^m (p^2)^n projection for every n (M3; F4 of design/LEVICIVITA.md): it is dropped, after
   checking that its projection vanishes on explicit 2-vectors. $wedgeDropped counts the integrands with a nonzero dropped part. *)
dropWedge[e0_, m_, n_, {nu_, p_}, {w_, q_}] := Module[{e = FRGscalarMSR`Internal`reduceWedge2D[e0], ew, odd, ex, pc, qc, vec},
   ew = e /. {wedge[p, q] -> wz, wedge[q, p] -> -wz};
   FRGAssert["Projection.wedgeOdd", FreeQ[ew, wedge] && Exponent[ew, wz] <= 1 && FreeQ[ew, Power[wz, _?Negative]],
    Cases[ew, _wedge | Power[wz, _], {0, Infinity}], "internal: a wedge other than +-wedge[p, q] to the power <= 1 is left after the 2D reduction"];
   odd = Coefficient[ew, wz, 1];
   If[odd =!= 0,
    vec = {p -> {pc[1], pc[2]}, q -> {qc[1], qc[2]}};
    ex = (pc[1] qc[2] - pc[2] qc[1]) odd /. scal[a_, b_] :> Expand[(a /. vec) . (b /. vec)];
    If[m > 0, ex = D[ex, {nu, m}]];
    ex = Nest[D[#, {pc[1], 2}] + D[#, {pc[2], 2}] &, ex, n];
    ex = ex /. {pc[_] -> 0, nu -> 0};
    FRGAssert["Projection.wedgeOdd", ExactlyZeroQ[ex, "Seeds" -> {1}], Short[ex],
     "internal: the part odd in wedge[p, q] has a nonzero projection (explicit 2-vectors)"];
    $wedgeDropped++];
   ew /. wz -> 0];

(* (d/d(I nu))^m / m! at nu = 0, and the (p^2)^n coefficient by directionalAverage *)
project[e0_, m_, n_, {nu_, p_}, {w_, q_}, dim_, deg_, eps_ : False] := Module[{e = e0},
   If[eps, e = dropWedge[e, m, n, {nu, p}, {w, q}]];
   If[m > 0, e = D[e, {nu, m}]];
   e = directionalAverage[e, n, p, q, dim] /. nu -> 0;
   e = e /. scal[q, q] -> LoopMomentumSquared /. deltaCompletion[deg];
   FRGAssert["Projection.complete", FreeQ[e, q | p | nu | projS | projC | wedge | Derivative[__][scal] | Derivative[__][\[CapitalDelta]]], e,
    "internal: external variables, vectors or Delta derivatives left after the projection"];
   e/(I^m m!)];

(* ---------- flows at generic p (design/MOMENTUMFUNCTIONS.md C1, C2) ---------- *)

(* vectors -> invariants: q^2, p^2, p.q *)
toInvariants[e_, p_, q_] := e /. {scal[p, p] -> ExternalMomentumSquared, scal[p, q] -> MomentumProduct, scal[q, q] -> LoopMomentumSquared};

(* frequency operator at nu = 0 on a midpoint integrand, in the invariants *)
genericFrequency[e0_, m_, {nu_, p_}, {w_, q_}, deg_] := Module[{e = e0},
   If[m > 0, e = D[e, {nu, m}]];
   e = toInvariants[e /. nu -> 0, p, q] /. deltaCompletion[deg];
   FRGAssert["Projection.complete", FreeQ[e, q | p | nu | wedge | _scal | Derivative[__][\[CapitalDelta]]], e,
    "internal: external variables, vectors or Delta derivatives left in a generic-p integrand"];
   e/(I^m m!)];

(* C2: subtract the orders p^(2j), j < n0, by their un-averaged Taylor coefficients at fixed q (p -> t p) and divide by p^(2 n0).
   Projection.functionalOrder: the odd orders t^(2j+1), j < n0, vanish (the q -> -q symmetrisation), so the result is O(1) as p -> 0 *)
subtractLowerOrders[e_, n0_, what_] := Module[{et, cl, odd},
   If[n0 == 0, Return[e]];
   et = e /. {ExternalMomentumSquared -> fnT^2 ExternalMomentumSquared, MomentumProduct -> fnT MomentumProduct};
   cl = Table[(D[et, {fnT, j}] /. fnT -> 0)/j!, {j, 0, 2 n0 - 1}] /. deltaCompletionAny;
   odd = cl[[2 ;; ;; 2]];
   FRGAssert["Projection.functionalOrder", ExactlyZeroQ[odd, "Seeds" -> {1}], {what, "odd orders" -> Short[odd]},
    "the generic-p integrand is not even in p at fixed q: its quotient by p^(2 n0) would be singular as p -> 0 (midpoint routing and q -> -q symmetrisation)"];
   Together[(e - Total[cl[[1 ;; ;; 2]]])/ExternalMomentumSquared^n0]];

(* derivatives of Delta left by the Taylor expansion, any degree <= 2 (the propagator degree is fixed later by the integrand) *)
deltaCompletionAny = Derivative[i_, j_][\[CapitalDelta]][w_, s_] :>
   (D[Sum[DeltaCoefficient[mm][xs] (I yw)^mm, {mm, 0, 2}], {yw, i}, {xs, j}] /. {yw -> w, xs -> s});

ProjectOnMonomials[loops_, symmetries_] := Block[{$FRGModule = "Projection"},
  Module[{t0 = AbsoluteTime[], v, prop, dim, deg, targets, lhsBad, byPair, res, nmax, nExplicit, checks = {}, warns = {}, stage, wq, eps,
    generic, genRes, lhsGen},
   validateStage[loops, "Loops"]; validateStage[symmetries, "Symmetries"];
   v = loops["Vertices"]; prop = loops["Propagator"];
   FRGAssert["Projection.sameAnsatz", v["Upstream"] === symmetries["Upstream"], {v["Upstream"], symmetries["Upstream"]},
    "the loops and the symmetries were built from different ansatz stages"];
   dim = symmetries["Dimension"]; deg = prop["DeltaDegree"]; wq = LoopVariables[loops];
   eps = Lookup[v["Ansatz"], "Epsilon", False];
   If[eps, AppendTo[checks, assertCheck["Projection.epsilonDimension", dim === 2, dim, "internal: eps (jac) needs the declared d = 2", "eps with d = 2"]]];
   $wedgeDropped = 0;
   targets = Select[symmetries["Classification"], #Class =!= "Structure" && #Class =!= "Functional" && IntegerQ[#Momentum] &];
   generic = Select[symmetries["Classification"], #Class === "Functional" || #Momentum === "AllMomenta" &];
   If[generic =!= {}, FRGAssert["Projection.functionalEpsilon", ! eps, generic[[All, "Pair"]],
     "jac (eps) together with momentum functions or \"AllMomenta\" tags is not supported (docs/FUTURE.md)"]];
   nmax = Max[0, targets[[All, "Momentum"]]];

   (* the d-normalisation of the momentum projection, against an independent closed form (LHS cannot see it: review E2) *)
   If[nmax >= 1, checkNormalisation[nmax, dim];
    AppendTo[checks, assertCheck["Projection.normalisation", True, Null, "",
      "(p.q)^(2k) (p^2)^(n-k) -> <c^(2k)> q^(2k) = q^(2k) Gamma(k+1/2) Gamma(d/2)/(Sqrt[Pi] Gamma(k+d/2)), n <= " <> ToString[nmax] <> ", d = " <> ToString[dim]]]];
   If[nmax >= 2, AppendTo[warns, infoCheck["ASSUMPTION", "Projection.shapeSmoothness",
      "p^{2n} flows contain shape derivatives up to order 2n (here " <> ToString[2 nmax] <> "): any explicit shape substituted later must be C^" <> ToString[2 nmax]]]];

   (* LHS check: the projection operator applied to Gamma_k (regulator off) reproduces the classified coefficient *)
   lhsBad = Select[targets, Expand[project[VertexValue[v, #Pair, {{nu, p}}], #Frequency, #Momentum, {nu, p}, {zz, zq}, dim, deg, eps] - #GammaK] =!= 0 &];
   AppendTo[checks, assertCheck["Projection.LHS", lhsBad === {}, lhsBad,
     "internal: the projection operator does not reproduce the Gamma_k coefficient (normalisation of the frequency or momentum projection)",
     "projection of Gamma_k = its coefficient, " <> ToString[Length[targets]] <> " monomials ((d/d(I nu))^m/m!, (p^2)^n coefficient; not the d-normalisation)"]];

   byPair = GroupBy[targets, #Pair &];
   res = Join @@ KeyValueMap[Function[{pair, ms}, Module[{loop = LoopIntegrand[loops, pair, nu, p]},
        Function[r, With[{ig = project[loop, r["Frequency"], r["Momentum"], {nu, p}, wq, dim, deg, eps]},
           (* n >= 2: the directional average against explicit vectors (n <= 1 is certified by the regression vs the old operator) *)
           If[r["Momentum"] >= 2, checkExplicitVectors[loop, ig, r["Frequency"], r["Momentum"], {nu, p}, wq, dim, deg]];
           Append[r, "Integrand" -> ig]]] /@ ms]], byPair];
   nExplicit = Count[res, r_ /; r["Momentum"] >= 2];
   If[nExplicit > 0, AppendTo[checks, assertCheck["Projection.explicitVectors", True, Null, "",
      ToString[nExplicit] <> " p^(2n >= 4) integrands = (Laplacian_p)^n/N_n on explicit vectors, d = " <>
       If[IntegerQ[dim], ToString[dim], "1, ..., " <> ToString[nmax + 1]]]]];
   (* generic p: functional slots (subtract, divide; plus the exact p -> 0 row of order n0, C4) and "AllMomenta" slots (whole slot) *)
   If[generic =!= {},
    lhsGen = Select[generic, Function[r, With[{g = genericFrequency[VertexValue[v, r["Pair"], {{nu, p}}], r["Frequency"], {nu, p}, {zz, zq}, deg]},
         If[r["Class"] === "Functional",
          ! ExactlyZeroQ[subtractLowerOrders[g, r["Momentum"], r["Pair"]] - (r["GammaK"] /. r["Function"] -> r["Function"][ExternalMomentumSquared])] ||
           Expand[project[VertexValue[v, r["Pair"], {{nu, p}}], r["Frequency"], r["Momentum"], {nu, p}, {zz, zq}, dim, deg] -
              (r["GammaK"] /. r["Function"] -> r["Function"][0])] =!= 0,
          ! ExactlyZeroQ[g - (r["GammaK"] /. pp2 -> ExternalMomentumSquared)]]]]];
    AppendTo[checks, assertCheck["Projection.LHSGeneric", lhsGen === {}, lhsGen,
      "internal: the generic-p operator does not reproduce the Gamma_k slot (functional: (slot - lower orders)/p^(2 n0) = ref F(p^2), and F(0) at p -> 0)",
      "generic-p operator on Gamma_k = its slot, " <> ToString[Length[generic]] <> " slot(s)"]];
    genRes = Function[r, Module[{mid = MidpointLoopIntegrand[loops, r["Pair"], nu, p], g},
        g = genericFrequency[mid, r["Frequency"], {nu, p}, wq, deg];
        If[r["Class"] === "Functional",
         Join[r, <|"Integrand" -> subtractLowerOrders[g, r["Momentum"], r["Pair"]],
           "AtZero" -> project[LoopIntegrand[loops, r["Pair"], nu, p], r["Frequency"], r["Momentum"], {nu, p}, wq, dim, deg]|>],
         Append[r, "Integrand" -> g]]]] /@ generic;
    AppendTo[checks, assertCheck["Projection.functionalOrder", True, Null, "",
      ToString[Length[genRes]] <> " generic-p integrand(s): midpoint routing, q -> -q symmetrised, d_s R on its line; functional ones minus their " <>
       "orders below n0, divided by p^(2 n0), O(1) as p -> 0; exact p -> 0 row of order n0 for each function"]];
    AppendTo[warns, infoCheck["INFO", "Projection.genericMomentum", "generic-p integrands in q^2 = LoopMomentumSquared (midpoint loop momentum), " <>
       "p^2 = ExternalMomentumSquared, p.q = MomentumProduct; Delta at |q -+ p/2|^2 (and q^2 for the subtracted orders)"]];
    res = Join[res, genRes]];
   AppendTo[checks, assertCheck["Projection.integrandsBuilt", Length[res] == Length[targets] + Length[generic], Length[res], "internal",
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
