(* ::Package:: *)

(* OldProjection.wl -- the p^2 projection operator used before design/MOMENTUMPROJECTION.md (the script's Laplacian trick, exact for
   n <= 1 only), kept test-only as the reference of the Invariant: the new operator must give Expand-identical integrands for every
   n <= 1 monomial. oldLapProject[e, m, n, {nu, p}, {w, q}, dim, deg, eps] = the former private `project`. *)

Begin["FRGscalarMSR`Private`"];

oldLapDimFactors[dim_, q_, lam_] := {Derivative[0, 1][scal][a_, b_] :> (a /. q -> lam), Derivative[1, 0][scal][a_, b_] :> (b /. q -> lam),
   Derivative[0, 2][scal][a_, b_] -> 0, Derivative[2, 0][scal][a_, b_] -> 0, Derivative[1, 1][scal][a_, b_] -> dim};

oldLapDropWedge[e0_, {nu_, p_}, {w_, q_}] := Module[{e = FRGscalarMSR`Internal`reduceWedge2D[e0], ew, wz},
   ew = e /. {wedge[p, q] -> wz, wedge[q, p] -> -wz};
   ew /. wz -> 0];

oldLapProject[e0_, m_, n_, {nu_, p_}, {w_, q_}, dim_, deg_, eps_ : False] := Module[{e = e0, lam},
   If[n > 1, Return[$Failed]];
   If[eps, e = oldLapDropWedge[e, {nu, p}, {w, q}]];
   If[m > 0, e = D[e, {nu, m}]];
   If[n == 1, e = D[e, p, p] /. oldLapDimFactors[dim, q, lam]];
   e = e /. p -> 0 /. nu -> 0;
   If[n == 1, e = (e /. lam -> 0) + (D[e, {lam, 2}] /. lam -> 0)/2 scal[q, q]];
   e = e /. scal[q, q] -> LoopMomentumSquared /. deltaCompletion[deg];
   e/(I^m m!)/If[n == 1, 2 dim, 1]];

(* the records of a Projection stage whose integrand differs (after Expand) from the old operator's; {} = the Invariant holds *)
oldLapMismatches[pr_Association] := Module[{loops = pr["Loops"], sy = pr["Symmetries"], v, dim, deg, wq, eps},
   v = loops["Vertices"]; dim = sy["Dimension"]; deg = loops["Propagator"]["DeltaDegree"]; wq = LoopVariables[loops];
   eps = Lookup[v["Ansatz"], "Epsilon", False];
   Select[pr["Integrands"], Expand[oldLapProject[LoopIntegrand[loops, #Pair, oldNu, oldP], #Frequency, #Momentum, {oldNu, oldP}, wq, dim, deg, eps] -
        (#Integrand)] =!= 0 &]];

End[];
