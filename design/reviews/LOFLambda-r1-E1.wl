(* the-enemy, LOFLAMBDA r1, E1: the Ward check at p = 0 is vacuous for the h-hbar entry: the whole loop integrand of
   d_s Gamma2_{h,hbar}(nu, p) vanishes POINTWISE (any nu, any w, q) at p = 0, because the vertex carrying the external h leg is
   lamk f(k_hbar^2) (p . k). So ExternalMomentum.ward cannot detect any error in the d/d(i nu) step, its sign or normalisation (F3).
   The f_lambda version keeps the same structure (vertex lamk flam(k_hbar^2) p.k), so the conclusion carries over.
   Usage: math -script design/reviews/LOFLambda-r1-E1.wl [repo root, default: cwd]  (run on a checkout of 6df4a91) *)
root = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], Directory[]];
Get[FileNameJoin[{root, "FRGscalarMSR.wl"}]];
Block[{Print}, Get[FileNameJoin[{root, "benchmarks", "KPZLO", "Pipeline.wl"}]]];
Do[With[{e = LoopIntegrand[kpzloLoops, pr, nuX, pX]},
   Print[pr, ": loop integrand (generic nu) at p = 0 identically zero? ",
    Together[Expand[e] /. scal[pX, _] -> 0] === 0,
    "   at generic p zero? ", Together[Expand[e]] === 0]], {pr, {{h, hbar}, {hbar, h}}}];
