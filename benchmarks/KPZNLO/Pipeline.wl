(* ::Package:: *)

(* KPZ at NLO (design/NLOPIPELINE.md): the stages of benchmarks/KPZLO (f_lambda = 1) and benchmarks/KPZLOLambda (flowing f_lambda),
   with the external frequency w kept in the loop: BuildLOEquations[..., "ExternalFrequency" -> True]. Functions of (p^2, w) on the
   left-hand side; every function in the loop at w = 0 (Kloss, Canet, Wschebor 2012, Eq. 43); w enters through G(omega + w, Q).
   kpznlo*: f_lambda = 1 (KCW "LO" = BMWflow "NLO"); kpzlamnlo*: flowing f_lambda (KCW NLO). Built once per kernel.
   Needs FRGscalarMSR and the two LO pipelines (benchmarks/KPZLO/Pipeline.wl, benchmarks/KPZLOLambda/Pipeline.wl) loaded. *)

If[! ValueQ[kpzlamnloAdim],
  kpznloEquations = BuildLOEquations[kpzloLoops, kpzloSymmetries, "ExternalFrequency" -> True];
  kpznloAdim = AdimensionalEquations[kpznloEquations, kpzloDims];
  kpzlamnloEquations = BuildLOEquations[kpzlamLoops, kpzlamSymmetries, "ExternalFrequency" -> True];
  kpzlamnloAdim = AdimensionalEquations[kpzlamnloEquations, kpzlamDims]];
