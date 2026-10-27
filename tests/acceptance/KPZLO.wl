(* ::Package:: *)

(* Suite hook of benchmarks/KPZLO (BMWflow notes/LO_flow_equations.tex, kpz_notation.tex): the library's KPZ LO stage (loop kernels
   of both functions, eta system, linear part, coupling flow) equals the hand-coded BMWflow LO equations exactly at symbolic d. The
   benchmark's logic lives in benchmarks/KPZLO; this file only wraps its checks in test[...]. *)

kpzloBenchDir = FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZLO"}];
Get[FileNameJoin[{testRoot, "reference", "KPZLOSetup.wl"}]];
Get[FileNameJoin[{kpzloBenchDir, #}]] & /@ {"Reference.wl", "Checks.wl"};

Scan[test["KPZ LO benchmark: " <> #Name <> " [" <> #Paper <> "]", #Pass, True] &, kpzloChecks[]];

(* the comparison discriminates: a 1/100 change of the A_d kernel of I_nu (A_d -> 101/100 A_d) is detected *)
test["KPZ LO benchmark: a perturbed I_nu is rejected",
  Module[{q2 = qh^2, Q2 = ph^2 + qh^2 + 2 ph qh ch, lq, lQ, kQ},
   lq = q2 (fnu[q2] + r[q2]); lQ = Q2 (fnu[Q2] + r[Q2]); kQ = fD[Q2] + r[Q2];
   ExactlyZeroQ[kpzloLibrary["Loop"][[2]] - (kpzloRefLoop["nu"] - 2 g/ph^2 (q2 + ph qh ch) (ph^2 + ph qh ch) kQ (1/100)/(4 lq (lq + lQ)^2) q2 (-etanu r[q2] - 2 q2 r'[q2]))]],
  False];
