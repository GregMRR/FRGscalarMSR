(* ::Package:: *)

(* Suite hook of benchmarks/KuramotoSivashinsky (Gosteva, Wschebor, Canet, arXiv:2607.15784): the library's KS flows equal the
   paper's Eqs. (34)-(37) exactly, including the p^4 projection (d_s tau), plus the time-reversal identities at d = 1. The benchmark's
   logic lives in benchmarks/KuramotoSivashinsky; this file only wraps its checks in test[...]. *)

ksBenchDir = FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KuramotoSivashinsky"}];
quiet[Get[FileNameJoin[{ksBenchDir, "Pipeline.wl"}]]];
Get[FileNameJoin[{ksBenchDir, #}]] & /@ {"Reference.wl", "Checks.wl"};

Scan[test["KS benchmark: " <> #Name <> " [" <> #Paper <> "]", #Pass, True] &, ksChecks[]];

(* the comparison discriminates: a 1/100 change of one coefficient of Eq. (37) (3 hn^4 sd -> 303/100 hn^4 sd) is detected *)
test["KS benchmark: a perturbed Eq. (37) is rejected",
  ExactlyZeroQ[ksLibrary["tau"] - ((ksRef["tau"] + Gh/(12 (Fh + Th qh^2 + rnu[qh^2])^2 qh^4) (3/100) sD[qh^2]) /. ksSRules)], False];
