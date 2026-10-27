(* ::Package:: *)

(* KPZ LO pipeline (functions of p^2) shared by unit and acceptance tests: the stages of benchmarks/KPZLO/Pipeline.wl, built once
   per run (kpzlo* symbols), plus a builder for variants of the LO ansatz used by the expectFail tests. *)

If[! ValueQ[kpzloAdim], quiet[Get[FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KPZLO", "Pipeline.wl"}]]]];

(* loEquationsOf[theory, density, couplings, functions, tags, regulator entries]: an LO pipeline up to BuildLOEquations *)
loEquationsOf[th_, dens_, cps_, funcs_, tags_List, regs_List] := Module[{an, sy, ve, re, pr, lo},
   an = DefineAnsatz[th, "Density" -> dens, "Couplings" -> cps, "Functions" -> funcs];
   sy = DeclareSymmetries[an, Sequence @@ tags];
   ve = BuildVertices[an]; re = DefineRegulator[an, Sequence @@ regs, "Shapes" -> {rn, rd}];
   pr = BuildPropagator[ve, re]; lo = BuildLoopIntegrands[ve, pr, re];
   BuildLOEquations[lo, sy]];
kpzloDensity = hbar (dt[h] - nuk fn[fnu, lap[h]] - lamk/2 grad[h].grad[h]) - Dk hbar fn[fD, hbar];
kpzloTags = {Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift symmetry"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"]};
kpzloRegs = {{hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2]};
