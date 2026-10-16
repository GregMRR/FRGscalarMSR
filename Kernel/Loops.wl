(* ::Package:: *)

(* Loops.wl -- BuildLoopIntegrands: the one-loop integrand of d_s Gamma2_ab(nu, p) at fixed loop frequency and momentum (w, q),
     d_s Gamma2_ab(nu, p) = Int dw/(2 Pi) d^dq/(2 Pi)^d { Tr[d_sR G(w,q) Gamma3_a(nu,p;w,q) G(w+nu,p+q) Gamma3_b(-nu,-p;w+nu,p+q) G(w,q)]
                                                       - 1/2 Tr[d_sR G(w,q) Gamma4_ab(nu,p;w,q) G(w,q)] },
   with the 1/(2 Pi) of the w-measure included and Int d^dq/(2 Pi)^d left to the numerics.
   The formula is hand-derived in docs/MODULES.md section "Loops" (E3/E4: not independently tested). *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

BuildLoopIntegrands::usage = "BuildLoopIntegrands[vertices, prop, regulator] returns the \"Loops\" stage; LoopIntegrand[loops, {a, b}, nu, p] \
is the integrand of d_s Gamma2_ab(nu, p) (leg a at (nu, p), leg b at (-nu, -p)) at loop point (w, q) = LoopVariables[loops].";
LoopIntegrand::usage = "LoopIntegrand[loops, {a, b}, nu, p] is the integrand of d_s Gamma2_ab(nu, p), with abstract d_s R entries dsR[c, d].";
LoopVariables::usage = "LoopVariables[loops] = {w, q}: the loop frequency and the loop momentum (a vector symbol) used by LoopIntegrand.";

Begin["`Private`"];

LoopVariables[l_Association] := {wL, qL};

gMat[prop_, w_, vec_] := prop["G"][w, scal[vec, vec]];

gamma4Matrix[v_, a_, b_, {nu_, p_}, {w_, q_}] := Table[VertexValue[v, {a, b, c, d}, {{nu, p}, {-nu, -p}, {w, q}}], {c, v["Fields"]}, {d, v["Fields"]}];

LoopIntegrand[l_Association, {a_, b_}, nu_, p_] := Module[{v = l["Vertices"], prop = l["Propagator"], dsr = l["DsR"], g0, gm, bubble, tadpole},
   g0 = gMat[prop, wL, qL];
   gm = gMat[prop, wL + nu, p + qL];
   bubble = Tr[dsr . g0 . Gamma3Matrix[v, a, {nu, p}, {wL, qL}] . gm . Gamma3Matrix[v, b, {-nu, -p}, {wL + nu, p + qL}] . g0];
   tadpole = If[l["HasGamma4"], Tr[dsr . g0 . gamma4Matrix[v, a, b, {nu, p}, {wL, qL}] . g0], 0];
   (bubble - tadpole/2)/(2 Pi)];

BuildLoopIntegrands[vertices_, prop_, regulator_] := Block[{$FRGModule = "Loops"},
  Module[{t0 = AbsoluteTime[], checks = {}, warns = {}, dsr, g4, stage},
   validateStage[vertices, "Vertices"]; validateStage[prop, "Propagator"]; validateStage[regulator, "Regulator"];
   FRGAssert["Loops.sameUpstream", prop["Upstream"] === {vertices["Fingerprint"], regulator["Fingerprint"]}, prop["Upstream"],
    "the propagator was not built from these vertices and this regulator"];
   dsr = DsRAbstractMatrix[regulator];
   AppendTo[checks, assertCheck["Loops.dsRSymmetric", dsr === Transpose[dsr], dsr, "internal: d_s R must be symmetric", "d_s R symmetric, abstract entries dsR[a, b]"]];
   g4 = AnyTrue[Keys[vertices["VertexTable"]], Length[#] == 4 &];
   If[g4, AppendTo[warns, infoCheck["WARN", "Loops.tadpole", "Gamma4 != 0: tadpole term -1/2 Tr[d_sR G Gamma4 G] included; hand-derived, untested"]]];
   AppendTo[warns, infoCheck["ASSUMPTION", "Loops.normalisation",
     "loop prefactor 1, tadpole -1/2 and 1/(2 Pi) hand-derived (docs/MODULES.md); the oracle and DE4plus_fast.wl share this formula, so neither verifies it"]];
   stage = stampStage[<|"Stage" -> "Loops", "Pipeline" -> "DerivativeExpansion",
      "Upstream" -> {vertices["Fingerprint"], prop["Fingerprint"], regulator["Fingerprint"]},
      "Vertices" -> vertices, "Propagator" -> prop, "Regulator" -> regulator, "Fields" -> vertices["Fields"], "DsR" -> dsr, "HasGamma4" -> g4,
      "Formula" -> "d_s Gamma2_ab(nu,p) = Int dw/(2Pi) [Tr(dsR G Gamma3_a G(w+nu,p+q) Gamma3_b G) - 1/2 Tr(dsR G Gamma4_ab G)], Int d^dq/(2Pi)^d left to the numerics",
      "Checks" -> checks, "Warnings" -> warns, "Timing" -> AbsoluteTime[] - t0|>];
   Print["Loops: ", stage["Formula"]];
   PrintCard[stage];
   stage]];

End[];
EndPackage[];
