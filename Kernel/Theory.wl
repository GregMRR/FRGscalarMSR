(* ::Package:: *)

(* Theory.wl -- DefineTheory: N physical scalars, their MSR response fields, the field order, the space dimension.
   Objective, assumptions and checks: docs/MODULES.md section "Theory". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

DefineTheory::usage = "DefineTheory[\"Scalars\" -> {phi1, phi2, ...}, \"Dimension\" -> Dim] declares N physical scalars; \
the response field of phi is created as the symbol phibar (same context). Field order: (phi1, phi1bar, phi2, phi2bar, ...). \
Returns the \"Theory\" stage.";
ResponseField::usage = "ResponseField[theory, phi] is the MSR response field of the scalar phi.";

Begin["`Private`"];

Options[DefineTheory] = {"Scalars" -> {}, "Dimension" -> None};

DefineTheory[opts : OptionsPattern[]] := Block[{$FRGModule = "Theory"}, Module[{t0 = AbsoluteTime[], phis, bars, dim, checks, stage},
    FRGAssert["Theory.options", MatchQ[{opts}, {(Rule | RuleDelayed)["Scalars" | "Dimension", _] ...}], {opts},
     "only the options \"Scalars\" and \"Dimension\" are accepted"];
    phis = OptionValue["Scalars"]; dim = OptionValue["Dimension"];
    FRGAssert["Theory.scalarList", ListQ[phis] && phis =!= {}, phis, "\"Scalars\" -> {H, Psi, ...} is required"];
    checks = {assertCheck["Theory.N", 1 <= Length[phis] <= 2, Length[phis],
       "v1 supports N = 1 or 2 scalars: the propagator denominator must have degree <= 2 in w (see docs/FUTURE.md)",
       "N = " <> ToString[Length[phis]] <> " (supported: 1, 2)"]};
    assertFreshSymbols["Theory.scalars", phis, "scalars"];
    bars = Symbol[Context[#] <> SymbolName[#] <> "bar"] & /@ phis;
    assertFreshSymbols["Theory.responseFields", bars, "response fields (auto-generated <name>bar)"];
    FRGAssert["Theory.dimension", Head[dim] === Symbol && dim =!= None, dim, "\"Dimension\" -> Dim is required: d stays a symbol"];
    assertFreshSymbols["Theory.dimension", {dim}, "the dimension symbol"];
    FRGAssert["Theory.dimensionDistinct", ! MemberQ[Join[phis, bars], dim], dim, "the dimension symbol must not be a field"];
    AppendTo[checks, assertCheck["Theory.symbols", True, Null, "", "fields and Dim are fresh and collision-free"]];
    stage = stampStage[<|"Stage" -> "Theory", "Pipeline" -> "DerivativeExpansion", "Upstream" -> {},
       "Scalars" -> phis, "ResponseFields" -> bars, "Fields" -> Flatten[Transpose[{phis, bars}]],
       "ResponseOf" -> AssociationThread[phis -> bars], "Dimension" -> dim,
       "Checks" -> checks, "Timing" -> AbsoluteTime[] - t0|>];
    Print["Theory: fields in order ", stage["Fields"], ", dimension ", dim, " (response fields are created automatically)"];
    PrintCard[stage];
    stage]];

ResponseField[theory_Association, phi_Symbol] := Block[{$FRGModule = "Theory"},
   FRGAssert["Theory.isScalar", MemberQ[theory["Scalars"], phi], phi, "not a declared physical scalar"];
   theory["ResponseOf"][phi]];

End[];
EndPackage[];
