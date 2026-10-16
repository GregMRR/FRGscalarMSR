(* ::Package:: *)

(* Stages.wl -- SaveStage / LoadStage: any stage to a .wdx file and back, with its fingerprint checked on load and, optionally,
   the stages it must have been built from (theory, ansatz, regulator, symmetries, ...).
   Objective, assumptions and checks: docs/MODULES.md section "Stages". *)

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

SaveStage::usage = "SaveStage[stage, file] writes a stage (with every stage it embeds) to file (.wdx) and returns the file name.";
LoadStage::usage = "LoadStage[file, \"Expect\" -> {stage1, ...}] reads a stage saved by SaveStage, checks its fingerprint and that it was \
built from every stage listed in \"Expect\" (e.g. the current theory, ansatz, regulator and symmetries); a mismatch aborts.";
StageFingerprints::usage = "StageFingerprints[stage] lists the fingerprints of the stage, of every stage it embeds and of their upstream stages.";

Begin["`Private`"];

StageFingerprints[s_Association] := Union[Flatten[{Lookup[s, "Fingerprint", {}], Lookup[s, "Upstream", {}],
     StageFingerprints /@ Select[Values[s], AssociationQ[#] && KeyExistsQ[#, "Stage"] &]}]];

SaveStage[stage_, file_String] := Block[{$FRGModule = "Stages"},
   FRGAssert["Stages.input", AssociationQ[stage] && KeyExistsQ[stage, "Stage"] && KeyExistsQ[stage, "Fingerprint"], Head[stage], "SaveStage takes a stage"];
   validateStage[stage, stage["Stage"]];
   Export[file, stage, "WDX"]];

Options[LoadStage] = {"Expect" -> {}};
LoadStage[file_String, OptionsPattern[]] := Block[{$FRGModule = "Stages"},
  Module[{s, fps, exp = OptionValue["Expect"], stale},
   FRGAssert["Stages.file", FileExistsQ[file], file, "no such file"];
   s = Import[file, "WDX"];
   FRGAssert["Stages.input", AssociationQ[s] && KeyExistsQ[s, "Stage"], Head[s], "the file does not contain a stage"];
   validateStage[s, s["Stage"]];
   FRGAssert["Stages.expect", AllTrue[exp, AssociationQ[#] && KeyExistsQ[#, "Fingerprint"] &], exp, "\"Expect\" -> {stage, ...}"];
   fps = StageFingerprints[s];
   stale = Select[exp, ! MemberQ[fps, #Fingerprint] &];
   FRGAssert["Stages.stale", stale === {}, #Stage & /@ stale,
    "the saved stage was not built from the current " <> StringRiffle[#Stage & /@ stale, ", "] <> " stage(s): rebuild it"];
   Print["LoadStage: ", s["Stage"], " stage, fingerprint ", StringTake[s["Fingerprint"], 8], If[exp === {}, "",
     ", built from the current " <> StringRiffle[#Stage & /@ exp, ", "]]];
   s]];

End[];
EndPackage[];
