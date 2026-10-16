(* ::Package:: *)

(* FRGscalarMSR -- Wetterich flow equations for N MSR scalars, derivative-expansion pipeline.
   Load from a notebook or a script with
       Get["~/sources/FRGscalarMSR/FRGscalarMSR.wl"]
   Modules are documented in docs/MODULES.md, conventions in docs/CONVENTIONS.md. *)

(* reloading: drop the previous definitions (several symbols are Protected) *)
If[MemberQ[Contexts[], "FRGscalarMSR`"],
  Unprotect["FRGscalarMSR`*", "FRGscalarMSR`Internal`*"];
  ClearAll["FRGscalarMSR`*", "FRGscalarMSR`Internal`*"]];

FRGscalarMSR`$FRGRoot = DirectoryName[$InputFileName];
FRGscalarMSR`$FRGVersion = "0.1.0 (milestone 5)";

(* load order matters: each module only uses the ones above it *)
FRGscalarMSR`$FRGModules = {"Core", "Theory", "Ansatz", "Vertices", "Regulator", "Propagator", "Symmetries", "Loops", "Projection", "FrequencyIntegration", "Closure", "Assemble", "Oracle", "Stages"};

Scan[Get[FileNameJoin[{FRGscalarMSR`$FRGRoot, "Kernel", # <> ".wl"}]] &, FRGscalarMSR`$FRGModules];

Print["FRGscalarMSR ", FRGscalarMSR`$FRGVersion, " -- modules: ", StringRiffle[FRGscalarMSR`$FRGModules, ", "]];
Print["WARNING: Ward identities are assumed (user-declared), not derived."];
