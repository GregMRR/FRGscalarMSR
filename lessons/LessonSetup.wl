(* ::Package:: *)

(* LessonSetup.wl -- helpers shared by the lesson notebooks (lessons/NN-Name.nb).
   lessonUpTo[stage] rebuilds the KPZ running example up to stage, with the cards hidden;
   showFailure[expr] evaluates expr, lets a failed library check print its FAIL block, and returns the check's name.
   Load after the library:  Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]] *)

(* the KPZ chain of tests/reference/KPZSetup.wl, in the load order of the modules (Symmetries is independent of
   Vertices..Propagator, so its place in the chain does not matter) *)
$lessonStages = {"Theory", "Ansatz", "Vertices", "Regulator", "Propagator", "Symmetries", "Loops", "Projection",
   "FrequencyIntegration", "Closure", "Equations"};

$lessonSteps = {
   Hold[kpzTheory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]],
   Hold[kpzAnsatz = DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2,
      "Couplings" -> {nuk, lamk, Dk}]],
   Hold[kpzVertices = BuildVertices[kpzAnsatz]],
   Hold[kpzRegulator = DefineRegulator[kpzAnsatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2],
      "Shapes" -> {rn, rd}]],
   Hold[kpzPropagator = BuildPropagator[kpzVertices, kpzRegulator]],
   Hold[kpzSymmetries = DeclareSymmetries[kpzAnsatz,
      Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
      Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
      Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
      "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]],
   Hold[kpzLoops = BuildLoopIntegrands[kpzVertices, kpzPropagator, kpzRegulator]],
   Hold[kpzProjected = ProjectOnMonomials[kpzLoops, kpzSymmetries]],
   Hold[kpzIntegrated = IntegrateFrequencies[kpzProjected]],
   Hold[kpzClosure = CheckClosure[kpzIntegrated, kpzSymmetries]],
   Hold[kpzEquations = AssembleEquations[kpzClosure]]};

Options[lessonUpTo] = {"Verbose" -> False};
lessonUpTo[stage_String, OptionsPattern[]] := Module[{n = First[FirstPosition[$lessonStages, stage, {0}]], r},
   If[n == 0, Print["lessonUpTo: unknown stage ", stage, "; use one of ", $lessonStages]; Return[$Failed, Module]];
   r = CheckAbort[
     If[TrueQ[OptionValue["Verbose"]], Scan[ReleaseHold, Take[$lessonSteps, n]], Block[{Print}, Scan[ReleaseHold, Take[$lessonSteps, n]]]];
     "done", $Aborted];
   If[r === $Aborted,
    Print["lessonUpTo: the KPZ chain stopped (check ", $FRGLastFailure, "); rerun lessonUpTo[\"", stage, "\", \"Verbose\" -> True] to see the FAIL block"];
    Return[$Failed, Module]];
   (* Extract[Hold[x = rhs], {1, 1}, Hold] is Hold[x]: the assigned variable, kept unevaluated *)
   Print["KPZ rebuilt up to ", stage, ": ", StringRiffle[ToString[Extract[#, {1, 1}, HoldForm]] & /@ Take[$lessonSteps, n], ", "],
    " (cards hidden)"];
   ReleaseHold[Extract[$lessonSteps[[n]], {1, 1}, Hold]]];

(* A failed check calls Exit[2] in a math -script run ($FRGExitOnFail); Block turns that into Abort[] for the duration of
   expr, and CheckAbort catches the abort, so a lesson cell can show a FAIL block without stopping the notebook. *)
SetAttributes[showFailure, HoldAll];
showFailure[expr_] := Block[{$FRGExitOnFail = False}, $FRGLastFailure = None;
   CheckAbort[expr; Print["showFailure: no check failed"]; None, $FRGLastFailure]];
