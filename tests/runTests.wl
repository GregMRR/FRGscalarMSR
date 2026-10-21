(* ::Package:: *)

(* runTests.wl -- unit and acceptance tests.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/tests/runTests.wl [--verbose] [--require-julia] [--only Name]
   Exit code 0 iff every test passes. Without --verbose the cards and FAIL blocks printed by the library are hidden.
   Tests that need julia or FlowRG are reported SKIP when those are missing; --require-julia turns a SKIP into a failure. *)

$HistoryLength = 0;
testRoot = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[testRoot], "FRGscalarMSR.wl"}]];
verboseQ = MemberQ[$CommandLine, "--verbose"];
requireJuliaQ = MemberQ[$CommandLine, "--require-julia"];
$FRGExitOnFail = False;   (* failed checks abort (caught below) instead of exiting the kernel *)

results = {}; skipped = {};
SetAttributes[quiet, HoldAll];
quiet[e_] := If[verboseQ, e, Block[{Print}, e]];
record[name_, ok_, msg_] := (AppendTo[results, {name, TrueQ[ok]}];
   If[! TrueQ[ok] || verboseQ, Print[If[TrueQ[ok], "  ok    ", "  FAIL  "], name, If[TrueQ[ok], "", "   " <> msg]]]);

(* skip[name, why]: a test that cannot run here (julia, FlowRG missing); a failure with --require-julia *)
skip[name_String, why_String] := (AppendTo[skipped, name]; Print["  SKIP  ", name, "   ", why];
   If[requireJuliaQ, record[name, False, "skipped with --require-julia: " <> why]]);

(* test[name, actual, expected]: actual must evaluate (without aborting) to expected *)
SetAttributes[test, HoldRest];
test[name_String, actual_, expected_] := Module[{a = CheckAbort[quiet[actual], $Aborted], e = expected},
   record[name, a === e, "got " <> ToString[Short[InputForm[a], 3]] <> ", expected " <> ToString[Short[InputForm[e], 3]]]];

(* expectFail[name, expr, check]: expr must abort through the library check named check *)
SetAttributes[expectFail, HoldRest];
expectFail[name_String, expr_, check_String] := Module[{r},
   $FRGLastFailure = None;
   r = CheckAbort[quiet[expr]; "returned", "aborted"];
   record[name, r === "aborted" && $FRGLastFailure === check,
    "outcome " <> r <> ", failed check " <> ToString[$FRGLastFailure] <> " (expected " <> check <> ")"]];

testFiles = Join[FileNames["*Test.wl", FileNameJoin[{testRoot, "unit"}]], FileNames["*.wl", FileNameJoin[{testRoot, "acceptance"}]]];
(* --only Name: run only the test files whose base name contains Name (e.g. --only Dimensions, --only KPZ) *)
With[{i = FirstPosition[$CommandLine, "--only"]}, If[! MissingQ[i] && Length[$CommandLine] > First[i],
   testFiles = Select[testFiles, StringContainsQ[FileBaseName[#], $CommandLine[[First[i] + 1]]] &]]];
Do[Module[{n0 = Length[results]},
   Print["--- ", FileBaseName[f]];
   CheckAbort[Get[f], record[FileBaseName[f] <> " (file aborted)", False, "uncaught abort"]];
   Print["    ", Count[Drop[results, n0], {_, True}], "/", Length[results] - n0, " passed"]], {f, testFiles}];

nFail = Count[results, {_, False}];
Print["=== ", Length[results] - nFail, "/", Length[results], " tests passed",
  If[nFail == 0, "", "; FAILED: " <> StringRiffle[Cases[results, {n_, False} :> n], ", "]],
  If[skipped === {}, "", "; " <> ToString[Length[skipped]] <> " SKIPPED"], " ==="];
Exit[If[nFail == 0 && results =!= {}, 0, 1]];
