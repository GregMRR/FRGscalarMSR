(* ::Package:: *)

(* runTests.wl -- unit and acceptance tests.
   Run:  math -script ~/sources/FRGscalarMSR/tests/runTests.wl [--verbose]
   Exit code 0 iff every test passes. Without --verbose the cards and FAIL blocks printed by the library are hidden. *)

$HistoryLength = 0;
testRoot = DirectoryName[$InputFileName];
Get[FileNameJoin[{ParentDirectory[testRoot], "FRGscalarMSR.wl"}]];
verboseQ = MemberQ[$CommandLine, "--verbose"];
$FRGExitOnFail = False;   (* failed checks abort (caught below) instead of exiting the kernel *)

results = {};
SetAttributes[quiet, HoldAll];
quiet[e_] := If[verboseQ, e, Block[{Print}, e]];
record[name_, ok_, msg_] := (AppendTo[results, {name, TrueQ[ok]}];
   If[! TrueQ[ok] || verboseQ, Print[If[TrueQ[ok], "  ok    ", "  FAIL  "], name, If[TrueQ[ok], "", "   " <> msg]]]);

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
Do[Module[{n0 = Length[results]},
   Print["--- ", FileBaseName[f]];
   CheckAbort[Get[f], record[FileBaseName[f] <> " (file aborted)", False, "uncaught abort"]];
   Print["    ", Count[Drop[results, n0], {_, True}], "/", Length[results] - n0, " passed"]], {f, testFiles}];

nFail = Count[results, {_, False}];
Print["=== ", Length[results] - nFail, "/", Length[results], " tests passed",
  If[nFail == 0, "", "; FAILED: " <> StringRiffle[Cases[results, {n_, False} :> n], ", "]], " ==="];
Exit[If[nFail == 0 && results =!= {}, 0, 1]];
