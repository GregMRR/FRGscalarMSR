(* ::Package:: *)

(* Core.wl -- assertions, exact zero tests, stage bookkeeping, summary cards, scal.
   Objective, assumptions and checks: docs/MODULES.md section "Core". *)

(* ---- internal helpers shared by the modules (not on the user's $ContextPath) ---- *)
BeginPackage["FRGscalarMSR`Internal`"];
assertCheck::usage = "assertCheck[name, cond, details, hint, summary] runs FRGAssert and returns a check record for a stage card.";
infoCheck::usage = "infoCheck[status, name, summary] is a non-failing card line (\"WARN\", \"INFO\", \"ASSUMPTION\").";
stampStage::usage = "stampStage[assoc] adds the \"Fingerprint\" key (hash of the stage content).";
validateStage::usage = "validateStage[stage, name] asserts the stage key and fingerprint of an upstream stage.";
reservedNames::usage = "reservedNames: symbol names the user may not use for fields or couplings.";
assertFreshSymbols::usage = "assertFreshSymbols[check, syms, what] asserts that user symbols are usable as fields / couplings.";
inputForm::usage = "inputForm[e, n] is a short one-line InputForm string of e.";
displayString::usage = "displayString[e] is a full one-line InputForm string of e without package contexts; scal[a, a] -> a^2, scal[a, b] -> a.b, \
wedge[a, b] -> wedge(a,b).";
reduceWedge2D::usage = "reduceWedge2D[e] rewrites every product of two wedges by the d = 2 identity \
wedge[a,b] wedge[c,d] = scal[a,c] scal[b,d] - scal[a,d] scal[b,c] until each term of the expanded e has at most one wedge.";
EndPackage[];

BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}];

FRGAssert::usage = "FRGAssert[name, cond, details, hint] returns name if cond is True; otherwise prints a red FAIL block \
(check, module, short offending expression, hint), sets $FRGLastFailure = name and aborts (Exit[2] in a math -script run).";
$FRGLastFailure::usage = "$FRGLastFailure is the name of the last check that failed (diagnostic only).";
$FRGModule::usage = "$FRGModule is the module name printed in FAIL blocks (set with Block by each constructor).";
$FRGExitOnFail::usage = "$FRGExitOnFail: if True a failed check calls Exit[2] instead of Abort[] (default: True in math -script runs).";
ExactlyZeroQ::usage = "ExactlyZeroQ[expr] is True iff expr vanishes exactly at several random rational points: every non-System \
symbol (couplings, q, k, Dim, ...) gets an independent random rational, and every remaining function atom r[y], r'[y], r''[y], \
f[...] gets its own independent random rational (same atom, same value). Option \"Seeds\" -> {1, 2, 3}.";
RandomRationalPoint::usage = "RandomRationalPoint[syms, seed] gives rules sym -> random rational in (0, 1), reproducible from seed.";
PrintCard::usage = "PrintCard[stage] prints the stage summary card: PASS/WARN/INFO/ASSUMPTION lines, LeafCount, timing, fingerprint.";
scal::usage = "scal[a, b] is the scalar product of the vectors a, b (bare symbols). Orderless; linear over sums and over \
explicit numbers only (not over couplings).";
wedge::usage = "wedge[a, b] = eps_ij a_i b_j (d = 2, eps_12 = +1) of the vectors a, b (bare symbols). Antisymmetric (canonical \
order, not Orderless), wedge[a, a] = 0; linear over sums and over explicit numbers only.";

Begin["`Private`"];

(* ---------------- failures ---------------- *)

$FRGLastFailure = None;
$FRGModule = "FRGscalarMSR";
$FRGExitOnFail = ($FrontEnd === Null && MemberQ[$CommandLine, "-script"]);

inputForm[e_, n_ : 3] := StringReplace[ToString[Short[InputForm[e], n]], "\n" -> " "];

displayString[e_] := StringReplace[
   StringReplace[ToString[e, InputForm], RegularExpression["FRGscalarMSR`([A-Za-z]+`)*"] -> ""],
   {RegularExpression["scal\\[(\\w+), \\1\\]"] -> "$1^2", RegularExpression["scal\\[(\\w+), (\\w+)\\]"] -> "$1.$2",
    RegularExpression["wedge\\[([^\\[\\],]+), ([^\\[\\],]+)\\]"] -> "wedge($1,$2)"}];

failBlock[name_, details_, hint_] := Module[{lines},
   lines = {"FAIL  check: " <> ToString[name] <> "   module: " <> ToString[$FRGModule],
     "      offending: " <> inputForm[details, 4],
     "      hint: " <> ToString[hint]};
   If[$FrontEnd === Null,
    Scan[Print, lines],
    Print[Framed[Column[Style[#, Red, If[# === First[lines], Bold, Plain]] & /@ lines], FrameStyle -> Red]]]];

SetAttributes[FRGAssert, HoldRest];
FRGAssert[name_, cond_, details_, hint_ : ""] :=
  If[TrueQ[cond], name,
   failBlock[name, details, hint];
   $FRGLastFailure = name;
   If[TrueQ[$FRGExitOnFail], Exit[2], Abort[]]];

SetAttributes[assertCheck, HoldRest];
assertCheck[name_, cond_, details_, hint_, summary_ : ""] := (FRGAssert[name, cond, details, hint];
   <|"Status" -> "PASS", "Check" -> name, "Detail" -> summary|>);
infoCheck[status_String, name_, summary_ : ""] := <|"Status" -> status, "Check" -> name, "Detail" -> summary|>;

(* ---------------- exact zero test (primer rule 7: never Simplify[...] == 0) ---------------- *)

RandomRationalPoint[syms_List, seed_Integer] := (SeedRandom[seed];
   Thread[syms -> RandomInteger[{2, 60}, Length[syms]]/RandomInteger[{61, 97}, Length[syms]]]);

systemSymbolQ[s_Symbol] := Context[s] === "System`";
userSymbols[e_] := Union[Cases[e, s_Symbol /; ! systemSymbolQ[s], {0, Infinity}, Heads -> False]];

valueAt[e_, seed_Integer] := Module[{pt, memo, atomRule},
   pt = RandomRationalPoint[userSymbols[e], seed];
   SeedRandom[seed + 7919];
   memo[x_] := memo[x] = RandomInteger[{2, 60}]/RandomInteger[{61, 97}];
   atomRule = x : (Derivative[__][_Symbol][___] | (h_Symbol /; ! systemSymbolQ[h])[___]) :> memo[x];
   Quiet[Check[Together[(e /. pt) /. atomRule], $Failed], {Power::infy, Infinity::indet}]];

Options[ExactlyZeroQ] = {"Seeds" -> {1, 2, 3}};
ExactlyZeroQ[l_List, opts : OptionsPattern[]] := AllTrue[Flatten[l], ExactlyZeroQ[#, opts] &];
ExactlyZeroQ[e_, OptionsPattern[]] := Module[{vals},
   vals = Table[Module[{v = valueAt[e, sd]}, If[v === $Failed, valueAt[e, sd + 1000], v]], {sd, OptionValue["Seeds"]}];
   Block[{$FRGModule = "Core"},
    FRGAssert["ExactlyZeroQ.numeric", AllTrue[vals, NumericQ[#] && Precision[#] === Infinity &], Select[vals, ! NumericQ[#] &],
     "after substituting random rationals the expression is not an exact number: an unknown head or a pattern is left"]];
   AllTrue[vals, PossibleZeroQ]];

(* ---------------- stages ---------------- *)

(* keys that do not define the content of a stage, hence are not hashed *)
volatileKeys = {"Fingerprint", "Checks", "Timing", "Warnings", "Assumptions"};
(* Normal re-evaluates the content: a field or coupling that got a value after the stage was built changes the hash *)
contentHash[stage_Association] := Hash[Normal[KeyDrop[stage, volatileKeys], Association], "SHA256", "HexString"];
stampStage[stage_Association] := Append[stage, "Fingerprint" -> contentHash[stage]];

validateStage[stage_, name_String] := (
   FRGAssert["Stage.input", AssociationQ[stage] && Lookup[stage, "Stage", None] === name, Lookup[If[AssociationQ[stage], stage, <||>], "Stage", stage],
    "expected the output of the " <> name <> " stage"];
   FRGAssert["Stage.fingerprint", contentHash[stage] === stage["Fingerprint"], stage["Stage"],
    "the stage was modified after it was built, or one of its symbols got a value since: rebuild it"];
   stage);

PrintCard[stage_Association] := Module[{chk = Lookup[stage, "Checks", {}], extra},
   extra = Join[Lookup[stage, "Warnings", {}], Lookup[stage, "Assumptions", {}]];
   Print["==== ", stage["Stage"], " card   fingerprint ", StringTake[Lookup[stage, "Fingerprint", "--------"], UpTo[8]],
    "   LeafCount ", LeafCount[KeyDrop[stage, volatileKeys]], "   ", Round[1000 Lookup[stage, "Timing", 0]], " ms ===="];
   Scan[Print[StringPadRight[#Status, 11], #Check, If[#Detail === "", "", "   " <> ToString[#Detail]]] &, Join[chk, extra]];
   Print["==== ", Count[chk, KeyValuePattern["Status" -> "PASS"]], " PASS, ", Length[extra], " warnings/assumptions ===="];];

(* ---------------- user symbols ---------------- *)

(* names the user writes or reads in library input/output with a fixed meaning: RG scale k, loop momentum q, external momentum p,
   shape argument y = q^2/k^2, time t, position x, frequencies; internal names (A, CC, X, ...) live in private contexts.
   Dim is user-chosen and therefore not reserved. *)
reservedNames = {"k", "q", "p", "y", "t", "x", "\[Omega]", "\[Nu]"};

assertFreshSymbols[check_String, syms_List, what_String] := Module[{lib = Join[Names["FRGscalarMSR`*"], Names["FRGscalarMSR`Internal`*"]]},
   FRGAssert[check <> ".symbols", AllTrue[syms, Head[#] === Symbol &], Select[syms, Head[#] =!= Symbol &],
    what <> " must be plain symbols without values (one of them evaluated to something else)"];
   FRGAssert[check <> ".system", NoneTrue[syms, systemSymbolQ], Select[syms, systemSymbolQ], what <> " must not be System symbols"];
   FRGAssert[check <> ".reserved", NoneTrue[syms, MemberQ[reservedNames, SymbolName[#]] &], Select[syms, MemberQ[reservedNames, SymbolName[#]] &],
    "these names are reserved (RG scale k, momenta, frequencies, Delta coefficients): rename them"];
   FRGAssert[check <> ".library", NoneTrue[syms, MemberQ[lib, SymbolName[#]] &], Select[syms, MemberQ[lib, SymbolName[#]] &],
    "these names collide with FRGscalarMSR symbols: rename them"];
   FRGAssert[check <> ".values", NoneTrue[syms, (OwnValues[#] =!= {} || DownValues[#] =!= {} || UpValues[#] =!= {}) &],
    Select[syms, (OwnValues[#] =!= {} || DownValues[#] =!= {} || UpValues[#] =!= {}) &], what <> " must have no definitions: ClearAll them"];
   FRGAssert[check <> ".duplicates", DuplicateFreeQ[syms], syms, what <> " are listed twice"];
   syms];

(* ---------------- scalar products (DE4plus_fast.wl section 1, E14: linear over numbers only) ---------------- *)

SetAttributes[scal, Orderless];
scal[(a_) + (b_), c_] := scal[a, c] + scal[b, c];
scal[(a_)*(b_), c_] /; NumberQ[a] := a*scal[b, c];
scal[-(a_), b_] := -scal[a, b];
scal[0, b_] := 0;

(* ---------------- 2D antisymmetric product eps_ij a_i b_j (design/LEVICIVITA.md F1, F2) ---------------- *)

wedge[0, _] := 0;
wedge[_, 0] := 0;
wedge[a_, a_] := 0;
wedge[(a_) + (b_), c_] := wedge[a, c] + wedge[b, c];
wedge[a_, (b_) + (c_)] := wedge[a, b] + wedge[a, c];
wedge[(a_)*(b_), c_] /; NumberQ[a] := a*wedge[b, c];
wedge[a_, (b_)*(c_)] /; NumberQ[b] := b*wedge[a, c];
wedge[b_, a_] /; ! OrderedQ[{b, a}] := -wedge[a, b];

(* F2, valid in d = 2 only *)
wedgePairRules = {
   wedge[a_, b_] wedge[c_, d_] :> scal[a, c] scal[b, d] - scal[a, d] scal[b, c],
   Power[wedge[a_, b_], n_Integer /; n >= 2] :> (scal[a, a] scal[b, b] - scal[a, b]^2) wedge[a, b]^(n - 2)};
reduceWedge2D[e_] := FixedPoint[Expand[# /. wedgePairRules] &, Expand[e]];

End[];
EndPackage[];
