(* ::Package:: *)

(* runNotebook.wl -- evaluates the Input cells of an example notebook in order, without a front end.
   Run:  math -script ~/sources/FRGscalarMSR/examples/runNotebook.wl ~/sources/FRGscalarMSR/examples/KPZ.nb
   Exit code 0 iff every cell evaluated without abort (a failed library check exits 2 by itself). *)

file = Last[$CommandLine];
inputs = Cases[Import[file, "NB"], Cell[s_String, "Input", ___] :> s, Infinity];
Print["--- ", FileNameTake[file], ": ", Length[inputs], " input cells"];
Do[Module[{r},
   Print["In[", i, "]:= ", StringReplace[inputs[[i]], "\n" -> " "] // StringTake[#, UpTo[150]] &];
   r = CheckAbort[ToExpression[inputs[[i]]], $Aborted];
   If[r === $Aborted || r === $Failed, Print["--- cell ", i, " failed"]; Exit[1]];
   If[r =!= Null, Print["Out[", i, "]= ", StringTake[ToString[r, InputForm], UpTo[200]]]]],
  {i, Length[inputs]}];
Print["--- all ", Length[inputs], " cells evaluated"];
Exit[0];
