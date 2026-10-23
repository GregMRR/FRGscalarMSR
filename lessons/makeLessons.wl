(* ::Package:: *)

(* makeLessons.wl -- writes lessons/NN-Name.nb and lessons/Glossary.nb from the sources lessons/src/NN-Name.md.
   Run:  math -script ~/Physics/Sources/FRGscalarMSR/lessons/makeLessons.wl
   Check one lesson headless:  math -script ~/Physics/Sources/FRGscalarMSR/examples/runNotebook.wl lessons/NN-Name.nb
   Edit the sources, not the .nb files, then regenerate. The source format is described in lessons/README.md.
   The build stops (exit 1) when a quoted symbol is missing from the code, when a [[glossary reference]] comes before the
   concept is explained, or when a concept is explained twice. *)

lessonsDir = DirectoryName[$InputFileName];
repoRoot = ParentDirectory[lessonsDir];
buildFail[msg__] := (Print["makeLessons: ", msg]; Exit[1]);
(* sources and code are UTF-8 (Greek letters, arrows); ReadString alone would decode them as Latin-1 *)
readUTF8[f_] := ByteArrayToString[ReadByteArray[f], "UTF-8"];

(* ---------- styles ---------- *)

(* inline code: raw text (a string in TextData is never typeset as input), monospace, coloured and shaded so it stands out *)
codeBox[c_String] := StyleBox[c, FontFamily -> "Source Code Pro", FontColor -> RGBColor[0.6, 0.12, 0.1], Background -> GrayLevel[0.92]];
boxStyle["why"] = {Background -> RGBColor[0.93, 0.93, 0.93], CellFrame -> {{4, 0}, {0, 0}}, CellFrameColor -> GrayLevel[0.55]};
boxStyle["mma"] = {Background -> RGBColor[0.91, 0.95, 1.], CellFrame -> {{4, 0}, {0, 0}}, CellFrameColor -> RGBColor[0.3, 0.5, 0.85]};
boxTitle["why"] = "Why it's built this way: ";
boxTitle["mma"] = "Mathematica: ";

(* ---------- glossary: concept -> {lesson label, one-line summary}, in order of first explanation ---------- *)

glossary = <||>;
currentLesson = "";

glossaryRef[c_String] := (If[! KeyExistsQ[glossary, c], buildFail["lesson ", currentLesson, ": [[", c, "]] is used before it is explained (no :::mma ", c, " block so far)"]];
   (* a list, flattened by inline: a Sequence returned from (...; ...) would splice into CompoundExpression and lose its first part *)
   {codeBox[c], StyleBox[" (Glossary, lesson " <> glossary[c][[1]] <> ")", FontSlant -> Italic, FontColor -> GrayLevel[0.4]]});

(* ---------- inline markup: `code`, `` code with a ` inside `` (Markdown's double backticks), **bold**, *word* (italic),
   [[glossary concept]] ---------- *)

inline[s_String] := Module[{parts},
   parts = StringSplit[s, {
      "[[" ~~ c : Except["]"] .. ~~ "]]" :> gref[c],
      "`` " ~~ c : Shortest[__] ~~ " ``" :> codeBox[c],
      "`" ~~ c : Except["`"] .. ~~ "`" :> codeBox[c],
      "**" ~~ c : Except["*"] .. ~~ "**" :> StyleBox[c, FontWeight -> Bold],
      "*" ~~ c : LetterCharacter .. ~~ "*" :> StyleBox[c, FontSlant -> Italic]}];
   Flatten[mathText /@ Flatten[parts /. gref[c_] :> glossaryRef[c]]]];

(* non-ASCII runs of prose (Greek, arrows, ∫, ∂, h̄) in Noto Sans: the default text font falls back to Mathematica's math font for
   them, drawn noticeably larger than the text around; a letter followed by a combining mark (h + overbar) stays in the run *)
unicodeRun = RegularExpression["(?:[A-Za-z](?=[\\x{0300}-\\x{036F}]))?[^\\x{00}-\\x{7F}]+"];
mathText[t_String] := StringSplit[t, u : unicodeRun :> StyleBox[u, FontFamily -> "Noto Sans"]];
mathText[e_] := e;

textCell[s_String, opts___] := Cell[TextData[inline[s]], "Text", opts];

(* ---------- code quotes, read from the source files at build time ---------- *)

fileLines[path_] := Module[{f = FileNameJoin[{repoRoot, path}]},
   If[! FileExistsQ[f], buildFail["lesson ", currentLesson, ": file not found: ", path]];
   StringSplit[readUTF8[f], "\n", All]];

programCell[path_, {s_, e_}, lines_] := Cell["(* " <> path <> ":" <> ToString[s] <> If[e > s, "-" <> ToString[e], ""] <> " *)\n" <>
    StringRiffle[lines[[s ;; e]], "\n"], "Program"];

(* the definition blocks of sym: a line at column 0 that defines sym (sym[..] :=, sym = .., Options[sym], SetAttributes[sym or
   {.., sym, ..}]) and the indented or blank lines after it; a comment directly above (no blank line in between) is included *)
definitionStarts[lines_, sym_] := Flatten[Position[lines, l_String /; StringMatchQ[l,
      (sym ~~ ("[" | " =" | " :=" | "=" | ":=") ~~ ___) | ("Options[" ~~ sym ~~ "]" ~~ ___) |
      ("SetAttributes[" ~~ (sym ~~ "," | ("{" ~~ ___ ~~ WordBoundary ~~ sym ~~ WordBoundary ~~ ___ ~~ "}")) ~~ ___)], {1}, Heads -> False]];
blockEnd[lines_, s_] := Module[{e = s},
   While[e < Length[lines] && (lines[[e + 1]] === "" || StringStartsQ[lines[[e + 1]], " " | "\t"]), e++];
   While[e > s && StringTrim[lines[[e]]] === "", e--];
   e];
commentAbove[lines_, s_] := Module[{b = s},
   While[b > 1 && StringTrim[lines[[b - 1]]] =!= "", b--];
   If[b < s && StringStartsQ[lines[[b]], "(*"] && StringEndsQ[StringTrim[lines[[s - 1]]], "*)"], b, s]];
mergeIntervals[iv_] := Fold[If[#1 =!= {} && #2[[1]] <= Last[#1][[2]] + 1, Append[Most[#1], {Last[#1][[1]], Max[Last[#1][[2]], #2[[2]]]}], Append[#1, #2]] &,
   {}, Sort[iv]];

quoteCells[path_, syms_List] := Module[{lines = fileLines[path], iv},
   iv = Flatten[Function[sym, Module[{st = definitionStarts[lines, sym]},
        If[st === {}, buildFail["lesson ", currentLesson, ": no definition of ", sym, " in ", path]];
        {commentAbove[lines, #], blockEnd[lines, #]} & /@ st]] /@ syms, 1];
   programCell[path, #, lines] & /@ mergeIntervals[iv]];

regionCells[path_, from_String, to_String] := Module[{lines = fileLines[path], s, e},
   s = FirstPosition[lines, l_String /; StringContainsQ[l, from], {0}, {1}, Heads -> False][[1]];
   If[s == 0, buildFail["lesson ", currentLesson, ": \"", from, "\" not found in ", path]];
   e = FirstPosition[lines[[s ;;]], l_String /; StringContainsQ[l, to], {0}, {1}, Heads -> False][[1]];
   If[e == 0, buildFail["lesson ", currentLesson, ": \"", to, "\" not found after \"", from, "\" in ", path]];
   {programCell[path, {s, s + e - 1}, lines]}];

fileCells[path_] := Module[{lines = fileLines[path], e},
   e = Length[lines]; While[e > 1 && StringTrim[lines[[e]]] === "", e--];
   {programCell[path, {1, e}, lines]}];

directiveCells[l_String] := Module[{w = StringSplit[StringDrop[l, 1]]},
   Switch[First[w],
    "quote", quoteCells[w[[2]], Drop[w, 2]],
    "file", fileCells[w[[2]]],
    "region", With[{p = StringTrim /@ StringSplit[StringDrop[l, StringLength["@region"]], "|"]}, regionCells @@ p],
    _, buildFail["lesson ", currentLesson, ": unknown directive ", l]]];

(* ---------- blocks ---------- *)

headingCell[l_String] := Module[{lev = StringLength[First[StringCases[l, StartOfString ~~ "#" ..]]], t},
   t = StringReplace[StringTrim[StringDrop[l, lev]], "`" -> ""];
   Cell[t, {"Title", "Section", "Subsection", "Subsubsection"}[[lev]]]];

(* a box (:::why, :::mma) is one cell: title, then its paragraphs and bullets; its code fences follow as Input cells *)
boxCells[kind_, args_, body_List] := Module[{cells = parseBlocks[body], texts, codes, head, concept, summary},
   If[kind === "mma",
    {concept, summary} = StringTrim /@ PadRight[StringSplit[args, "|", 2], 2, ""];
    If[summary === "", buildFail["lesson ", currentLesson, ": :::mma ", concept, " needs \"| one-line summary\""]];
    If[KeyExistsQ[glossary, concept], buildFail["lesson ", currentLesson, ": ", concept, " is explained twice (first in lesson ",
      glossary[concept][[1]], "); use [[", concept, "]] instead"]];
    glossary[concept] = {currentLesson, summary}; head = concept,
    head = args];
   texts = Cases[cells, Cell[TextData[t_], st : "Text" | "Item", ___] :> If[st === "Item", Join[{"\n\[Bullet] "}, Flatten[{t}]], Join[{"\n\n"}, Flatten[{t}]]]];
   (* the box's code examples carry its colours, so they read as part of it *)
   codes = Cases[cells, Cell[c_, st : "Input" | "Program", ___] :> Cell[c, st, Sequence @@ boxStyle[kind]]];
   If[Length[texts] + Length[codes] =!= Length[cells], buildFail["lesson ", currentLesson, ": only text, bullets and code may appear in a ::: ", kind, " box"]];
   (* an mma box opens with its one-line summary, the same line the Glossary shows *)
   If[kind === "mma", texts = Prepend[texts, Join[{"\n"}, inline[summary]]]];
   Prepend[codes, Cell[TextData[Join[{StyleBox[boxTitle[kind] <> head, FontWeight -> Bold]}, Flatten[texts]]], "Text", Sequence @@ boxStyle[kind]]]];

exerciseCells[args_, body_] := Prepend[parseBlocks[body], Cell["Exercise: " <> args, "Subsection"]];
solutionCells[body_] := {Cell[CellGroupData[Prepend[parseBlocks[body], Cell["Solution", "Subsubsection"]], Closed]]};

closingLine[lines_, i_] := Module[{j = i + 1},
   While[j <= Length[lines] && StringTrim[lines[[j]]] =!= ":::", If[StringStartsQ[lines[[j]], ":::" ~~ LetterCharacter], buildFail["lesson ", currentLesson,
      ": ::: blocks do not nest (line ", lines[[i]], ")"]]; j++];
   If[j > Length[lines], buildFail["lesson ", currentLesson, ": unclosed ", lines[[i]]]]; j];
fenceEnd[lines_, i_] := Module[{j = i + 1},
   While[j <= Length[lines] && ! StringStartsQ[lines[[j]], "```"], j++];
   If[j > Length[lines], buildFail["lesson ", currentLesson, ": unclosed code fence ", lines[[i]]]]; j];

(* a closed Solution group must end its exercise: text after it would be folded into it when the front end regroups cells *)
checkAfterSolution[lines_, j_] := Module[{n = j + 1},
   While[n <= Length[lines] && StringTrim[lines[[n]]] === "", n++];
   If[n <= Length[lines] && ! StringStartsQ[lines[[n]], "#" | ":::exercise"],
    buildFail["lesson ", currentLesson, ": a :::solution block must be followed by a heading or another :::exercise, not: ", lines[[n]]]]];

parseBlocks[lines_List] := Module[{i = 1, n = Length[lines], out = {}, para = {}, flush, l, j, kind, args, item},
   flush[] := (If[para =!= {}, AppendTo[out, textCell[StringRiffle[para, " "]]]]; para = {});
   While[i <= n,
    l = lines[[i]];
    Which[
     StringTrim[l] === "", flush[]; i++,
     StringStartsQ[l, "```"], flush[]; j = fenceEnd[lines, i];
      AppendTo[out, Cell[StringRiffle[lines[[i + 1 ;; j - 1]], "\n"], If[StringTrim[StringDrop[l, 3]] === "text", "Program", "Input"]]]; i = j + 1,
     StringStartsQ[l, ":::"], flush[]; j = closingLine[lines, i];
      {kind, args} = PadRight[StringSplit[StringDrop[l, 3], Whitespace, 2], 2, ""];
      AppendTo[out, Switch[kind,
        "why" | "mma", boxCells[kind, args, lines[[i + 1 ;; j - 1]]],
        "exercise", exerciseCells[args, lines[[i + 1 ;; j - 1]]],
        "solution", checkAfterSolution[lines, j]; solutionCells[lines[[i + 1 ;; j - 1]]],
        _, buildFail["lesson ", currentLesson, ": unknown block ", l]]]; i = j + 1,
     StringStartsQ[l, "@"], flush[]; AppendTo[out, directiveCells[l]]; i++,
     StringStartsQ[l, "#"], flush[]; AppendTo[out, headingCell[l]]; i++,
     StringStartsQ[l, "- "], flush[]; item = StringDrop[l, 2]; i++;
      While[i <= n && StringStartsQ[lines[[i]], "  "] && StringTrim[lines[[i]]] =!= "", item = item <> " " <> StringTrim[lines[[i]]]; i++];
      AppendTo[out, Cell[TextData[inline[item]], "Item"]],
     True, AppendTo[para, StringTrim[l]]; i++]];
   flush[];
   Flatten[out]];

(* ---------- notebooks ---------- *)

(* CellContext -> Notebook: each lesson's variables live in a context of its own, so lessons opened in one kernel do not interfere *)
nb[cells_] := Notebook[cells, StyleDefinitions -> "Default.nb", CellContext -> Notebook];
inputCount[cells_] := Count[cells, Cell[_, "Input", ___], Infinity];

sources = Sort[Select[FileNames["*.md", FileNameJoin[{lessonsDir, "src"}]], StringMatchQ[FileNameTake[#], RegularExpression["\\d\\d-.+\\.md"]] &]];
If[sources === {}, buildFail["no sources in lessons/src"]];

Do[Module[{name = FileBaseName[src], cells},
   currentLesson = StringTake[name, 2];
   cells = parseBlocks[StringSplit[readUTF8[src], "\n", All]];
   Export[FileNameJoin[{lessonsDir, name <> ".nb"}], nb[cells]];
   Print["wrote ", name, ".nb (", inputCount[cells], " inputs)"]],
  {src, sources}];

glossaryCells = Join[
   {Cell["Mathematica glossary", "Title"],
    textCell["Every Mathematica command, tool or practice used in the lessons, in the order the lessons explain it. Each entry gives \
a one-line summary and the lesson whose blue \"Mathematica:\" box explains it in full. Generated by lessons/makeLessons.wl from the \
lesson sources, so it always matches them."]},
   Flatten[KeyValueMap[Function[{lesson, entries},
       {Cell["Lesson " <> lesson <> " \[LongDash] " <> StringDrop[FileBaseName[SelectFirst[sources, StringStartsQ[FileBaseName[#], lesson] &]], 3], "Section"],
        Cell[TextData[{codeBox[#[[1]]], "  \[LongDash]  ", Sequence @@ inline[#[[2]]]}], "Item"] & /@ entries}],
      GroupBy[KeyValueMap[{#1, #2[[2]], #2[[1]]} &, glossary], Last -> Most]]]];
Export[FileNameJoin[{lessonsDir, "Glossary.nb"}], nb[glossaryCells]];
Print["wrote Glossary.nb (", Length[glossary], " concepts)"];
