# Lesson 15 - JuliaExport

**Physics step.** None: the dimensionless flow is written, at one integer space dimension d, as a Julia file that FlowRG loads and integrates. The file format is a contract owned by the other repository, `~/Physics/Sources/FlowRG.jl/docs/equations-interface.md`: `COUPLINGS`, `ETAS`, `SHAPES`, `MUST_FIX`, and the functions `betal`, `betanl`, `etaMl`, `etaMnl`, `etabl`, `etabnl` with positional arguments, plus `causal_poles`.

**Stage.** Input: the Adimensional stage. Output: a file (no new stage).

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
equations = lessonUpTo["Equations"];
dims = Block[{Print}, DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}, "Scale" -> k]];
adim = Block[{Print}, AdimensionalEquations[equations, dims]];
```

## 1. Using it

```wl
file = ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{$TemporaryDirectory, "lessonKPZ", "Equations2D.jl"}]];
```

```wl
ReadString[file]
```

A header that says where the file comes from and how to regenerate it; three constants; one Julia function per non-zero part. In the loop functions the dimensionless atoms are computed first, as local variables (`atom_h_hbar_0_0 = 1 + rn(aq^2)`), so the long expressions stay readable. `Drn(aq^2)` is r_n′(aq²).

`JuliaParts` gives the same content as Mathematica expressions, before the translation to text:

```wl
parts = JuliaParts[adim, 2];
{Keys[parts], parts["betal"], parts["causalPoles"]}
```

## 2. The parts at one dimension

@quote Kernel/JuliaExport.wl JuliaParts derivativeRule

`Dim` is replaced by the requested d, the loop measure is multiplied into every loop part, and each shape derivative `Derivative[1][rn][y]` becomes an applied function `Drn[y]` ([[D and unknown functions]]). Only r, r′ and r″ have names in the FlowRG interface; a p⁴ flow, which needs r‴ and r⁗ (lesson 08), aborts here with `JuliaExport.shapeDerivative`.

## 3. From expressions to Julia text

@quote Kernel/JuliaExport.wl toJulia atomName symTuple tupleString

:::mma Expressions as text | `ToString[e, InputForm, PageWidth -> Infinity]` prints e on one line in input syntax; `StringReplace[s, {"a" -> "b", ...}]` rewrites the text; `RegularExpression["..."]` matches by a regex.
Mathematica's input syntax is close to Julia's: `+ - * / ^` and parentheses are the same. The translation turns brackets into parentheses, renames constants (`Pi` → `pi`, through a placeholder symbol `juliaPi`, so that no other `Pi` in a name is touched), and removes every ``FRGscalarMSR`...`` context prefix with a regular expression, as `displayString` does (lesson 00).
```wl
{ToString[x^2 + Sin[y], InputForm], FRGscalarMSR`Private`toJulia[aq^2 rn[aq^2] + Pi]}
```
:::

`StringReplace` with a list of rules makes **one pass** over the string: at each position it tries the rules in order, replaces, and moves on; the result of one rule is never seen by another. That explains a latent problem in `toJulia`: the rule `"Sqrt(" -> "sqrt("` can never match, because at that point the text still reads `Sqrt[`, and the `"[" -> "("` rule only produces the parenthesis in the output.

```wl
{StringReplace["Sqrt[x]", {"[" -> "(", "]" -> ")", "Sqrt(" -> "sqrt("}], FRGscalarMSR`Private`toJulia[Sqrt[aq]]}
```

A `Sqrt` in an exported expression would therefore be written `Sqrt(aq)`, which Julia does not know. It does not happen in the current models, whose exported expressions are rational in aq, but `JuliaExport.systemHeads` allows `Sqrt`, `Exp`, `Log`, `Abs`, so the check would not catch it. Two passes (brackets first, then names) would fix it.

## 4. Writing the file

@region Kernel/JuliaExport.wl | (* one Julia definition; compact loop parts | Riffle[defs, ""]], "\n"] <> "\n", "Text"];

The atoms are written through placeholder symbols `atomJ1Z`, `atomJ2Z`, ..., which `toJulia` prints like any symbol, and the placeholders are renamed to `atom_h_hbar_0_0` in the final string: a symbol name cannot contain the characters the Julia name needs, but a string can.

:::mma Building text | `Riffle[list, x]` inserts x between the elements; `StringRiffle[list, sep]` joins strings; `DateString[...]` formats the date; `Export[file, string, "Text"]` writes text; `CreateDirectory` makes a folder.
`Riffle[defs, ""]` puts an empty line between the function definitions. `Quiet[CreateDirectory[dir], CreateDirectory::eexist]` ([[Check and Quiet]]) creates the folder and silences only the "already exists" message.
```wl
{Riffle[{a, b, c}, x], StringRiffle[{"a", "b"}, ", "], DateString[{"Year", "-", "Month", "-", "Day"}]}
```
:::

## 5. The checks

@region Kernel/JuliaExport.wl | (* check 7: no unit, no dimensionful loop momentum | ": load with fixed_etas = ("

Before writing, the module checks that nothing dimensionful is left (no k, reference, q, symbolic `Dim`), that every symbol is a coupling, an η, `aq`, a shape or an atom, and that the only System functions are those with a Julia translation. A file that uses r″ gets a WARN, because FlowRG must bind `D2<shape>` (`docs/handoff/FlowRG-MUST_FIX.md`).

```wl
showFailure[ExportJulia[adim, "Dimension" -> d, "File" -> FileNameJoin[{$TemporaryDirectory, "x.jl"}]]]
```

:::why A file format owned by someone else
The library generates code for a program it does not control. The contract is written down in the FlowRG repository, and any change to it goes through a hand-off spec in `docs/handoff/` (CLAUDE.md, Where to change X), so the two sides cannot drift silently. On this side, `tests/reference/JuliaRoundTrip.wl` evaluates each exported file in plain Julia with smooth test shapes and compares it with Mathematica at a random point (relative 10⁻¹¹), then loads it with FlowRG's `load_model`; without Julia those tests print SKIP, and `runTests.wl --require-julia` turns a SKIP into a failure (lesson 16). One file per dimension, because FlowRG takes d as fixed.
:::

## Exercises

:::exercise Explicit form
Export the same flow with `"Form" -> "Explicit"` to another file and compare the two `etaMnl` definitions.
:::

:::solution
```wl
file2 = ExportJulia[adim, "Dimension" -> 2, "Form" -> "Explicit", "File" -> FileNameJoin[{$TemporaryDirectory, "lessonKPZ", "Equations2DExplicit.jl"}]];
StringCases[ReadString[file2], StartOfLine ~~ "etaMnl(" ~~ Shortest[___] ~~ "\n\n"]
```
No atoms: the shapes appear directly, and the expressions are longer. `StringCases` with `StartOfLine` and `Shortest` extracts the text from the line that starts the definition to the next blank line.
:::

:::exercise Another dimension
Export for d = 3 and compare the measure in the header of the two files.
:::

:::solution
```wl
file3 = ExportJulia[adim, "Dimension" -> 3, "File" -> FileNameJoin[{$TemporaryDirectory, "lessonKPZ", "Equations3D.jl"}]];
{JuliaParts[adim, 2]["measure"], JuliaParts[adim, 3]["measure"]}
```
aq/(2π) in d = 2 and aq²/(2π²) in d = 3: `LoopMeasure` at the given d (lesson 14).
:::

:::exercise The loop part of λ
Why is there no `betanl` function in the KPZ file?
:::

:::solution
```wl
{parts["betanl"], StringContainsQ[ReadString[file], "betanl("]}
```
λ is NonRenormalised (lesson 06): its β has no loop part, and an all-zero part is left out of the file (only `betal` is always written).
:::
