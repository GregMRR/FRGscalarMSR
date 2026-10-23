# Lesson 01 - Theory

**Physics step.** Declare the N physical scalars φ_a, create their MSR response fields φ̄_a, fix the field order (φ1, φ̄1, φ2, φ̄2, ...) that indexes every matrix downstream (Γ2, R, G), and fix the space dimension d, which stays a symbol.

**Stage.** Input: two options. Output: the `"Theory"` stage. It is the root of the pipeline: its `"Upstream"` list is empty.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]]
```

## 1. Using it

KPZ has one scalar, the height `h`:

```wl
theory = DefineTheory["Scalars" -> {h}, "Dimension" -> Dim]
```

The card has two PASS lines, one per `assertCheck` in the code we read below. The returned stage is an [[Association]]:

```wl
{Keys[theory], theory["Fields"], theory["ResponseOf"], Context[hbar]}
```

`hbar` did not exist before the call: `DefineTheory` created it, in the same context as `h`. With two scalars the field order interleaves fields and response fields:

```wl
theory2 = DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim];
theory2["Fields"]
```

## 2. The constructor, line by line

@quote Kernel/Theory.wl DefineTheory

Every constructor in the library has this skeleton: `Block[{$FRGModule = "..."}, Module[{t0 = AbsoluteTime[], ...}, checks...; stage = stampStage[<|...|>]; Print...; PrintCard[stage]; stage]]`. The `Block` names the module in FAIL blocks (see [[Block and Module]]), the `Module` holds the local variables, `AbsoluteTime[]` (the current time in seconds) measures the `"Timing"` shown on the card, and the stage is the value of the last expression.

### Options, and why they are checked

The function takes only options: `opts : OptionsPattern[]` names the whole sequence of options `opts`, so that `{opts}` is the list of rules that were passed (see [[Options]] and [[Alternatives and named patterns]]).

:::mma MatchQ and repeated patterns | `MatchQ[e, patt]` is True if e matches patt; `p ..` matches one or more repetitions of p, `p ...` zero or more.
`MatchQ[{opts}, {(Rule | RuleDelayed)["Scalars" | "Dimension", _] ...}]` reads: a list of zero or more rules (written with `->` or `:>`) whose key is `"Scalars"` or `"Dimension"`. `Rule[a, b]` is the full form of `a -> b`: a pattern can match on any head, including `Rule`.
```wl
{MatchQ[{1, 2, 3}, {_Integer ..}], MatchQ[{}, {_Integer ..}], MatchQ[{}, {_Integer ...}], FullForm["a" -> 1]}
```
:::

Without this check, `"Scalar" -> {h}` (a typo) would be silently ignored: `OptionValue["Scalars"]` would return the default `{}`, and the error would surface later with a message about an empty list. The check names the real problem:

```wl
showFailure[DefineTheory["Scalar" -> {h}, "Dimension" -> Dim]]
```

### Checks as a list of records

`checks = {assertCheck[...]}` starts the list of card lines and `AppendTo[checks, ...]` extends it. Some checks are plain `FRGAssert` calls instead: they abort on failure like the others, but do not add a PASS line to the card. The rule in the library is that a check worth reading on the card is an `assertCheck`; a precondition whose success is obvious is an `FRGAssert`.

:::mma AppendTo | `AppendTo[x, e]` is `x = Append[x, e]`: it rebinds the variable x to a longer list.
Lists, like associations, are values: `Append` returns a new list. `AppendTo` assigns it back to the variable, which works only on a variable (here a `Module` local).
:::

`Theory.N` allows one or two scalars. With two, the propagator denominator Δ(ω) is a polynomial of degree 2 in ω, the largest degree v1 integrates in closed form (lesson 09).

```wl
showFailure[DefineTheory["Scalars" -> {u1, u2, u3}, "Dimension" -> Dim]]
```

### The symbols must be fresh

`assertFreshSymbols` (lesson 00, section 6) runs on the scalars, then on the generated response fields, then on the dimension symbol. Note what happens with a symbol that has a value:

```wl
uu = 1;
showFailure[DefineTheory["Scalars" -> {uu}, "Dimension" -> Dim]]
```

```wl
Clear[uu]
```

The offending expression is `{1}`, not `{uu}`: Mathematica evaluated the arguments before `DefineTheory` saw them, so the library never knew the name `uu`. That is why `assertFreshSymbols` first checks that every element is still a `Symbol` (the `.symbols` check), and only then looks at definitions (`.values`), which catch symbols that are functions (`uu[x_] := ...`) and therefore do not evaluate away.

### Making the response fields

:::mma Symbol and SymbolName | `Symbol["name"]` returns the symbol with that name, creating it if needed; `SymbolName[s]` and `Context[s]` give its short name and its context.
`Symbol[Context[#] <> SymbolName[#] <> "bar"] & /@ phis` turns `{h}` into `{hbar}`. Passing the full name `` "Global`hbar" `` puts the new symbol in the same context as `h`, whatever `$Context` is when the library runs. `Context` has the attribute `HoldFirst` ([[Hold attributes]]): it looks at the expression as written. `Context[#] &` works because the pure function puts `h` itself in place of `#` before `Context` runs; on a computed expression, write `Context @@ {expr}` (CLAUDE.md, Gotchas).
```wl
{Attributes[Context], Context @@ {First[{h}]}, Quiet[Check[Context[First[{h}]], "error: Context saw First[{h}], not h"]]}
```
:::

:::mma Apply | `f @@ expr` replaces the head of expr by f: `f @@ {a, b}` is `f[a, b]`.
`Context @@ {First[{h}]}` first evaluates the list to `{h}`, then turns it into `Context[h]`. `Apply` is used throughout the library to change heads: `List @@ (a + b)` is `{a, b}`, `Plus @@ {a, b}` is `a + b`.
:::

:::mma Transpose and Flatten | `Transpose[{l1, l2}]` pairs up two lists element by element; `Flatten` removes the inner brackets.
`Flatten[Transpose[{{H, Psi}, {Hbar, Psibar}}]]` gives `{{H, Hbar}, {Psi, Psibar}}` and then `{H, Hbar, Psi, Psibar}`: the field order.
:::

:::why Response fields are created, not declared
You declare `h` and the library makes `hbar`. You cannot forget a response field, give it an unexpected name, or pair it with the wrong scalar, and every later module can find the partner of a field through `"ResponseOf"`. The naming rule `<name>bar` is part of the interface: you write `hbar` in the ansatz, which is why `assertFreshSymbols` also checks the generated names (a coupling you already called `hbar` would be caught here, not later).
:::

:::why The field order is fixed once, here
Every matrix of the pipeline (Γ2, the regulator R, the propagator G) is indexed by the field order. Fixing it in the root stage, and passing the stage along, means no module can choose a different order. The order interleaves each field with its response field because the MSR structure is then visible as a block structure (causality forces Γ2_{φφ} = 0, `docs/CONVENTIONS.md`).
:::

### The dimension

The dimension must be either a fresh symbol or `FixedDimension[2]`.

:::mma SameQ and Equal | `a === b` (SameQ) asks whether two expressions are identical and is always True or False; `a == b` (Equal) asks for mathematical equality and may stay undecided.
`dim === FixedDimension[2]` compares the structure of the expressions. `==` would be wrong for this: it answers questions about numbers, and on symbols it often returns itself unevaluated. `2 == 2.` is True (equal numbers), `2 === 2.` is False (an integer and a machine real are different expressions).
```wl
{2 == 2., 2 === 2., nuk == 2, nuk === 2}
```
:::

:::mma Inert heads | A symbol with no definitions, used as a head, is just a label: `FixedDimension[2]` evaluates to itself and carries its meaning as data.
`FixedDimension` has no rules at all; `DefineTheory` recognises it with `===`. The attribute `Protected` ([[Attributes and Protected]]) only stops anyone from giving it rules. The library uses inert heads everywhere to represent things that are not to be computed by Mathematica: `dt`, `grad`, `FieldFactor`, `scal` before its rules (lesson 02).
:::

```wl
theory2D = DefineTheory["Scalars" -> {A1, A2}, "Dimension" -> FixedDimension[2]];
theory2D["Dimension"]
```

A plain `2` is refused:

```wl
showFailure[DefineTheory["Scalars" -> {h}, "Dimension" -> 2]]
```

:::why The dimension stays a symbol
The flow equations are derived for arbitrary d, as rational functions of the symbol `Dim`, and exported to FlowRG with d as a parameter. Fixing d early would lose that generality without anyone noticing, so a plain number is refused: `"Dimension" -> 2` and a `Dim` that got the value 2 are indistinguishable once evaluated, and both abort. The only exception is the Levi-Civita operator `jac`, which exists only in d = 2; there you declare it explicitly with the inert head `FixedDimension[2]`, and the card says "2 (declared)".
:::

### The stage

The stage is assembled with `stampStage` ([[Hash]]). `AssociationThread[phis -> bars]` builds the map `<|h -> hbar|>`. `"Pipeline" -> "DerivativeExpansion"` names the only pipeline v1 implements (`docs/FUTURE.md` lists others).

## 3. ResponseField

@quote Kernel/Theory.wl ResponseField

A small accessor with its own check: asking for the response field of something that is not a physical scalar aborts.

```wl
{ResponseField[theory, h], showFailure[ResponseField[theory, hbar]]}
```

## 4. The running example from now on

From the next lesson on, each lesson starts by rebuilding the KPZ example up to the stage it needs, with the cards hidden:

@quote lessons/LessonSetup.wl $lessonStages $lessonSteps lessonUpTo

`Hold[expr]` keeps an expression unevaluated ([[Hold attributes]]); `ReleaseHold` evaluates it. The chain is stored as a list of held assignments, and `lessonUpTo` evaluates the first n of them. `Block[{Print}, ...]` gives `Print` no definitions for the duration, so it prints nothing: the cards are hidden, but a failing check still aborts.

```wl
kpz = lessonUpTo["Theory"];
kpz["Fields"]
```

## Exercises

:::exercise A reserved name
Try to declare a scalar called `k`. Which check fails, and why is `k` reserved?
:::

:::solution
```wl
showFailure[DefineTheory["Scalars" -> {k}, "Dimension" -> Dim]]
```
`k` is the RG scale; the library writes it into every flow equation (lesson 00, section 6).
:::

:::exercise A response field that already exists
Give `Phibar` a value, then declare the scalar `Phi`. Which check fails? Clean up afterwards.
:::

:::solution
```wl
Phibar = 0;
showFailure[DefineTheory["Scalars" -> {Phi}, "Dimension" -> Dim]]
```
```wl
Clear[Phibar]
```
`Theory.responseFields.symbols`, with offending value `{0}`: the generated name is checked like a user symbol. `Symbol["...Phibar"]` returns the existing symbol `Phibar`, which then evaluates to its value 0, exactly as `uu` did in section 2.
:::

:::exercise Read the stage
Using only `theory2` and the functions of this lesson, build the list of pairs `{{H, Hbar}, {Psi, Psibar}}` in at least two different ways.
:::

:::solution
```wl
{Transpose[{theory2["Scalars"], theory2["ResponseFields"]}], Normal[theory2["ResponseOf"]] /. Rule -> List, {#, ResponseField[theory2, #]} & /@ theory2["Scalars"]}
```
The second way uses `Normal`, which turns an association into a list of rules, then replaces the head `Rule` by `List` (`a -> b` is `Rule[a, b]`).
:::
