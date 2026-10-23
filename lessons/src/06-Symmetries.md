# Lesson 06 - Symmetries

**Physics step.** Decide, *before* computing any loop, what happens to every monomial (iω)^m (q²)^n of every Γ2 entry up to the truncation orders: it is a coupling of the ansatz (its flow is an equation), protected by a symmetry (its flow must vanish), neglected (its flow is computed, reported and dropped), or zero by structure. Ward identities are declared by you, not derived: the card says so.

**Stage.** Input: the Ansatz stage, tags `Protected[...]`, `Neglected[...]`, `NonRenormalised[...]`, and `"TruncationOrders"`. Output: the `"Symmetries"` stage, with one record per monomial. It is used by the projection (lesson 08) and the closure check (lesson 10).

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
ansatz = lessonUpTo["Ansatz"];
```

## 1. Using it

For KPZ: the coefficient of iω in Γ2_{h̄h} is protected by the Galilean Ward identity, the mass term by the shift symmetry h → h + c, λ is not renormalised (Galilean invariance), and every other monomial up to (iω)² (q²)¹ is neglected.

```wl
symmetries = DeclareSymmetries[ansatz,
   Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
   Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
   "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>];
```

The table has 18 rows: 3 pairs × 3 frequency powers × 2 momentum powers. Each row is a record:

```wl
Select[symmetries["Classification"], #Class === "Declared" &]
```

## 2. Tags are expressions

@region Kernel/Symmetries.wl | Neglected::usage | NonRenormalised::usage

@region Kernel/Symmetries.wl | (* Protected is a System symbol | SetAttributes[{Neglected, NonRenormalised}, Protected];

`Neglected` and `NonRenormalised` are [[Inert heads]]. `Protected` is more unusual: it is the System symbol of the attribute ([[Attributes and Protected]]). As a *function* it has no definitions, so `Protected[{hbar, h}, ...]` stays unevaluated, and the library reuses it as a tag.

```wl
{Attributes[Protected], Protected[{a, b}, "x" -> 1]}
```

:::why Declarations as data
A tag like `Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "..."]` is an inert expression: the library can read it, match it, print it on the card, and store it in the stage. Requiring a `"Reason"` makes every assumption about the physics explicit and readable later on the card, next to the PASS lines. The alternative, a set of options or flags, would carry the same information without the explanation, and would be easy to forget.
:::

## 3. Reading the tags

@quote Kernel/Symmetries.wl specMatchQ pairMatchQ tagRecord tagMatchQ

`specMatchQ` dispatches on the form of the spec ([[Several definitions of one function]]): `All`, a list of orders, or one integer. `pairMatchQ` matches a pair pattern against the unordered pair, so `{hbar, _}` matches `{h, hbar}` too. In `tagRecord`, `Association[Rest[List @@ tag]]` turns the options of a tag into an association: `List @@ tag` ([[Apply]]) gives `{pattern, opt1, opt2, ...}` and `Rest` drops the pattern.

@region Kernel/Symmetries.wl | DeclareSymmetries[ansatz_, args___] | nonren = Cases[tags, _NonRenormalised];

:::mma Multiple assignment | `{a, b} = {1, 2}` assigns element by element.
`{mf, mp} = {trunc["Frequency"], trunc["Momentum"]}` reads both truncation orders at once. Later, `{cl, reason} = Which[...]` gets a class and its reason from one [[Which]].
:::

## 4. Classifying every monomial

@region Kernel/Symmetries.wl | (* Gamma_k coefficients (regulator off) | ToString[Length[cls]] <> " monomials classified"]];

@quote Kernel/Symmetries.wl gammaCoefficient

The module rebuilds the vertices itself, silently (`Block[{Print}, BuildVertices[ansatz]]`): it needs the coefficients of Γ_k (without regulator), and taking them from its own ansatz input avoids a second stage argument that could come from a different ansatz. `gammaCoefficient` extracts the coefficient of (iν)^m (p²)^n with two nested [[Coefficient and Exponent]].

:::mma Iterators that depend on each other | `Table[expr, {pr, pairs}, {m, 0, mf}, {n, 0, mp}]` runs over every combination; `Flatten[t, 2]` removes two levels of nesting.
The result of a three-iterator `Table` is a list of lists of lists. `Flatten[..., 2]` flattens the first two levels only, leaving a flat list of records (each record is an association and must not be flattened further). `Flatten[l, n]` with a level is the safe form whenever the elements themselves are lists or associations.
```wl
{Table[{i, j}, {i, 2}, {j, 2}], Flatten[Table[{i, j}, {i, 2}, {j, 2}], 1]}
```
:::

The order of the [[Which]] is the precedence of the classes:
- **Causality**: an entry between two physical fields; Γ_φφ = 0, and its flow will be checked to vanish.
- **Structure**: an odd frequency power in a diagonal entry; Γ2_aa(ω) = Γ2_aa(−ω), so it vanishes identically and is not computed.
- **Protected** wins over everything else, but a Protected monomial that carries a coupling is a contradiction, and aborts.
- **DeclaredNumeric**: a number in the ansatz, like the 1 in front of iω; its flow must vanish.
- **Declared**: a coefficient with couplings; it gets a flow equation.
- **Neglected**; anything left is **Untagged**, and `Symmetries.allTagged` aborts.

```wl
showFailure[DeclareSymmetries[ansatz, Protected[{hbar, h}, "Frequency" -> 0, "Momentum" -> 1, "Reason" -> "a mistake"], Neglected[_, "Reason" -> "rest"],
  "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]]
```

:::why Classify before computing
The loops of lessons 07–10 are the expensive part of the pipeline. A forgotten monomial, a contradictory tag or a truncation that misses a term of the ansatz are mistakes in the *declaration*, which can be found by looking only at Γ_k. Doing it here means they abort in a fraction of a second, with a table that names the monomial, instead of after minutes of loop computation, or not at all (`docs/MODULES.md`: "the classification needs no flow, so it can fail early"). Aborting on an untagged monomial forces every truncation decision to be written down.
:::

## 5. Checks on the declaration as a whole

@region Kernel/Symmetries.wl | (* declared monomials of the ansatz beyond | AppendTo[warns, infoCheck["ASSUMPTION", "Symmetries.wardIdentities"

A truncation lower than the ansatz would silently drop a coupling's equation:

```wl
showFailure[DeclareSymmetries[ansatz, Neglected[_, "Reason" -> "rest"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 0|>]]
```

A tag that matches no monomial, because a tag of higher precedence took them all, gets a WARN. `NonRenormalised` produces an ASSUMPTION line: in v1 it cannot be checked, because that would need the flow of Γ⁽³⁾. And one ASSUMPTION line always reminds that the Protected tags are declared Ward identities.

## 6. The table

@quote Kernel/Symmetries.wl ShowClassification

The same structure as `ShowMonomials` (lesson 02): a `Grid` in a notebook, padded text in a script.

## Exercises

:::exercise A forgotten Neglected
Declare only the two Protected tags, without `Neglected`. Read the FAIL block: which monomials are untagged, and why are those of `{h, h}` not among them?
:::

:::solution
```wl
showFailure[DeclareSymmetries[ansatz, Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean"],
  Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]]
```
The six (iω)^m (q²)^n of h̄h and h̄h̄ that are neither declared nor protected. The `{h, h}` entries are classified as Causality before any tag is looked at, and the odd powers of `{hbar, hbar}` as Structure.
:::

:::exercise Precedence
Add `Neglected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "test"]` to the KPZ declaration, before the general `Neglected[_]`. What does the card say about it, and why?
:::

:::solution
```wl
DeclareSymmetries[ansatz,
   Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
   Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "h -> h + c (shift symmetry)"],
   Neglected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "test"],
   Neglected[_, "Reason" -> "beyond the ansatz"], NonRenormalised[lamk, "Reason" -> "Galilean invariance"],
   "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>];
```
A WARN `Symmetries.unusedTag`: the monomial (iω)¹(q²)⁰ of h̄h is already Protected, and Protected takes precedence over Neglected, so the new tag matches nothing.
:::

:::exercise Query the classification
List, as `{pair, m, n}`, the monomials whose flow must vanish exactly (Protected, DeclaredNumeric or Causality).
:::

:::solution
```wl
{#Pair, #Frequency, #Momentum} & /@ Select[symmetries["Classification"], MemberQ[{"Protected", "DeclaredNumeric", "Causality"}, #Class] &]
```
:::
