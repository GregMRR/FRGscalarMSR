# Lesson 02 - Ansatz

**Physics step.** Γ_k = ∫dt d^dx 𝓛. The ansatz is the density 𝓛, written by you with the operators `dt`, `grad`, `div`, `lap` (and `jac` in d = 2). The module rewrites it as a list of **monomials**: a coefficient (couplings only) times a product of derivatives of single fields. That is the form the next module, Vertices, can Fourier-transform term by term (∂t → +iω, ∇ → −iq, `docs/CONVENTIONS.md`).

**Stage.** Input: the Theory stage, `"Density" -> 𝓛`, `"Couplings" -> {...}`. Output: the `"Ansatz"` stage, whose key `"Monomials"` is the parsed list.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
theory = lessonUpTo["Theory"]
```

## 1. Using it

The KPZ ansatz, Γ_k = ∫ h̄ (∂t h − ν ∇²h − λ/2 (∇h)²) − D h̄²:

```wl
ansatz = DefineAnsatz[theory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2, "Couplings" -> {nuk, lamk, Dk}];
```

The table is printed before the card so you can check the parse by eye. Each row is one record of `ansatz["Monomials"]`:

```wl
ansatz["Monomials"][[4]]
```

A factor `FieldFactor[φ, nt, {indices}]` is ∂t^nt ∂_i ∂_j ... φ. Repeated `SpaceIndex[n]` are summed: `FieldFactor[h, 0, {SpaceIndex[1], SpaceIndex[1]}]` is ∇²h, and two factors carrying `SpaceIndex[1]` once each are ∇h·∇h.

## 2. A small language of inert operators

@region Kernel/Ansatz.wl | dt::usage | lap::usage

@region Kernel/Ansatz.wl | SetAttributes[{dt, grad | systemOps =

`dt`, `grad`, `div`, `lap`, `jac` have usage messages and no definitions: they are [[Inert heads]], protected so nobody gives them one. `opHeads` and `systemOps` store **patterns** in variables, built with [[Alternatives and named patterns]]: `FreeQ[e, opHeads]` reads "e contains none of the five operators".

:::why The ansatz has its own operators
Mathematica's `Dt` is the total derivative and knows nothing about which symbols are couplings: it differentiates them too. Its `Laplacian` needs explicit coordinates and, applied to a symbol that does not depend on them, gives zero before the library can see anything (`docs/MODULES.md`, Ansatz caveat). With inert operators, the library alone decides what a derivative means: it acts on fields, never on couplings, which depend on k only.
```wl
{Dt[nuk h], Laplacian[h, {x, y}]}
```
:::

:::mma FullForm | `FullForm[e]` shows the tree Mathematica actually stores, with every head written out.
What you type is notation. `grad[h].grad[h]` is stored as `Dot[grad[h], grad[h]]`, `a - b` as `Plus[a, Times[-1, b]]`, `x^2` as `Power[x, 2]`. Patterns match the stored tree, so this is the form to think in when you read the parser.
```wl
{FullForm[grad[h].grad[h]], FullForm[hbar (dt[h] - nuk lap[h])]}
```
:::

## 3. Checks before any work

@region Kernel/Ansatz.wl | DefineAnsatz[theory_, opts | "a coupling has the name of a field or of the dimension"];

`validateStage` comes first ([[Hash]]): the theory must be a genuine, unmodified Theory stage. Then the options, then the couplings ([[OwnValues and DownValues]] through `assertFreshSymbols`).

@region Kernel/Ansatz.wl | (* E7: System derivative heads | "jac with d = 2 (declared)"]]];

:::mma Complement and Intersection | `Complement[all, l1, l2, ...]` gives the elements of all that are in none of the li; `Intersection[l1, l2]` the common ones.
`badSyms` collects every non-System symbol of the density (heads included: `Heads -> True`, see [[Cases and level specifications]]) and removes the fields, the couplings, `Dim` and the operators. Anything left is a symbol you forgot to declare.
:::

These checks look at the raw input, so their hints can be precise:

```wl
showFailure[DefineAnsatz[theory, "Density" -> hbar dt[nuk h], "Couplings" -> {nuk}]]
```

```wl
showFailure[DefineAnsatz[theory, "Density" -> hbar Dt[h] + g hbar^2, "Couplings" -> {}]]
```

:::why Validate the input before transforming it
A malformed density caught after the expansion would show up as a strange monomial, far from what you wrote. Checked first, it is reported as you wrote it, with a hint in the vocabulary of the input (`write nk dt[H], not dt[nk H]`). Each check also documents an assumption of the module: couplings depend on k only, the density is polynomial, only package operators are allowed.
:::

## 4. Vector rank: is this a scalar or a vector?

Before expanding, the module checks the vector algebra: `grad` makes a vector, `div` and `.` make scalars, a sum cannot mix the two.

@quote Kernel/Ansatz.wl rank rankOf needRank

:::mma Which | `Which[c1, v1, c2, v2, ..., True, v]` returns the value of the first condition that is True.
A chain of `If`s written flat. The final `True` is the "otherwise" case.
:::

:::mma Several definitions of one function | A function can have many rules; Mathematica tries the more specific patterns first and the most general last.
`rankOf` has one rule per kind of expression: `grad[a_]`, `dt[a_]`, `e_Plus`, `e_Times`, `Power[a_, n_]`, and finally `rankOf[e_, fld_]`, which matches anything. Mathematica orders the rules by specificity, so the catch-all only runs when nothing else matched. Here it is a failed check: an expression the parser does not know is an error, not something to guess about. `DownValues` ([[OwnValues and DownValues]]) lists the rules in the order they are tried; the catch-all is last. Each left-hand side stays wrapped in `HoldPattern`, which stops it from being evaluated (evaluating `rankOf[e_, fld_]` would run the catch-all). This "dispatch on the head" is the main way the library structures its code; you will see it in every module.
```wl
{Length[DownValues[FRGscalarMSR`Private`rankOf]], DownValues[FRGscalarMSR`Private`rankOf][[-1, 1]]}
```
:::

`rankOf[e_Plus, ...]` computes the rank of every term (`List @@ e`, [[Apply]]) and requires them all equal (`Equal @@ r` is `Equal[r1, r2, ...]`).

```wl
{FRGscalarMSR`Private`rank[grad[h].grad[h], h | hbar], FRGscalarMSR`Private`rank[grad[h], h | hbar]}
```

A product of two different gradients without a dot is refused:

```wl
showFailure[FRGscalarMSR`Private`rank[grad[h] grad[hbar], h | hbar]]
```

But try the same with the same gradient twice:

```wl
showFailure[FRGscalarMSR`Private`rank[grad[h] grad[h], h | hbar]]
```

:::mma Automatic simplification of Times and Plus | `Times` and `Plus` sort their arguments and combine equal ones as soon as an expression is built: `x x` becomes `x^2`, `a + a` becomes `2 a`.
`Times` has the attributes `Flat`, `Orderless` and `OneIdentity`, and built-in rules that collect powers. So `grad[h] grad[h]` is already `Power[grad[h], 2]` when `rank` sees it, and the failure comes from the `Power` rule ("only scalars can be raised to a power"), not from the product rule. Code that walks expressions must be written for the stored form, never for the typed one.
```wl
{FullForm[grad[h] grad[h]], Attributes[Times]}
```
:::

## 5. Expanding the operators: the Leibniz rule

Next, every operator is pushed down until it acts on a single field, and every field becomes a `FieldFactor`.

@quote Kernel/Ansatz.wl derive toFactors

:::mma Curried definitions | `f[a_][x_] := ...` defines a function of a whose result is itself a function of x.
`derive[op]` is "the derivation in direction op": `derive["t"]` is ∂t, `derive[SpaceIndex[...]]` is ∂_i. Each rule implements one property of a derivation: linearity (`e_Plus`), the Leibniz rule (`e_Times`), the chain rule for powers, the action on a single field (one more `dt`, or one more index), and zero on constants. Such rules are stored in `SubValues[derive]` rather than `DownValues`.
:::

:::mma MapAt, Table and Sum | `MapAt[f, e, i]` applies f to the i-th part of e only; `Table[expr, {i, n}]` and `Sum[expr, {i, n}]` build a list or a sum for i = 1..n.
`derive[op][e_Times] := Sum[MapAt[derive[op], e, i], {i, Length[e]}]` is the Leibniz rule written literally: the sum over i of the product with only the i-th factor differentiated. `{i, n}` is an **iterator**; `Table`, `Sum`, `Do` all take them.
```wl
FRGscalarMSR`Private`derive["t"][FieldFactor[h, 0, {}] FieldFactor[hbar, 0, {}]]
```
:::

`toFactors` defines two local functions, `S` (convert a scalar expression) and `V` (convert a vector expression, given the index it carries), with a rule per operator: `lap[x]` is two derivatives along one fresh index, `div[v]` contracts a fresh index of `v`, a dot product gives both vectors the same fresh index.

:::mma Local functions | Symbols created by `Module` can be given definitions like any symbol; they are private to that call and can use the other local variables.
`S` and `V` exist only inside one call of `toFactors`, and they see `fld` (the pattern of the fields) and `new` (the index generator) without these being passed around. `S[x_] /; FreeQ[x, fld | opHeads] := x` is a conditional definition ([[Patterns]]): a field-free factor is left as it is. As with [[Several definitions of one function]], the final `S[x_]` is the catch-all failure.
:::

:::mma With | `With[{j = value}, body]` evaluates value once and inserts it into body as a constant.
`With[{j = new[]}, Expand[derive[j][derive[j][S[x]]]]]` draws one fresh index and uses the *same* one for both derivatives: that is the Laplacian ∂_j∂_j. Writing `new[]` twice would give two different indices. `With` is the lightest of the three scoping constructs ([[Block and Module]]): it creates no variable at all.
:::

`new[] := SpaceIndex[{"dummy", ++counter}]` (in `DefineAnsatz`) increments the local `counter` (`++counter` adds 1 and returns the new value), so every call gives an index never used before. `Expand` multiplies out products of sums.

```wl
Module[{c = 0}, FRGscalarMSR`Private`toFactors[lap[h hbar], h | hbar, SpaceIndex[{"dummy", ++c}] &]]
```

This is ∇²(h h̄) = h̄ ∇²h + 2 ∇h·∇h̄ + h ∇²h̄. The vector case of a product, `V[x_Times, i_]`, uses `Pick[fs, r, 0]`, which keeps the factors whose rank is 0: in `H grad[Psi]` the scalar `H` is converted with `S` and the vector `grad[Psi]` with `V`.

## 6. Canonical monomials

After expansion the density is a sum of products of `FieldFactor`s with dummy indices. Two monomials that differ only by the names of their dummy indices, or by the order of their factors, are the same term and must be combined.

@region Kernel/Ansatz.wl | FRGAssert["Ansatz.scalarDensity" | mons = splitTerm /@ terms;

@quote Kernel/Ansatz.wl splitTerm canonMinForm canonical

`splitTerm` separates each term into its coefficient (the factors free of `FieldFactor`) and its list of field factors, writing a power `f^2` as two copies of `f`.

:::mma Sort and canonical order | `Sort[list]` puts any expressions in Mathematica's canonical order; `Permutations[list]` lists all orderings.
`Sort` works on arbitrary expressions, not only numbers: every expression has a place in one total order. `canonMinForm` tries every numbering of the summed indices (`Permutations`), sorts the indices of each factor (derivatives commute) and then the factors, and keeps the smallest result (`First[Sort[...]]`). Equal monomials give the same minimum, whatever their dummy names. `Function[perm, ...]` is a [[Pure functions]] with a named argument.
```wl
{Permutations[{1, 2, 3}], Sort[{FieldFactor[hbar, 0, {}], FieldFactor[h, 1, {}], FieldFactor[h, 0, {}]}]}
```
:::

:::mma Except | `Except[p]` matches anything that does not match p.
`Cases[..., Except[_EpsIndex]]` collects the summed indices but not the slots of ε.
:::

:::why Canonical form by brute force
Trying every permutation of the dummy indices costs n! steps, which looks wasteful. But a monomial of the derivative expansion carries few contractions (one at order p², two at p⁴), and the result is a class invariant: two monomials are equal exactly when their minima are, for any number of contractions (review E9 of `design/MOMENTUMPROJECTION.md`). A cleverer canonicalisation would be faster and much harder to prove right. Because of it, `lap[lap[h]]` and `div[grad[lap[h]]]` are recognised as the same term (exercise 2).
:::

For `jac`, the two slots of ε_ij are swapped as well: a monomial equal to its own swap is zero, and when the swapped form is the smaller one the sign is −1.

## 7. Combining, recording, sorting

@region Kernel/Ansatz.wl | mons = Function[m, With[{c = canonical | mons = SortBy[mons,

:::mma GroupBy and KeyValueMap | `GroupBy[list, key -> value]` collects the elements with the same key into an association of lists; `KeyValueMap[f, assoc]` applies f[key, value] to each entry.
`GroupBy[mons, Last -> First]` groups the pairs `{coefficient, factors}` by their factors (the last part) and keeps the coefficients (the first part). `KeyValueMap[{Together[Total[#2]], #1} &, ...]` adds the coefficients of each group. Equal monomials are combined in one pass, without comparing them pairwise.
```wl
GroupBy[{{1, x}, {2, y}, {3, x}}, Last -> First]
```
:::

:::mma DeleteCases, SortBy and FirstPosition | `DeleteCases[list, patt]` drops the matching elements; `SortBy[list, f]` sorts by the value of f (a list of criteria breaks ties in turn); `FirstPosition[list, x]` is the position of the first x.
`DeleteCases[mons, {0, _}]` removes the monomials whose coefficients cancelled; their number is reported as an INFO line. `SortBy` then orders the monomials by number of fields, then by the fields, then by time order, then by momentum order: the order of the printed table.
:::

`monomialRecord` turns each pair into the record you saw in section 1:

@quote Kernel/Ansatz.wl monomialRecord

:::why One simple record per monomial
Every later module reads `ansatz["Monomials"]`, never the density you typed. The records state explicitly what they would otherwise recompute: the fields, how many are response fields, the time and momentum orders. The parser is the only code that has to understand the ansatz language; the rest of the pipeline works on a flat list of records. The printed table is the same records shown to you, so the step where your physics turns into data can be checked by eye.
:::

## 8. Checks on the monomials, and the stage

@region Kernel/Ansatz.wl | epsQ = ! FreeQ | declared but absent from the density:

These are the physical assumptions of v1: every monomial contains a response field (MSR structure), there are no terms linear in the fields, at most quartic monomials, and a quadratic `jac` term, a total derivative, is refused rather than silently doing nothing. `Count[list, patt]` counts the elements that match ([[Patterns]]).

```wl
showFailure[DefineAnsatz[theory, "Density" -> hbar dt[h] + nuk h lap[h], "Couplings" -> {nuk}]]
```

@region Kernel/Ansatz.wl | stage = stampStage[<|"Stage" -> "Ansatz" | PrintCard[stage];

:::mma Sequence | `Sequence[a, b]` is spliced into whatever contains it: `f[x, Sequence[a, b]]` is `f[x, a, b]`; `Sequence[]` disappears.
`Sequence @@ If[epsQ, {"Epsilon" -> True}, {}]` adds the key `"Epsilon"` only when a monomial carries ε. Without ε the key is absent, and later modules read it with `Lookup[ansatz, "Epsilon", False]`.
:::

The Ansatz stage contains the whole Theory stage under `"Theory"`, and `"Upstream"` lists its fingerprint. Every stage carries its ancestors, so a later stage is self-contained.

## 9. Printing the table

@quote Kernel/Ansatz.wl ShowMonomials fieldString MonomialString

`MonomialString` writes a record back in the input language. `Nest[f, x, n]` applies f n times (`Nest["dt[" <> # <> "]" &, "h", 2]` is `"dt[dt[h]]"`), `StringRiffle[list, sep]` joins strings with a separator, `Tally` counts repetitions, `MapThread[f, {l1, l2}]` applies f to pairs of elements. In a notebook the table is a `Grid`; in a script it is padded text, chosen with the same `$FrontEnd === Null` test as the FAIL block ([[Abort, CheckAbort, Exit]]).

```wl
MonomialString /@ ansatz["Monomials"]
```

## 10. jac, in d = 2

`jac[a, b]` = ε_ij ∂_i a ∂_j b exists only in d = 2, so its theory must declare `FixedDimension[2]` (lesson 01). The slots of ε are tagged `EpsIndex[1]` and `EpsIndex[2]`:

```wl
theory2D = DefineTheory["Scalars" -> {A1, A2}, "Dimension" -> FixedDimension[2]];
ansatz2D = DefineAnsatz[theory2D, "Density" -> A1bar dt[A1] + A2bar dt[A2] + g A1bar jac[A1, A2] + g2 A1bar jac[A2, A1], "Couplings" -> {g, g2}];
Last[ansatz2D["Monomials"]]
```

`jac[A2, A1]` is `−jac[A1, A2]`: the two terms combine into one monomial with coefficient `g - g2`. A quadratic `jac` is a total derivative and aborts:

```wl
showFailure[DefineAnsatz[theory2D, "Density" -> A1bar dt[A1] + g jac[A1bar, A1], "Couplings" -> {g}]]
```

The momentum projection of `jac` terms needs the d = 2 identity of `wedge` (lesson 00, section 7); it reappears in lesson 08.

## Exercises

:::exercise Predict the parse
Write down, by hand, the monomials of `hbar div[h grad[h]]`. Then check with `DefineAnsatz` (no couplings are needed).
:::

:::solution
∇·(h ∇h) = ∇h·∇h + h ∇²h, so two cubic monomials, each with coefficient 1:
```wl
DefineAnsatz[theory, "Density" -> hbar div[h grad[h]], "Couplings" -> {}];
```
:::

:::exercise Two ways to write ∇⁴
Give `lap[lap[h]]` and `div[grad[lap[h]]]` two different couplings in one density. How many monomials do you expect, and with which coefficient?
:::

:::solution
```wl
DefineAnsatz[theory, "Density" -> c1 hbar lap[lap[h]] + c2 hbar div[grad[lap[h]]], "Couplings" -> {c1, c2}];
```
One monomial, `hbar lap[lap[h]]` with coefficient `c1 + c2`, of order p⁴: the canonical form of section 6 recognises the two as the same term.
:::

:::exercise Query the records
From `ansatz`, select the monomials of momentum order 2 and print them in the input language with their coefficients.
:::

:::solution
```wl
{#Coefficient, MonomialString[#]} & /@ Select[ansatz["Monomials"], #MomentumOrder == 2 &]
```
`#MomentumOrder` reads a key of the record ([[Association]]); `Select` keeps the records that pass the test.
:::
