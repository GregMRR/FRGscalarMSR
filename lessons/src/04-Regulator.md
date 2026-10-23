# Lesson 04 - Regulator

**Physics step.** Declare the regulator R_k(q) that is added to Γ2, entry by entry, with **abstract** shape functions r(y), y = q²/k², and derive its scale derivative ∂s R = k ∂k R by the chain rule. The couplings depend on k, so ∂s R contains `dsCoupling[c]` = k dc/dk, which later becomes −η c.

**Stage.** Input: the Ansatz stage, rules `{a, b} -> R_ab`, `"Shapes" -> {...}`. Output: the `"Regulator"` stage. Regulator does not need the Vertices stage: the two are independent and meet in the Propagator (lesson 05).

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
ansatz = lessonUpTo["Ansatz"];
```

## 1. Using it

KPZ: R_{h̄h}(q) = ν q² r_n(y) and R_{h̄h̄}(q) = −2D r_d(y).

```wl
regulator = DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
```

Each entry is a **Γ2-level** entry: what is added to Γ2_ab, including the factor 2 of a diagonal pair (`docs/MODULES.md`, convention E10), so that (Γ2 + R)_{h̄h̄} = −2D(1 + r_d). The accessors give matrices in field order at q² = s:

```wl
{RegulatorMatrix[regulator, s], DsRAbstractMatrix[regulator]}
```

```wl
DsRRules[regulator, s]
```

`dsR[h, hbar]` is a placeholder for the entry of ∂s R: the loop integrands of lesson 07 are written with it, and it is replaced by the derived expression only at assembly (lesson 11), so that intermediate expressions stay small.

## 2. Arguments: rules and an option

@region Kernel/Regulator.wl | Options[DefineRegulator] | "give at least one pair {a, b} -> R_ab"];

:::mma Sequence patterns | `x__` matches a sequence of one or more arguments, `x___` of zero or more; inside the body, `{x}` is the list of them.
`DefineRegulator[ansatz_, args___]` accepts any number of arguments after the ansatz. `OptionsPattern[]` would not do here: the entries `{hbar, h} -> ...` are rules too, but their keys are lists, not option names. The constructor sorts `{args}` itself with `Cases` ([[Cases and level specifications]]): rules with a two-element list as key are entries, rules with a string key are options, and anything else is an error (`Regulator.arguments`).
:::

`Lookup[opts, "Shapes", {}]` works on a list of rules as well as on an association ([[Association]]).

## 3. Which pairs are allowed

@region Kernel/Regulator.wl | (* pairs: fields of the theory | "each pair given once"]];

`SubsetQ[fields, pair]` tests that both fields belong to the theory, and `IntersectingQ[pair, responseFields]` that the pair contains a response field. A φφ entry would give Γ2_{φφ} ≠ 0 and break causality:

```wl
showFailure[DefineRegulator[ansatz, {h, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}]]
```

Each pair is then put in field order (`SortBy[#[[1]], pos]`, [[DeleteCases, SortBy and FirstPosition]]): R is even in q and ω-independent, hence symmetric, so `{hbar, h}` and `{h, hbar}` are the same entry and may be given only once.

## 4. Your q and your k

@region Kernel/Regulator.wl | (* the user writes q (loop momentum) | "R depends on q^2 only"]];

@quote Kernel/Regulator.wl symbolsNamed toS

The entries are written with *your* symbols `q` and `k`. In these lessons they live in the notebook's own context (lesson 00), in a script in ``Global`` : the library cannot know in advance which symbol `q` is, so `symbolsNamed` finds them by their short name (`SymbolName`, [[Symbol and SymbolName]]). `k` remains the scale symbol in every later stage, and your `q` is stored as `"LoopMomentum"`.

:::mma PatternTest | `x_?test` matches x only if `test[x]` is True; `n_Integer?EvenQ` matches even integers.
`toS` replaces `q^n` with n even by `sq^(n/2)`: from here on, R is a function of q² only. An odd power of `q` survives the replacement, and `Regulator.evenInQ` catches it. `p_?test` differs from `p_ /; cond` ([[Patterns]]) only in form: the test is a function applied to the match.
:::

```wl
showFailure[DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q/k], "Shapes" -> {rn}]]
```

`sq` is a private helper symbol standing for q². It is the example of the shared-private-context gotcha in CLAUDE.md: every module's helpers live in ``FRGscalarMSR`Private`` (lesson 00, section 2), so a helper called `sq` in another module would be the same symbol.

```wl
Context[FRGscalarMSR`Private`sq]
```

## 5. What may appear in an entry

@region Kernel/Regulator.wl | allowed = Join[cps, {kk, sq}]; | "no zero entry"]];

Three separate checks: only couplings, `k`, q² and shapes as symbols (no frequency, since R is ω-independent); only the declared shapes as functions (`(h_Symbol /; ...)[___] :> h` collects the heads of all function calls); and every shape only in the form r[q²/k²], with no derivative and no other argument.

```wl
showFailure[DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2] w, "Shapes" -> {rn}]]
```

:::why Shapes stay abstract
The shape functions are never given a formula. The flow equations come out as functionals of r(y), r′(y), r″(y), ..., valid for every regulator, and the choice (Litim, exponential, ...) is made in FlowRG, after the export (lesson 15). Abstract shapes also make the symbolic checks stronger: `ExactlyZeroQ` gives every atom r(y), r′(y) its own random value (lesson 00, section 5), so an identity that passes holds for arbitrary shapes, not by luck for one of them. The price is that the shapes must appear in one fixed form, r[q²/k²], which is what the checks of this section enforce.
:::

## 6. The scale derivative, two ways

@region Kernel/Regulator.wl | (* d_s R by the chain rule | "d_s R = sum_c dR/dc dsCoupling[c]

:::mma D and unknown functions | `D[e, x]` differentiates; on a function with no definition it applies the chain rule and writes the derivative as `Derivative[1][r][arg]`, displayed `r'[arg]`.
`k D[rn[q^2/k^2], k]` gives −2 q²/k² r_n′(q²/k²), that is −2 y r′(y): Mathematica knows the chain rule without knowing r. `Derivative[1][rn]` is itself an expression, the derivative of the function `rn`, which is why `ExactlyZeroQ` treats `Derivative[1][rn][y]` as one atom.
```wl
k D[rn[q^2/k^2], k]
```
:::

:::mma Map at a level | `Map[f, matrix, {2}]` applies f to every element at depth 2, that is, to every entry of a matrix.
`Map[f, {{1, 2}, {3, 4}}, {2}]` is `{{f[1], f[2]}, {f[3], f[4]}}`. Without the level specification, `Map` would apply f to the rows.
:::

The first route differentiates at fixed couplings and adds Σ_c ∂R/∂c · `dsCoupling[c]`. The second route makes each coupling an explicit function `cf[c][k]`, takes the total derivative k d/dk, and replaces `cf[c]'[k]` by `dsCoupling[c]/k`. `Regulator.dsRChainRule` requires the two to agree, with `ExactlyZeroQ`.

## 7. The stage stores functions

@region Kernel/Regulator.wl | stage = stampStage[<|"Stage" -> "Regulator" | "Checks" -> checks, "Warnings" -> warns

@quote Kernel/Regulator.wl RegulatorMatrix DsRegulatorMatrix DsRAbstractMatrix DsRRules

:::mma Function and Evaluate | `Function[{x}, body]` holds its body unevaluated; `Function[{x}, Evaluate[body]]` evaluates it first, so the stored function contains the computed result.
`"RMatrix" -> Function[Evaluate[{sq}], Evaluate[rmat]]` stores the matrix as a function of q²: `regulator["RMatrix"][s]` substitutes `s` for `sq`. `Function` has the attribute `HoldAll` ([[Hold attributes]]); without `Evaluate`, the body would be the *symbol* `rmat`, a `Module` variable that is gone after the call. The first `Evaluate` is needed for the same reason: `{sq}` must be the list of the actual private symbol.
```wl
{Function[{x}, Expand[(x + 1)^2]], Function[{x}, Evaluate[Expand[(x + 1)^2]]], regulator["RMatrix"]}
```
:::

`DsRRules` pairs each abstract `dsR[a, b]` with its derived entry (`MapThread[Rule, {m1, m2}, 2]` builds rules entry by entry) and removes the zero entries and the duplicates (`DeleteDuplicates`: R is symmetric, so `dsR[h, hbar]` appears twice).

:::why Matrices as functions of q²
The later modules need R and ∂s R at the loop momentum q, and at shifted momenta such as (q + p)² for the momentum projection (lesson 08). Storing a function of q² lets each caller pass whatever it needs, `scal[q, q]` or `scal[q + p, q + p]`, without the regulator knowing about it. Since R depends on q only through q² (section 4), one argument is enough.
:::

## Exercises

:::exercise ∂s R by hand
Compute k ∂k of R_{h̄h} = ν(k) q² r_n(q²/k²) by hand and compare with `DsRRules`.
:::

:::solution
k ∂k gives (k ν′) q² r_n + ν q² (−2 q²/k²) r_n′ = q² `dsCoupling[nuk]` r_n(y) − 2 ν q⁴/k² r_n′(y), the first rule of:
```wl
DsRRules[regulator, s]
```
:::

:::exercise An unused shape
Regulate only the h̄h entry but declare both `rn` and `rd`. What does the card say? Is the result usable?
:::

:::solution
```wl
DefineRegulator[ansatz, {hbar, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn, rd}];
```
A WARN `Regulator.unusedShapes`; the stage is valid (a regulator need not touch every entry), the warning only flags a likely typo.
:::

:::exercise Recognition by name
Write the regulator with the momentum called `p` instead of `q`. Which check fails, and why?
:::

:::solution
```wl
showFailure[DefineRegulator[ansatz, {hbar, h} -> nuk p^2 rn[p^2/k^2], "Shapes" -> {rn}]]
```
`symbolsNamed` finds no symbol called `q`, so `p` is just an unknown symbol: `Regulator.knownSymbols`. The loop momentum must be called `q` (it is a reserved name, lesson 00, section 6).
:::
