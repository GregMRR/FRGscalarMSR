# Lesson 13 - Dimensions

**Physics step.** Derive the scaling dimension of every coupling, every field and of time from the requirement that Γ_k = ∫dt dᵈx 𝓛 is dimensionless, with [x] = 1/k. The dimensions are written in a basis you choose: the RG scale k and N + 1 **reference couplings** (N = number of scalars), whose running defines the anomalous dimensions, ∂s ref = −η ref. This module needs only the ansatz, so it can run right after lesson 02, before any loop.

**Stage.** Input: the Ansatz stage and `"References"`. Output: the `"Dimensions"` stage, used by Adimensional (lesson 14).

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
ansatz = lessonUpTo["Ansatz"];
```

## 1. Using it

KPZ has one scalar, so two references: ν and D.

```wl
dims = DeriveDimensions[ansatz, "References" -> {nuk -> "nu", Dk -> "D"}, "Scale" -> k];
```

The table gives each coupling's dimension, its dimensionless name and the linear part of its flow. λ has dimension k^{1−d/2} ν^{3/2}/√D, so the dimensionless coupling is `alamk` = λ/(k^{1−d/2} ν^{3/2} D^{−1/2}), whose flow starts with −1 + d/2 − η_D/2 + 3η_ν/2. The frequency scales as k² ν, so z = 2 − η_ν.

```wl
{DimensionOf[dims, lamk], DimensionOf[dims, h], DimensionOf[dims, "Frequency"], LinearPart[dims, lamk]}
```

:::why Why "Scale" -> k is passed explicitly here
The option's default is `"Scale" :> Symbol["Global`k"]`, the symbol `k` of the ``Global`` context. In a plain notebook that is the `k` you typed in the regulator. In these lesson notebooks, each with its own context (lesson 00), your `k` is a different symbol, so the default would name the wrong one; passing `"Scale" -> k` names yours. The delayed rule `:>` in `Options` makes the default evaluate only when it is used ([[Rules and replacement]]). The same assumption appears twice more in the library: the dimensionless loop momentum `aq` is created as ``Global`aq`` (Dimensions, Adimensional), and Assemble falls back to ``Global`q`` when the regulator has no momentum.
:::

## 2. One linear equation per monomial term

@quote Kernel/Dimensions.wl basisForm linearPart couplingPowers termsOf monomialEquations equationRow

A dimension is an **exponent vector** over the basis (k, ref₁, ..., ref_{N+1}): k^x ν^{y₁} D^{y₂} is the vector {x, y₁, y₂}. For a monomial term c × fields × derivatives, "Γ dimensionless" is the linear equation

[c] + Σ[fields] − (number of ∂t)·[t] + (number of spatial indices)·e_k + [t] − d·e_k = 0,

where the last two terms are the measure dt dᵈx. `couplingPowers` reads the power of each coupling in a term as c ∂term/∂c / term, a number for a monomial.

```wl
dims["Equations"][[1]]
```

## 3. Solving the system

@region Kernel/Dimensions.wl | basis = Prepend[refs, kk]; nb = Length[basis]; | "references independent"]];

:::mma Linear algebra | `NullSpace[m]` gives a basis of the solutions of m.x = 0; `MatrixRank[m]` the rank; `LinearSolve[m, b]` solves m.x = b.
The unknowns are [t], the field dimensions and the non-reference couplings. Before solving, the module counts the **freedoms**: with the references treated as unknowns too, the null space of the homogeneous system must have exactly N + 1 vectors, the changes of units that leave Γ invariant besides k. The references must then fix all of them (`MatrixRank` of their rows is N + 1).
```wl
{NullSpace[{{1, 1, 0}, {0, 1, 1}}], MatrixRank[{{1, 2}, {2, 4}}], LinearSolve[{{2, 0}, {0, 4}}, {2, 8}]}
```
:::

An amplitude written as a product, `nuk Ak lap[h]`, leaves one unit free too many:

```wl
ansatzA = Block[{Print}, DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk Ak lap[h]) - Dk hbar^2, "Couplings" -> {nuk, Ak, Dk}]];
showFailure[DeriveDimensions[ansatzA, "References" -> {nuk -> "nu", Dk -> "D"}, "Scale" -> k]]
```

The hint suggests the fix: `"Dimensionless" -> {Ak}`.

@region Kernel/Dimensions.wl | (* solve each basis component: | (overdetermined system consistent)"]];

The system is overdetermined (a coupling may appear in several monomials), so it is solved through the normal equations mᵀm u = −mᵀ c, and then the residual is checked to vanish exactly. A nonzero residual lists the monomials that disagree; `Nothing` is an element that removes itself from a list, so `If[ok, Nothing, label]` keeps only the conflicts.

```wl
{{1, Nothing, 2}, UnitVector[3, 2]}
```

`UnitVector[n, i]` is the dimension of the i-th basis element, used for the references.

## 4. An independent check: rescale the action

@quote Kernel/Dimensions.wl evalDensity rescalingCheckQ

The derived dimensions say how every quantity changes when the units change. The check makes such a change, with random factors λ_j, and verifies that the action is invariant: explicit polynomial test fields are put into the **density as you wrote it**, evaluated with a small interpreter `ev` of the operators (`dt` → ∂t, `grad` → list of ∂_i, ...), before and after the rescaling.

:::mma Unique and LCM | `Unique["x"]` creates a new symbol never used before (x1, x2, ...); `LCM` is the least common multiple.
The test coordinates are fresh symbols, so they cannot collide with anything in the density. The random factors are raised to the least common denominator of the exponents, so that fractional dimensions like k^{1−d/2} still give rational numbers.
```wl
{Unique["x"], Unique["x"], LCM[4, 6]}
```
:::

:::why Linear algebra on the parse, and a check on the raw density
The linear system uses the monomial table of lesson 02, so a parser error would give consistent but wrong dimensions. The rescaling check does not use the table at all: it reads the density through its own interpreter. This is the strategy of lesson 03 (plane waves against the parse) applied to dimensions. Checking at two values of d (2 and 3) and two random seeds makes an accidental agreement negligible.
:::

## 5. Names, η symbols, the stage

@region Kernel/Dimensions.wl | (* names of the dimensionless couplings and of the etas *) | PrintCard[stage];

A dimensionful coupling c gets the dimensionless name `a<>c` (λ → `alamk`), created next to it (`Symbol[Context[#] <> "a" <> SymbolName[#]]`, [[Symbol and SymbolName]]). A reference names its η: `nuk -> "nu"` gives `etanu`. Every generated name is checked against collisions.

## 6. Inspection

@quote Kernel/Dimensions.wl DimensionOf LinearPart CheckDimensions

`CheckDimensions` compares hand-derived dimensions with the derived ones. To turn a product of powers back into an exponent vector, it takes `PowerExpand[Log[e]]` and reads the coefficient of each `Log[basis element]`.

:::mma PowerExpand and Log | `PowerExpand[Log[x^a y^b]]` gives `a Log[x] + b Log[y]`, assuming positive quantities.
`Log` turns products of powers into sums, which `Coefficient` ([[Coefficient and Exponent]]) can read. `PowerExpand` makes the expansion unconditional, valid here because dimensions are positive scales.
```wl
PowerExpand[Log[k^2 nuk^(1/2)]]
```
:::

```wl
CheckDimensions[dims, <|"Frequency" -> k^2 nuk|>]
```

```wl
showFailure[CheckDimensions[dims, <|h -> k Dk|>]]
```

`EmbeddedStage` and `AnsatzOf` reach a stage embedded at any depth, for example the ansatz inside a loaded Equations stage (lesson 12):

@quote Kernel/Dimensions.wl EmbeddedStage AnsatzOf

`FirstCase[list, patt, default]` returns the first element matching patt.

## Exercises

:::exercise The KPZ coupling
Canet's coupling is g = λ² D/ν³ (lesson 11). Using `DimensionOf`, check that its dimension is k^{2−d}.
:::

:::solution
```wl
DimensionOf[dims, lamk]^2 DimensionOf[dims, Dk]/DimensionOf[dims, nuk]^3
```
The dimensions multiply like the quantities. So k^{d−2} g is dimensionless, the combination used in Canet's flow.
:::

:::exercise Other references
Choose ν and λ as references instead of ν and D. What is the dimension of D now, and what is its linear part?
:::

:::solution
```wl
dims2 = DeriveDimensions[ansatz, "References" -> {nuk -> "nu", lamk -> "l"}, "Scale" -> k];
{DimensionOf[dims2, Dk], LinearPart[dims2, Dk]}
```
D ∝ k^{2−d} ν³/λ², with linear part −2 + d − 2η_l + 3η_ν: the same physics in another basis.
:::

:::exercise Too few references
Give only one reference. Which check fails, and why does KPZ need two?
:::

:::solution
```wl
showFailure[DeriveDimensions[ansatz, "References" -> {nuk -> "nu"}, "Scale" -> k]]
```
`Dimensions.referenceCount`. With N = 1 scalar there are N + 1 = 2 free units: one fixes the dimensions of h (and so h̄), the other the dimension of time.
:::
