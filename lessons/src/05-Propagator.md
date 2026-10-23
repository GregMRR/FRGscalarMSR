# Lesson 05 - Propagator

**Physics step.** Invert Γ2 + R to get the propagator G(ω, q), using the MSR block structure: the (φ, φ) block vanishes by causality, so G has the retarded response G_{φ̄φ} = B(ω)⁻¹, the advanced one, and the correlation G_{φφ} = −(B(−ω)ᵀ)⁻¹ C B(ω)⁻¹, with B the (φ, φ̄) block and C the noise block. The poles are the roots of Δ(ω) = det B(ω); the module states the conditions for them to be causal.

**Stage.** Input: the Vertices and Regulator stages, built from the same Ansatz. Output: the `"Propagator"` stage.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
regulator = lessonUpTo["Regulator"]; vertices = kpzVertices;
```

`lessonUpTo` returns the last stage it built and leaves the earlier ones in variables (`kpzTheory`, `kpzAnsatz`, `kpzVertices`, ...), as its printout says.

## 1. Using it

```wl
propagator = BuildPropagator[vertices, regulator];
```

Γ2 + R and G, in **atoms**:

```wl
{FullMatrix[propagator, w, s], PropagatorMatrix[propagator, w, s]}
```

`GammaAtom[h, hbar, 0][s]` stands for the coefficient of (iω)⁰ in (Γ2 + R)_{h h̄} at q² = s, that is ν s (1 + r_n(s/k²)). `Δ[w, s]` stands for Δ(ω) at q² = s. Their values are kept aside:

```wl
AtomRules[propagator]
```

Substituting them gives the familiar KPZ propagator (here, and only here, on a small expression, `Simplify` is harmless):

```wl
Simplify[PropagatorMatrix[propagator, w, s] /. DeltaRules[propagator] /. AtomRules[propagator]]
```

:::why Abstract atoms until the end
The explicit entries contain couplings and shape functions; G built from them would contain them many times over, and every later expression (loop integrands, their projections, the frequency integrals) would grow accordingly. With one atom per entry coefficient, G stays a small rational function of a handful of symbols, the heavy algebra of lessons 07–10 runs on it, and the atoms are replaced by their values only once, at assembly (lesson 11). This is rule 1 of the computer-algebra rules in the shallow-water project primer, "abstract until the end" (`~/sync_cnrs/Projects/shallow_water/AGENT_PRIMER.md` §7, cited in `design/FRGDESIGN.md`, module 5). The atoms are curried, `GammaAtom[a, b, m][s]` ([[Curried definitions]]), because the same coefficient is needed at different momenta, q² and (q + p)²; `ExactlyZeroQ` treats each `GammaAtom[...][s]` as one atom with its own random value (lesson 00, section 5).
:::

## 2. From Γ2 + R to atoms

@region Kernel/Propagator.wl | BuildPropagator[vertices_, regulator_] | ToString[Length[atoms]] <> " real entry atoms"]];

`Propagator.sameAnsatz` compares the `"Upstream"` fingerprints of the two inputs ([[Hash]]): both must point to the same Ansatz stage.

`Gamma2Matrix` (lesson 03) and `RegulatorMatrix` (lesson 04) are evaluated with private symbols `ww` for ω and `ss` for q².

:::mma Coefficient and Exponent | `Coefficient[e, x, m]` is the coefficient of x^m in the polynomial e; `Exponent[e, x]` its degree in x.
`Coefficient[e, ww, m] / I^m` is the coefficient of (iω)^m. For each entry and each power, the module either keeps a number (`NumericQ[c]` is True for exact or approximate numbers, like the 1 in front of iω) or creates the atom `GammaAtom[a, b, m]` and stores its value as a function of q² ([[Function and Evaluate]]) in the local association `atomVals`.
```wl
{Coefficient[3 + 2 x + 5 x^2, x, 2], Exponent[3 + 2 x + 5 x^2, x], NumericQ[3/4], NumericQ[nuk]}
```
:::

Only the entries (a, b) with a before b in field order get atoms; the entry (b, a) is the same polynomial at −ω, which is the factor `I If[i <= j, 1, -1] ww`. `Propagator.atomForm` substitutes the values back and requires the result to equal Γ2 + R exactly.

:::mma Assigning to an association | `assoc[key] = value` adds or replaces one entry of an association stored in a variable.
Unlike `Append` ([[Association]]), this modifies the variable in place, like `AppendTo` for lists ([[AppendTo]]). It is used here inside a `Module`, on a local variable: the stage itself is still built once and never modified.
:::

## 3. The block structure and Δ

@region Kernel/Propagator.wl | (* block structure in the order | "internal"];

:::mma Part with lists | `m[[{1, 3}, {2}]]` extracts the rows 1, 3 and the column 2 of a matrix; the same form on the left of `=` assigns a block.
`abs[[iphi, ibar]]` is the (φ, φ̄) block B: the rows at the positions of the physical fields, the columns at those of the response fields. `Det` is the determinant. `B[w_] := ... /. ww -> w` is a local function giving B at any frequency.
```wl
mm = ConstantArray[0, {3, 3}]; mm[[{1, 3}, {1, 3}]] = {{1, 2}, {3, 4}}; {mm, mm[[{1, 3}, {3}]]}
```
:::

`Propagator.deltaDegree` requires 1 ≤ deg_ω Δ ≤ 2: degree 0 means there is no ∂t term and the frequency integrals diverge; degree 3 or more is beyond v1 (`docs/FUTURE.md`). With two scalars, deg Δ = 2 is the common case; KPZ has degree 1.

## 4. The inverse in closed form

@region Kernel/Propagator.wl | (* G = M^-1 with M | "G_{phibar phibar} = 0"]];

@quote Kernel/Propagator.wl adjugate

The inverse of B is its adjugate divided by its determinant. The adjugate is computed from cofactors: `Drop[m, {i}, {j}]` removes row i and column j, `Det` of what remains is the minor. Instead of the determinant itself, the formula uses the inert symbol `Δ[w, s]` (and `Δ[-w, s]` for the advanced block), so G is written in atoms and Δ only.

:::mma Matrix algebra | `.` multiplies matrices (and vectors); `Transpose`, `Det`, `IdentityMatrix[n]` work on nested lists.
Matrices are lists of rows, nothing more: `{{a, b}, {c, d}} . {{e, f}, {g, h}}` is the matrix product. `-dinv . C . binv` is the correlation block −(B(−ω)ᵀ)⁻¹ C B(ω)⁻¹, and `Map[Together, gpp, {2}]` ([[Map at a level]], [[Exact numbers]]) puts each entry over a common denominator.
```wl
{{{a, b}, {c, d}} . {{1, 0}, {0, 2}}, Det[{{a, b}, {c, d}}], Drop[{{a, b}, {c, d}}, {1}, {2}]}
```
:::

`Propagator.inverse` and `Propagator.inverseRight` multiply G by Γ2 + R in both orders, with Δ replaced by its value, and require the identity matrix exactly ([[Exact numbers]] through `ExactlyZeroQ`): every atom, ω and q² gets a random rational.

:::why Use the structure, then check the result
A generic symbolic `Inverse` of a 4×4 matrix of polynomials works, but produces a large expression with the determinant expanded everywhere and no visible structure. Using the MSR block form gives a short closed formula in which each block has a physical meaning and Δ appears explicitly, as needed by the frequency integration (lesson 09). Because the formula is written by hand, it is then checked against the definition of an inverse, in both orders, which costs a few rational evaluations.
:::

## 5. Validity of the poles: an assumption, not a check

@region Kernel/Propagator.wl | (* validity: roots of Delta | " at every loop momentum q (not checked

With u = iω, a root of Δ in the upper half ω-plane is a root u with negative real part. For degree 1, Δ = c0 + c1 u has its root at u = −c0/c1, so the condition is c0/c1 > 0; for degree 2 it is c1/c2 > 0 and c0/c2 > 0. The conditions are expressions like `DeltaCoefficient[0][s]/DeltaCoefficient[1][s] > 0`, which stay unevaluated as long as their atoms have no values: they are data ([[If and TrueQ]] explained why such an inequality cannot be decided).

```wl
ValidityConditions[propagator][s]
```

:::why Why an ASSUMPTION line instead of a check
Whether the poles are causal depends on the values of the couplings at every q, which are only known while the flow is integrated numerically. The library cannot decide it, so it does not pretend to: the card shows an ASSUMPTION, Adimensional (lesson 14) rewrites the conditions in dimensionless form, and ExportJulia (lesson 15) writes them as `causal_poles` for FlowRG to check during integration (`design/CAUSALITY.md`). An honest card distinguishes what was verified (PASS) from what was assumed.
:::

## 6. The stage and its accessors

@region Kernel/Propagator.wl | stage = stampStage[<|"Stage" -> "Propagator" | "Checks" -> checks, "Warnings" -> warns

@quote Kernel/Propagator.wl PropagatorMatrix FullMatrix AtomRules DeltaRules ValidityConditions

`DeltaRules` returns rules with patterns on the left, `Δ[w_, s_] :> ...`: one rule replaces Δ at every frequency and momentum where it occurs ([[Rules and replacement]]). `Δ` is a symbol whose name is a Greek letter (written `\[CapitalDelta]` in the source file); Mathematica allows such names like any other.

## Exercises

:::exercise The correlation function
From the result of section 1, write G_{hh}(ω, q) and check that it is 2D(1 + r_d) / (ω² + ν² q⁴ (1 + r_n)²).
:::

:::solution
```wl
gExplicit = PropagatorMatrix[propagator, w, s] /. DeltaRules[propagator] /. AtomRules[propagator];
ExactlyZeroQ[gExplicit[[1, 1]] - 2 Dk (1 + rd[s/k^2])/(w^2 + nuk^2 s^2 (1 + rn[s/k^2])^2)]
```
:::

:::exercise Why G_{h̄h̄} = 0
Using the block form, explain why G_{φ̄φ̄} always vanishes, and check it on `propagator`.
:::

:::solution
The inverse of [[0, B(ω)], [B(−ω)ᵀ, C]] has a zero block where the original has C: multiplying out shows that the (φ̄, φ̄) block of the inverse must be 0 for the product to be the identity. This is the statement that ⟨φ̄ φ̄⟩ = 0 in MSR.
```wl
PropagatorMatrix[propagator, w, s][[2, 2]]
```
:::

:::exercise Mixing two ansätze
Build a second KPZ ansatz without the λ term, a regulator on it, and pass that regulator with the original `vertices` to `BuildPropagator`. Which check fails?
:::

:::solution
```wl
ansatz2 = Block[{Print}, DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h]) - Dk hbar^2, "Couplings" -> {nuk, Dk}]];
regulator2 = Block[{Print}, DefineRegulator[ansatz2, {hbar, h} -> nuk q^2 rn[q^2/k^2], "Shapes" -> {rn}]];
showFailure[BuildPropagator[vertices, regulator2]]
```
`Propagator.sameAnsatz`: the two stages point to different Ansatz fingerprints, even though the propagators would happen to agree. A stage only accepts inputs with a common history.
:::
