# Lesson 03 - Vertices

**Physics step.** Fourier-transform the parsed density (∂t → +iω, ∇ → −iq, all momenta incoming) and build every vertex Γ⁽ⁿ⁾, n = 2..4, as the symmetrised n-th functional derivative of Γ_k with the δ function stripped. Γ⁽²⁾ gives the matrix Γ2_ab(ω, q); Γ⁽³⁾ the cubic vertices of the one-loop flow.

**Stage.** Input: the Ansatz stage. Output: the `"Vertices"` stage, whose `"VertexTable"` maps each field multiset to its vertex.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
ansatz = lessonUpTo["Ansatz"];
```

## 1. Using it

```wl
vertices = BuildVertices[ansatz];
```

The table starts with the conventions (`docs/CONVENTIONS.md`), then Γ2 and the higher vertices. For KPZ: Γ2_{h h̄}(ω, q) = iω + ν q², Γ2_{h̄ h̄} = −2D, and one cubic vertex Γ⁽³⁾_{h h h̄}(ω1, q1; ω2, q2) = λ q1·q2. Three accessors evaluate vertices in the leg convention:

```wl
{Gamma2Matrix[vertices, w, q], VertexValue[vertices, {hbar, h, h}, {{w1, p1}, {w2, p2}}], Gamma3Matrix[vertices, hbar, {nu, p}, {w, q}]}
```

`VertexValue[v, {a1, ..., an}, {{ω1, q1}, ..., {ω(n−1), q(n−1)}}]` puts leg j at (ωj, qj) and the last leg at minus their sum. In `{hbar, h, h}` above, the last `h` carries −(p1 + p2), so λ p2·(−p1 − p2) = −λ(p1·p2 + p2²).

## 2. What is stored

```wl
vertices["VertexTable"]
```

The keys are field **multisets** sorted in field order: `{h, h, hbar}` stands for every ordering of the legs. The values are the vertices with *all n* momenta independent, written with private leg symbols `wv[l]`, `qv[l]` (frequency and momentum of leg l), before the momentum conservation is used.

:::why Store one vertex per multiset, evaluate any ordering on demand
Γ⁽³⁾_{h h h̄}, Γ⁽³⁾_{h h̄ h} and Γ⁽³⁾_{h̄ h h} are the same function with the legs relabelled. Storing the fully symmetrised vertex once, as a function of all n momenta, and relabelling at evaluation time keeps the table small (a quartic vertex of two fields would otherwise need up to 4! entries) and makes the symmetry impossible to get wrong: there is only one expression. The last leg is fixed by momentum conservation only when a vertex is evaluated, so the stored form never depends on which leg was chosen as "last".
:::

## 3. Fourier symbols of a monomial

@quote Kernel/Vertices.wl factorSymbol monomialAt

:::mma Complex numbers | `I` is the exact imaginary unit; `Im`, `Re`, `Conjugate` work on exact and symbolic expressions.
`(I wv[l])^nt (-I)^Length[ix]` is the symbol of ∂t^nt ∂_i... at leg l: one factor iω per time derivative and one −i per spatial index (the vector itself, `qv[l]`, is attached to each index by the rule `Thread[ix -> qv[l]]`). `I` is exact like a rational number, so `(-I)^2` is exactly `-1`.
:::

`factorSymbol` returns two things for each factor: its scalar prefactor, and rules saying which momentum each index carries. `monomialAt` then groups the rules by index (`GroupBy[..., First -> Last]`, [[GroupBy and KeyValueMap]]): each `SpaceIndex` appears on two factors, so it collects two momenta, and `scal @@ {qa, qb}` contracts them. The ε pair becomes a `wedge` ([[FixedPoint]] showed its algebra).

:::mma KeySelect and Values | `KeySelect[assoc, test]` keeps the entries whose key passes test; `Values[assoc]` lists the values.
`Values[KeySelect[g, FreeQ[#, EpsIndex] &]]` takes the momentum pairs of the ordinary indices and leaves out the ε slots, which are handled separately.
:::

For the cubic KPZ monomial −λ/2 h̄ ∇h·∇h, with h̄ at leg 1 and the two h at legs 2 and 3:

```wl
m = ansatz["Monomials"][[4]];
{FRGscalarMSR`Private`factorSymbol[m["Factors"][[2]], 2], FRGscalarMSR`Private`monomialAt[m, {1, 2, 3}]}
```

(−λ/2)(−i)² q2·q3 = λ/2 q2·q3.

## 4. Symmetrisation: summing over assignments of legs

@quote Kernel/Vertices.wl vertexExpr

The functional derivative of ∫ h̄ ∇h·∇h with respect to h̄(1) h(2) h(3) is the sum over the two ways of attaching the two h factors to legs 2 and 3. `vertexExpr` does exactly this: among all `Permutations` of the legs ([[Sort and canonical order]]) it keeps those that put each factor on a leg of the same field, and adds the symbols. Identical fields therefore give the multiplicity factors automatically (2 for h̄², 3! for φ̄φ³). Here the two terms are equal and give λ q2·q3, the entry of the table.

## 5. Building the table

@region Kernel/Vertices.wl | BuildVertices[ansatz_] := | table = Association[DeleteCases

:::mma Tuples | `Tuples[list, n]` gives every list of n elements drawn from list, with repetition and in every order.
`Tuples[{h, hbar}, 3]` has 8 elements. Sorting each in field order and taking the `Union` leaves the 4 multisets {h,h,h}, {h,h,h̄}, {h,h̄,h̄}, {h̄,h̄,h̄}. `Join @@ Table[..., {n, 2, nmax}]` concatenates the lists for n = 2, 3, 4 ([[Apply]] with head `Join`).
```wl
{Tuples[{h, hbar}, 2], Union[Sort /@ Tuples[{h, hbar}, 3]]}
```
:::

`sortF` is a local function ([[Local functions]]) sorting by position in the field order ([[DeleteCases, SortBy and FirstPosition]]); `Sort` alone would sort alphabetically, which is not the field order in general. `Association[DeleteCases[{key -> vertex, ...}, _ -> 0]]` drops the multisets with no vertex: the pattern `_ -> 0` matches any rule whose right-hand side is 0.

## 6. Evaluating in the leg convention

@quote Kernel/Vertices.wl fieldPos VertexValue Gamma2Matrix Gamma3Matrix

:::mma Ordering | `Ordering[list]` gives the positions that would sort list: `list[[Ordering[list]]]` is sorted.
`VertexValue` receives the fields in your order, say `{hbar, h, h}`. `ord = Ordering[...]` finds which of your legs goes to each slot of the canonical multiset `{h, h, hbar}`; the stored `wv[j]`, `qv[j]` are then replaced by the frequency and momentum of your leg `ord[[j]]`.
```wl
{Ordering[{3, 1, 2}], {c, a, b}[[Ordering[{c, a, b}]]]}
```
:::

:::mma Key | `Key[k]` wraps a key so that it is never read as something else; `Lookup[assoc, Key[{a, b}], 0]` looks up the key `{a, b}`.
The keys of the vertex table are lists. `Lookup[table, {h, hbar}]` would mean "look up the keys `h` and `hbar`" (a list of keys) and return two missing values. `Key` says "this list is one key".
```wl
{Lookup[vertices["VertexTable"], {h, hbar}, 0], Lookup[vertices["VertexTable"], Key[{h, hbar}], 0]}
```
:::

`Table[expr, {a, list}, {b, list}]` with two iterators ([[MapAt, Table and Sum]]) builds a matrix whose entries run over the field order: that is all `Gamma2Matrix` and `Gamma3Matrix` do.

## 7. An independent check: plane waves into the raw density

The card's first line says that the vertices agree with a second derivation that shares no code with the first:

@quote Kernel/Vertices.wl planeWaveCoefficient symbolAsPlaneWave

Each field is replaced by a sum of plane waves, φ_a(t, x) = Σ_l ε_l e^{i(ω_l t − q_l·x)}, with explicit 3-component vectors (2 in a declared d = 2). The package operators are replaced by Mathematica's own derivatives, acting on explicit coordinates `xx[1], xx[2], xx[3]`. The coefficient of ε_1⋯ε_n is then the vertex, and it must equal the table entry written in components.

:::mma Array and explicit vectors | `Array[f, n]` is `{f[1], ..., f[n]}`; `.` (Dot) contracts explicit lists; `D`, `Grad`, `Div`, `Laplacian` differentiate with respect to explicit variables.
With `X = Array[xx, 3]` and `kq[l] = {kc[l, 1], kc[l, 2], kc[l, 3]}`, `kq[l] . X` is the explicit sum kc[l,1] xx[1] + ... . `Grad[e, X]` is the gradient with respect to those three variables. These are the System operators the ansatz forbids: here they are safe, because the coordinates are explicit symbols that the fields really depend on.
```wl
Grad[Exp[I (w t - {k1, k2} . {x1, x2})], {x1, x2}]
```
:::

:::mma Replacing heads by functions | A rule `dt -> (D[#, tt] &)` turns `dt[e]` into `(D[#, tt] &)[e]`, which evaluates to `D[e, tt]`.
Replacing the *head* of `dt[...]` by a pure function ([[Pure functions]]) gives the inert operator its meaning in one rule. `D[e, Sequence @@ Table[eps[l], {l, n}]]` then differentiates once with respect to each ε_l ([[Sequence]]); setting all ε to 0 afterwards keeps exactly the coefficient of ε_1⋯ε_n.
:::

```wl
FRGscalarMSR`Private`planeWaveCoefficient[ansatz["Density"], {h, hbar}, {h, h, hbar}, 3]
```

This is λ (k1·k2) written in components, as the table says.

:::why Two derivations instead of one
The first derivation (sections 3–4) is fast and produces the structured form the pipeline needs, but it depends on the parser of lesson 02 and on the Fourier rules. The second goes back to the density you typed, uses only Mathematica's own calculus, and knows nothing about `FieldFactor`, canonical forms or symbol rules. If both agree on every multiset (including those that should have no vertex: nothing may be lost), the parser and the transform are right. This is the library's general strategy: when a step can be computed in two independent ways, compute both and compare exactly. The second way is too slow for the later steps, but cheap enough here.
:::

## 8. The other checks

@region Kernel/Vertices.wl | stage = <|"Fields" -> fields | checks = Join[checks

`Gamma2Matrix` is first called on a provisional association that has only the keys it needs (`"Fields"`, `"VertexTable"`): the accessors are plain functions of associations ([[Association]]), not tied to a finished stage.

The three "(implementation check)" lines are near-tautological: by construction they cannot fail unless the code above is broken (review E11). They are labelled as such so that the card does not pretend to test more than it does. `PolynomialQ[e, {w, ss}]` tests whether e is a polynomial in those variables; `ConstantArray[0, dims]` is a zero matrix of the given dimensions.

The last check is physical. The free response function G_{φ̄φ}(ω) = 1/Γ2_{φφ̄}(ω) must have its pole in the upper half ω-plane (causality, with the Fourier sign of `docs/CONVENTIONS.md`).

:::mma Solve | `Solve[eqn, x]` gives the solutions as lists of rules `{{x -> s1}, ...}`; `x /. Solve[...]` extracts the values.
The couplings get random positive rationals (`RandomRationalPoint`, lesson 00), q² is set to 1/2, and `Solve` finds the poles exactly. `Exponent[e, w]` is the degree of e in w; `Im` takes the imaginary part. The result is checked with [[Which]]: no pole (WARN), a pole below the real axis (FAIL), a pole on the axis (WARN, no damping), otherwise PASS.
:::

A damping term with the wrong sign is caught:

```wl
badAnsatz = Block[{Print}, DefineAnsatz[lessonUpTo["Theory"], "Density" -> hbar (dt[h] + nuk lap[h]) - Dk hbar^2, "Couplings" -> {nuk, Dk}]];
showFailure[BuildVertices[badAnsatz]]
```

:::why A physical convention check
A sign error in the Fourier convention gives the time-mirrored theory: every formula is consistent, every algebraic check passes, and the flow equations are wrong. Only a physical requirement can catch it: the response function must be retarded. The check costs one `Solve` per scalar and is the first place in the pipeline where a convention is tested against physics rather than against itself (review E1).
:::

## 9. Display

@quote Kernel/Vertices.wl legSymbols ShowVertex ShowGamma2

The display uses its own symbols, created in yet another context, ``FRGscalarMSR`Vertices`Private``, with names `ω1`, `q1`, .... `displayString` (lesson 00's `Core`) strips every ``FRGscalarMSR`…`` context when printing, so they print as plain `ω1`, but they can never collide with a symbol of yours called `q1`. `Collect[e, w]` groups the terms of e by powers of w. `Do[expr, {a, n}, {b, n}]` is a loop with two iterators, for side effects only.

```wl
ShowVertex[vertices, {hbar, h, h}]
```

## Exercises

:::exercise Γ2 by hand
From the KPZ density, derive Γ2_{h̄ h}(ω, q) by hand (leg h̄ at (ω, q), leg h at (−ω, −q)) and compare with `Gamma2Matrix`.
:::

:::solution
h̄ ∂t h gives iω at the h leg, here at −ω: −iω. −ν h̄ ∇²h gives −ν (−i)² (−q)² = +ν q². So Γ2_{h̄ h}(ω, q) = −iω + ν q²:
```wl
Gamma2Matrix[vertices, w, q][[2, 1]]
```
:::

:::exercise The symmetry of identical legs
Check with `VertexValue` that Γ⁽³⁾_{h h h̄} does not change when the two h legs exchange their momenta, and that Γ⁽³⁾_{h h̄ h} at legs (p1, p2) is a relabelling of Γ⁽³⁾_{h h h̄}.
:::

:::solution
```wl
{Expand[VertexValue[vertices, {h, h, hbar}, {{w1, p1}, {w2, p2}}] - VertexValue[vertices, {h, h, hbar}, {{w2, p2}, {w1, p1}}]],
 Expand[VertexValue[vertices, {h, hbar, h}, {{w1, p1}, {w2, p2}}] - VertexValue[vertices, {h, h, hbar}, {{w1, p1}, {-w1 - w2, -p1 - p2}}]]}
```
In the second, leg 3 of `{h, hbar, h}` (an h at −(p1+p2)) moves to slot 2, and h̄ moves to the last slot, where it automatically sits at −(p1 + (−p1 − p2)) = p2.
:::

:::exercise A new vertex, and a check that cannot run
Add the term `g hbar grad[h].grad[hbar]` to the KPZ density (with a new coupling `g`) and look at the new entry of the vertex table. Then remove `hbar dt[h]` too: what does the card say?
:::

:::solution
```wl
a2 = Block[{Print}, DefineAnsatz[lessonUpTo["Theory"], "Density" -> hbar (dt[h] - nuk lap[h]) - Dk hbar^2 + g hbar grad[h].grad[hbar], "Couplings" -> {nuk, Dk, g}]];
Block[{Print}, BuildVertices[a2]]["VertexTable"]
```
```wl
a3 = Block[{Print}, DefineAnsatz[lessonUpTo["Theory"], "Density" -> - nuk hbar lap[h] - Dk hbar^2, "Couplings" -> {nuk, Dk}]];
BuildVertices[a3];
```
The new cubic vertex is `{h, hbar, hbar}`. Without the `dt` term, Γ2_{h h̄} has no ω, the pole check has nothing to test, and the card shows the WARN `Vertices.retardedPole[h]`.
:::
