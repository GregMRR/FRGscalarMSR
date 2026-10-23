# Lesson 08 - Projection

**Physics step.** The flow of a coupling is the coefficient of its monomial (iν)^m (p²)^n in ∂sΓ2_ab(ν, p). The module extracts it from the loop integrand: (∂/∂(iν))^m/m! at ν = 0, and the (p²)^n coefficient at p = 0 as the direction average of the s^{2n} Taylor coefficient along p = s e (`design/MOMENTUMPROJECTION.md`, M1, valid for any n). The external variables are set to zero first ("pipeline 1"), so the only frequency poles left are those of Δ(±ω).

**Stage.** Input: the Loops and Symmetries stages. Output: the `"Projection"` stage, with one ω-integrand per classified monomial (every class except Structure), as a function of ω and q² = `LoopMomentumSquared`.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
loops = lessonUpTo["Loops"]; symmetries = kpzSymmetries;
```

## 1. Using it

```wl
projected = ProjectOnMonomials[loops, symmetries];
```

16 integrands: the 18 monomials of lesson 06 minus the 2 of class Structure. The one that gives the flow of ν:

```wl
ShowIntegrand[projected, nuk]
```

It is a rational function of the atoms at q² (written `LoopMomentumSquared`), of Δ(±ω), and of the loop frequency `wL`. `ShowIntegrand[projected, nuk]` finds the record by its Γ_k coefficient; `ShowIntegrand[projected, {h, hbar}, {0, 1}]` would find it by pair and orders.

## 2. The operator

@quote Kernel/Projection.wl project

`project` takes the integrand at general (ν, p), differentiates m times in ν, takes the p-projection, sets ν = 0, renames q² and completes the derivatives of Δ (section 4). `eps_ : False` is an optional argument with a default ([[Patterns]]).

:::mma Repeated derivatives | `D[e, {x, n}]` is the n-th derivative of e in x; `D[e, {x, n}, {y, m}]` mixes variables.
`D[e, {nu, m}]` followed by `/I^m m!` and ν → 0 is the coefficient of (iν)^m in the Taylor series. The same `D[e, {projS, 2 n}]/(2 n)!` at s = 0 gives the s^{2n} coefficient in `directionalAverage`.
:::

## 3. The momentum projection: a direction average

@quote Kernel/Projection.wl projCAverage directionalAverage

The loop depends on p only through `scal[p, p]` and `scal[p, q]` (`scal` is bilinear, so (p + q)² is already expanded, lesson 00, section 7). Write p = s e with e a unit vector at angle c = cos θ to q: `scal[p, p]` → s², `scal[p, q]` → s c. The s^{2n} coefficient is a polynomial in c, and averaging over the directions e uses ⟨c^{2k}⟩ = (2k−1)!!/(d(d+2)⋯(d+2k−2)) times |q|^{2k}; odd powers average to zero.

:::mma Special functions and products | `k!!` is the double factorial, `Product[f, {j, a, b}]` a product, `Pochhammer[a, n]` = a(a+1)⋯(a+n−1), `Gamma[x]` the Γ function; `FunctionExpand` rewrites special values in elementary form.
These turn the direction average into a closed rational function of d, symbolic to the end ([[Iterators that depend on each other]] for `Product`'s iterator).
```wl
{7!!, Pochhammer[d/2, 2], FunctionExpand[Gamma[5/2]]}
```
:::

```wl
{FRGscalarMSR`Private`directionalAverage[scal[p, q]^2, 1, p, q, Dim], FRGscalarMSR`Private`directionalAverage[scal[p, q]^4, 2, p, q, Dim]}
```

⟨(p·q)²⟩ = p² q²/d gives q²/d as the (p²)¹ coefficient; at order p⁴, 3 q⁴/(d(d+2)).

:::why A direction average instead of a Laplacian
The usual trick takes (∇²_p)ⁿ at p = 0 and divides by its value on (p²)ⁿ. It needs explicit vectors or a hand-coded Laplacian on `scal` expressions, and the previous version of the library got it right only for n ≤ 1. The direction average works with `scal[p, p]` and `scal[p, q]` as two scalars, any n, with one table of ⟨c^{2k}⟩. The old operator is kept as a test-only reference (`tests/reference/OldProjection.wl`), and the new one is required to agree with it on every existing integrand (`design/MOMENTUMPROJECTION.md`).
:::

## 4. Derivatives of Δ

Differentiating `Δ[nu + wL, s]` with respect to ν produces the derivative of the inert function Δ:

```wl
D[Δ[nu + w, s], nu]
```

:::mma Derivative of several arguments | `Derivative[i, j][f][x, y]` is ∂ⁱ_x ∂ʲ_y f; a pattern `Derivative[i_, j_][Δ][w_, s_]` matches any such derivative.
`Δ` has no definition, so `D` writes its derivatives with `Derivative` ([[D and unknown functions]]). `deltaCompletion` is a rule that replaces any derivative of Δ by the derivative of its polynomial form Σ_m c_m(s) (iω)^m, so derivatives come out in the atoms `DeltaCoefficient[m][s]` and their s-derivatives.
:::

@quote Kernel/Projection.wl deltaCompletion

```wl
D[Δ[nu + w, s], nu] /. FRGscalarMSR`Private`deltaCompletion[1]
```

## 5. Building the integrands, and the checks

@region Kernel/Projection.wl | ProjectOnMonomials[loops_, symmetries_] | Append[r, "Integrand" -> ig]]] /@ ms]], byPair];

The integrands are computed per pair: `GroupBy[targets, #Pair &]` ([[GroupBy and KeyValueMap]]) groups the monomials, and `LoopIntegrand` (lesson 07) is called once per pair, then projected once per monomial.

Three checks guard the operator, each independent of the others:
- **`Projection.LHS`**: the same `project`, applied to Γ_k itself (regulator off), must give back each classified coefficient. It tests the frequency normalisation and the order in s; Γ_k depends on p only through p², so it cannot see the d-normalisation (review E2).
- **`Projection.normalisation`**: the projection of (p·q)^{2k} (p²)^{n−k} against the closed form q^{2k} Γ(k+½)Γ(d/2)/(√π Γ(k+d/2)), written with `Gamma` instead of `projCAverage` (review E17).
- **`Projection.explicitVectors`**, for n ≥ 2: the result against (∇²_p)ⁿ/N_n computed on explicit d-component vectors, at the declared d, or at d = 1, …, n+1 for a symbolic d.

@quote Kernel/Projection.wl checkNormalisation lapNorm explicitProjection

:::mma Nest | `Nest[f, x, n]` applies f to x n times: `f[f[f[x]]]` for n = 3.
`Nest[Function[g, Sum[D[g, {x, 2}], {x, pc}]], e, n]` applies the explicit Laplacian n times. (`Nest` appeared in lesson 02's `fieldString` on strings; here it iterates an operator.)
```wl
{Nest[f, x, 3], FRGscalarMSR`Private`explicitProjection[scal[p, q]^2, 0, 1, {nu, p}, {w, q}, 3, 1]}
```
:::

In 3 explicit dimensions, (∇²_p)(p·q)² / N_1 = 2q²/6 = q²/3, the same as q²/d at d = 3.

:::why Why the explicit check is needed only for n ≥ 2, at n+1 dimensions
A p^{2n} projection is a rational function of d whose numerator can be written with n+1 independent d-structures, so agreeing at n+1 integer dimensions pins it down for every d (review E4). For n ≤ 1 the regression against the old operator certifies the result already, and the KPZ example of these lessons has n ≤ 1, so its card has no `explicitVectors` line. The p⁴ flows are tested against the literature in `benchmarks/KuramotoSivashinsky`.
:::

## 6. With ε (d = 2)

@quote Kernel/Projection.wl dropWedge

With `jac`, products of wedges are first reduced to scalar products (`reduceWedge2D`, lesson 00, section 7). Every term is then wedge[p, q]^e f(p², p·q, q², ν) with e = 0 or 1, and the e = 1 part has a vanishing projection for every m and n: reflecting e about q keeps c and flips the sign of wedge[e, q] (M3, `design/LEVICIVITA.md` F4). It is dropped, after checking on explicit 2-vectors that its projection really vanishes.

`$wedgeDropped` is one of the few global variables of the library: a private counter, reset by `ProjectOnMonomials` and read for its card line. `$wedgeDropped++` adds 1 to it.

## 7. The stage

@region Kernel/Projection.wl | stage = stampStage[<|"Stage" -> "Projection" | PrintCard[stage];

@quote Kernel/Projection.wl ShowIntegrand

:::mma SelectFirst and Missing | `SelectFirst[list, test]` gives the first element that passes test, or `Missing["NotFound"]`; `MissingQ` tests for it.
`ShowIntegrand` has two definitions ([[Several definitions of one function]]): one by pair and orders, one by Γ_k coefficient, which searches the records with `SelectFirst` and then calls the first.
```wl
{SelectFirst[{1, 2, 3}, # > 1 &], SelectFirst[{1, 2, 3}, # > 5 &], MissingQ[SelectFirst[{1, 2, 3}, # > 5 &]]}
```
:::

## Exercises

:::exercise The frequency projection by hand
Apply `project` (with m = 1, n = 0) to `I nu + 3 nu^2 + nuk scal[p, p]`. What do you expect, and why does the `scal[p, p]` term not contribute?
:::

:::solution
```wl
FRGscalarMSR`Private`project[I nu + 3 nu^2 + nuk scal[p, p], 1, 0, {nu, p}, {w, q}, Dim, 1]
```
1: the coefficient of (iν)¹. The ν² term has no linear part, and the p² term is removed by the n = 0 projection (p = 0).
:::

:::exercise The p⁴ direction average against the closed form
Check `directionalAverage[scal[p, q]^4, 2, p, q, Dim]` against q⁴ Γ(5/2)Γ(d/2)/(√π Γ(2 + d/2)).
:::

:::solution
```wl
Together[FunctionExpand[FRGscalarMSR`Private`directionalAverage[scal[p, q]^4, 2, p, q, Dim] - scal[q, q]^2 Gamma[5/2] Gamma[Dim/2]/(Sqrt[Pi] Gamma[2 + Dim/2])]]
```
0: this is the k = 2 case of `Projection.normalisation`.
:::

:::exercise Integrands of the protected monomials
Which records of `projected["Integrands"]` have a class whose flow must vanish? Are their integrands zero already? When does the vanishing get checked?
:::

:::solution
```wl
{#Pair, #Frequency, #Momentum, #Class, LeafCount[#Integrand]} & /@ Select[projected["Integrands"], MemberQ[{"Protected", "Causality", "DeclaredNumeric"}, #Class] &]
```
Some are already 0 (LeafCount 1), most are not: the vanishing holds only after the frequency integral (and, for Ward identities, sometimes only after the q integral). `CheckClosure` (lesson 10) tests it after the frequency integration.
:::
