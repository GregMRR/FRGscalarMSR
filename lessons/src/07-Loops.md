# Lesson 07 - Loops

**Physics step.** The one-loop flow of Γ2, from two field derivatives of the Wetterich equation ∂sΓ = ½ Tr[∂sR (Γ⁽²⁾ + R)⁻¹]:

∂sΓ2_ab(ν, p) = ∫dω/(2π) ∫dᵈq/(2π)ᵈ { Tr[∂sR G(ω,q) Γ3_a(ν,p; ω,q) G(ω+ν, p+q) Γ3_b(−ν,−p; ω+ν, p+q) G(ω,q)] − ½ Tr[∂sR G Γ4_ab G] }.

The module builds the integrand at a fixed loop point (ω, q), with the 1/(2π) of the frequency measure included; the integral over q is left to the numerics. The derivation (steps 1–6) is written out in `docs/MODULES.md`, section Loops.

**Stage.** Input: the Vertices, Propagator and Regulator stages. Output: the `"Loops"` stage, which holds its ingredients; the integrand of each entry is computed on request.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
propagator = lessonUpTo["Propagator"];
```

## 1. Using it

```wl
loops = BuildLoopIntegrands[kpzVertices, propagator, kpzRegulator];
```

The card has no PASS line that tests the formula, and an ASSUMPTION line that says what was tested elsewhere: the bubble against Canet's KPZ flow (`tests/acceptance/KPZ.wl`), the tadpole not at all.

```wl
{wq = LoopVariables[loops], LeafCount[LoopIntegrand[loops, {hbar, hbar}, nu, p]]}
```

The loop variables are private symbols (`wL`, `qL`): they appear in the integrand, and only the next modules, which integrate over them, need to know them.

## 2. The whole module

@quote Kernel/Loops.wl LoopVariables gMat gamma4Matrix LoopIntegrand BuildLoopIntegrands

`gMat` evaluates the propagator's stored function ([[Function and Evaluate]]) at a frequency and at the square of a momentum vector, `scal[vec, vec]`: with `vec = p + qL`, the regulator and the atoms inside G are automatically at (p + q)² ([[Evaluation as rewriting]] expands `scal[p + qL, p + qL]`). `Gamma3Matrix` (lesson 03) gives the vertex matrices; their product is a [[Matrix algebra]] product.

:::mma Tr and exact constants | `Tr[m]` is the trace of a matrix; `Pi`, `E`, `I` are exact constants, never rounded.
`(bubble - tadpole/2)/(2 Pi)` keeps π as a symbol: the integrand is exact, and the final equations carry exactly the factors of 2π the derivation puts there.
```wl
{Tr[{{a, b}, {c, d}}], InputForm[Pi/2], N[Pi/2]}
```
:::

The tadpole is included only when there is a quartic vertex (`"HasGamma4"`), and then with a WARN: its prefactor −½ is derived by hand but has never been tested against an independent result.

:::why The integrand as a function of the entry
The stage stores the vertices, the propagator, the regulator and the abstract ∂sR matrix, and `LoopIntegrand[loops, {a, b}, nu, p]` assembles the integrand only when it is asked for. The projection (lesson 08) needs it for a few entries only, and expanded around different external momenta; building all entries in advance would compute and store expressions nobody uses. Holding the full upstream stages inside the Loops stage also keeps it self-contained: a later module can reach the vertices or the propagator through it, and its fingerprint covers them.
:::

:::why Say on the card what has not been tested
The loop formula is the one part of the pipeline that is only derived by hand, never computed twice (review E3/E4). The library cannot check it by itself; an acceptance test fixes the bubble's prefactor, sign and 2π against an independent result from the literature (Canet, arXiv:2509.01472, Eq. (20)), and nothing yet tests the tadpole. Instead of a PASS line that would suggest otherwise, the card carries an ASSUMPTION and, when Γ⁽⁴⁾ ≠ 0, a WARN. The benchmarks in `benchmarks/` exist for the same reason: they test the pipeline against papers.
:::

## 3. What the integrand looks like

```wl
LoopIntegrand[loops, {hbar, hbar}, nu, p]
```

Everything is still abstract: the atoms `GammaAtom[...][...]` and `Δ[...]` (lesson 05) at q² and at (p + q)², the regulator derivatives `dsR[...]` (lesson 04), and the couplings of the vertices. Each Δ factor carries its frequency and momentum, which is what the frequency integration (lesson 09) needs to find the poles.

## Exercises

:::exercise The poles of the integrand
List the different `Δ[...]` factors of the h̄h̄ integrand. Which are evaluated at the loop momentum q, which at p + q?
:::

:::solution
```wl
Union[Cases[LoopIntegrand[loops, {hbar, hbar}, nu, p], _Δ, Infinity]]
```
Δ(±ω) at q² come from the two G(ω, q); Δ(±(ω + ν)) at (p + q)² from G(ω + ν, p + q). `_Δ` matches any expression with head `Δ` ([[Patterns]]).
:::

:::exercise A quartic vertex
Add the quartic term `g hbar h^3` to the KPZ density (coupling `g`) and build the chain up to the loops. What does the card say?
:::

:::solution
```wl
Block[{Print},
  a4 = DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h] - lamk/2 grad[h].grad[h]) - Dk hbar^2 + g hbar h^3, "Couplings" -> {nuk, lamk, Dk, g}];
  v4 = BuildVertices[a4];
  r4 = DefineRegulator[a4, {hbar, h} -> nuk q^2 rn[q^2/k^2], {hbar, hbar} -> -2 Dk rd[q^2/k^2], "Shapes" -> {rn, rd}];
  p4 = BuildPropagator[v4, r4]];
BuildLoopIntegrands[v4, p4, r4];
```
A WARN `Loops.tadpole`: Γ⁽⁴⁾_{h̄hhh} ≠ 0, so the tadpole term is included, and it is untested. (`Block[{Print}, ...]` hides the cards of the intermediate stages, as `lessonUpTo` does.)
:::

:::exercise A propagator from elsewhere
Pass `kpzVertices`, the propagator `p4` of the previous exercise, and `kpzRegulator` to `BuildLoopIntegrands`. Which check fails?
:::

:::solution
```wl
showFailure[BuildLoopIntegrands[kpzVertices, p4, kpzRegulator]]
```
`Loops.sameUpstream`: the propagator's `"Upstream"` must be exactly the fingerprints of these vertices and this regulator.
:::
