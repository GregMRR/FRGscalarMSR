# Lesson 09 - FrequencyIntegration

**Physics step.** Integrate every projected integrand over the real ω axis, exactly, by residues. Each integrand is a sum of terms ω^n Δ(ω)^{−a} Δ(−ω)^{−b}; the roots of Δ(ω) lie in the upper half plane (the propagator's validity conditions), those of Δ(−ω) in the lower one, so closing the contour upwards picks the residues of Δ(ω) only. The result is a rational function of the coefficients c0, c1, c2 of Δ: no ω, no Δ, no square roots.

**Stage.** Input: the Projection stage. Output: the `"Integrated"` stage, whose `"Flows"` hold, for each monomial, either the symbolic flow (`"Value"`, class Declared) or exact values at random points (`"Samples"`, the other classes).

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
projected = lessonUpTo["Projection"];
```

## 1. Using it

```wl
integrated = IntegrateFrequencies[projected];
```

```wl
ShowFlow[integrated, {h, hbar}, {0, 1}]
```

```wl
ShowFlow[integrated, {h, hbar}, {1, 0}]
```

The first is the flow of ν (class Declared): symbolic, in atoms. The second is protected by the Galilean Ward identity: it is not integrated symbolically, only evaluated exactly at three random points, and all three values are 0 (the closure check of lesson 10 will assert it).

## 2. One integral, by residues

@quote Kernel/FrequencyIntegration.wl FrequencyIntegral

The residue of a pole of order a at z is the (a−1)-th derivative of (ω − z)^a f(ω), at z, divided by (a−1)! ([[Repeated derivatives]]). With degree 2, the roots are z1,2 = iα ± β; the sum of the two residues is symmetric under z1 ↔ z2, hence even in β, and replacing β² by (4c0c2 − c1²)/(4c2²) makes the result rational. The same formula covers the overdamped (β² < 0) and the underdamped (β² > 0) cases.

```wl
{FrequencyIntegral[0, 1, 1, 1], FrequencyIntegral[0, 1, 1, 2]}
```

`FrequencyIntegralC[m]` stands for c_m. A numerical check against Mathematica's own integrator, with Δ(ω) = 2 + 3iω:

```wl
{N[FrequencyIntegral[0, 1, 1, 1] /. {FrequencyIntegralC[0] -> 2, FrequencyIntegralC[1] -> 3}], NIntegrate[1/((2 + 3 I w) (2 - 3 I w)), {w, -Infinity, Infinity}]}
```

:::mma Numerator, Denominator and radicals | `Numerator[r]` and `Denominator[r]` split a rational expression; a square root is `Power[x, 1/2]`, so `FreeQ[r, Power[_, _Rational]]` tests for radicals.
`r` is made even in β by multiplying numerator and denominator by the denominator with β → −β; then both contain only even powers of β, and β → √x followed by x → its value leaves no root. `FrequencyIntegration.rational` asserts that no `Power[_, _Rational]` is left (CLAUDE.md, Conventions: no radicals).
```wl
{Numerator[(a + b)/c2], Denominator[(a + b)/c2], FullForm[Sqrt[x]], FreeQ[Sqrt[x], Power[_, _Rational]]}
```
:::

A divergent integral aborts:

```wl
showFailure[FrequencyIntegral[2, 1, 1, 1]]
```

:::mma Memoization of a function | `f[n_, a_] := f[n, a] = ...` stores each computed value as a new definition of f, so the next call with the same arguments returns it at once.
This is [[Memoization]] applied to a top-level function: `FrequencyIntegral[0, 1, 1, 1]` is computed once per session, and every integrand that needs it afterwards reuses the stored rule. The stored values are visible in `DownValues[FrequencyIntegral]` ([[OwnValues and DownValues]]) as the rules without patterns; the card counts them.
```wl
Length[Select[DownValues[FrequencyIntegral], FreeQ[#, Pattern] &]]
```
:::

:::why A memoised formula instead of a lookup table
The original notebook used a hand-made table of integrals, which had to be extended whenever a new combination (n, a, b) appeared, and silently missed any combination that was not in it. Here each integral is computed from one residue formula the first time it is needed and remembered: nothing is skipped, and nothing depends on the order of rules. The memo is a cache of a pure function of (n, a, b, deg), so it is the one kind of global state that cannot make results depend on history.
:::

## 3. Integrating an integrand

@quote Kernel/FrequencyIntegration.wl integrate

Δ(ω) → 1/X and Δ(−ω) → 1/Y turn the integrand into a polynomial in ω, X, Y with coefficients free of ω.

:::mma CoefficientRules | `CoefficientRules[poly, {x, y, ...}]` lists the monomials of a polynomial as `{exponents} -> coefficient`.
Each rule `{n, a, b} -> coefficient` is one term ω^n X^a Y^b, so the integral is `Total` of coefficient × `FrequencyIntegral[n, a, b, deg]`. This dispatches on the structure of the integrand, term by term, without pattern rules that could miss a term.
```wl
CoefficientRules[3 w^2 x + 5 y, {w, x, y}]
```
:::

## 4. Sampled flows

@quote Kernel/FrequencyIntegration.wl samplePoint sampleValue validPointQ validSample

A flow that must vanish (Protected, DeclaredNumeric, Causality) is not integrated symbolically: an exact zero is zero at every point, so three random rational points prove it with the reliability of `ExactlyZeroQ` (lesson 00, section 5). Every atom, `dsR` and `dsCoupling` gets a random value, and the coefficients of Δ are computed *from* these atoms, so that Δ = det B still holds at the point.

A Neglected flow is reported as a size, so its point must be physical: q² > 0 and the validity conditions of the propagator satisfied. `validSample` redraws points until one is valid.

:::mma Return from a loop | `Do[..., {j, 0, 199}]` is a loop; `Return[value, Module]` leaves the enclosing `Module` at once with that value.
`validSample` tries up to 200 seeds and returns the first valid sample; if the loop ends, the next line, a failing `FRGAssert`, runs. `And @@ list` ([[Apply]]) combines a list of conditions into one.
:::

## 5. The constructor

@region Kernel/FrequencyIntegration.wl | IntegrateFrequencies[projected_] | PrintCard[stage];

`Join[KeyDrop[r, "Integrand"], out, <|"Seconds" -> ...|>]` replaces the integrand of each record by its result and adds the time it took; `Join` merges associations ([[Association]]).

:::why Symbolic where needed, sampled where enough
Only the Declared flows become equations, so only they are integrated symbolically. The others only have to be shown zero, or measured; doing that at exact random points is faster by orders of magnitude and loses nothing: a zero test at random rationals is as reliable as a symbolic simplification, and much more so than `Simplify` on large expressions. The card says how it was done (the INFO line), so that a reader knows what "0" means.
:::

## Exercises

:::exercise A residue by hand
For degree 1, Δ(ω) = c0 + i c1 ω has its root at ω = i c0/c1. Compute ∫ dω / (Δ(ω) Δ(−ω)) by residues and compare with `FrequencyIntegral[0, 1, 1, 1]`.
:::

:::solution
Δ(ω)Δ(−ω) = c0² + c1² ω²: the residue at ω = i c0/c1 of 1/(c1²(ω − i c0/c1)(ω + i c0/c1)) is 1/(c1² · 2i c0/c1) = 1/(2i c0 c1); times 2πi gives π/(c0 c1).
```wl
FrequencyIntegral[0, 1, 1, 1]
```
:::

:::exercise Degree 2 against NIntegrate
Take Δ(ω) = 2 + 3(iω) + (iω)² (roots in the upper half plane). Compare `FrequencyIntegral[2, 2, 2, 2]` with a numerical integral of ω²/(Δ(ω)² Δ(−ω)²).
:::

:::solution
```wl
dd[x_] := 2 + 3 I x - x^2;
{N[FrequencyIntegral[2, 2, 2, 2] /. {FrequencyIntegralC[0] -> 2, FrequencyIntegralC[1] -> 3, FrequencyIntegralC[2] -> 1}],
 NIntegrate[w^2/(dd[w]^2 dd[-w]^2), {w, -Infinity, Infinity}]}
```
:::

:::exercise Which flows are sampled
List the classes of the records of `integrated["Flows"]` and how many have a `"Value"` or `"Samples"`.
:::

:::solution
```wl
Tally[{#Class, If[KeyExistsQ[#, "Value"], "Value", "Samples"]} & /@ integrated["Flows"]]
```
`Tally` counts the distinct pairs: only the two Declared monomials (ν and D) are symbolic.
:::
