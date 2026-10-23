# Lesson 11 - Assemble

**Physics step.** Write each flow equation ∂s(combination) = ∫dᵈq/(2π)ᵈ (integrand) in two forms. The **compact** form uses scale-stripped atoms (`ScaleAtom`, e.g. 1 + r_n(y)) and the abstract ∂sR entries: short, readable, close to what one writes by hand. The **explicit** form substitutes everything: couplings, `dsCoupling[c]`, the shapes and their derivatives, `Dim`. That is the q-integrand handed to the dimensionless step and to FlowRG.

**Stage.** Input: the Closure stage. Output: the `"Equations"` stage, the end of the dimensionful pipeline.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
closure = lessonUpTo["Closure"];
```

## 1. Using it

```wl
equations = AssembleEquations[closure];
```

```wl
ShowEquation[equations, nuk]
```

`ShowEquation` prints the compact form, the atoms it uses, the ∂sR entries and where the equation is valid (the causal-pole condition of lesson 05, now in explicit form). The explicit integrand of the D equation:

```wl
EquationOf[equations, Dk]["Integrand"]
```

## 2. Scale-stripped atoms

@quote Kernel/Assemble.wl stripScale atomSystem

Each propagator atom v(s), s = q², is written v = c sⁿ X, where c is the number and couplings common to all its terms, n the smallest power of s, and X what is left. For R_{h̄h}: ν s (1 + r_n) gives c = ν, n = 1, X = 1 + r_n(y).

```wl
FRGscalarMSR`Private`stripScale[nuk s + nuk s rn[s/k^2], s, {nuk, lamk, Dk}, {rn, rd}]
```

When X is a sum it becomes the atom `ScaleAtom[a, b, m, 0]`; its s-derivatives become `ScaleAtom[a, b, m, j]`. `atomSystem` returns three things: rules from `GammaAtom` to atom forms, the atom values, and the forms, so that `Assemble.atomSync` can check every atom and derivative against the propagator.

```wl
{equations["Atoms"], equations["AtomForms"]}
```

:::why Generate the atoms the hand calculation would use
The original script had hand-written `toAtoms`/`atomVals`: names like xE = A + r and xE1 = A + r + y r′ chosen by a person for readability. Here the same atoms are derived from the propagator by one rule (strip the common scale, name what is left), so a new ansatz gets them for free and they cannot drift from the propagator they stand for (`Assemble.atomSync`). In shallow water the rule reproduces the script's atoms exactly (`docs/MODULES.md`, Assemble).
:::

## 3. From the closure to the equations

@region Kernel/Assemble.wl | AssembleEquations[closure_] | "d_s R linear in dsCoupling[c]"]];

:::mma Replace and ReplaceAll | `Replace[e, rule]` tries the rule on the whole of e only; `e /. rule` (`ReplaceAll`) tries every part.
`Replace[reg["LoopMomentum"], None -> Symbol["Global`q"]]` replaces the value only if it *is* `None`. With `/.`, a `None` buried anywhere inside would be replaced too. `Replace[e, rules, {1}]`, seen in lesson 02's `splitTerm`, applies at level 1 only.
```wl
{Replace[None, None -> 1], Replace[{None}, None -> 1], {None} /. None -> 1}
```
:::

The Δ coefficients and their q²-derivatives are first expressed in `GammaAtom` atoms (`deltaRules`), then the atoms in scale-stripped form. The loop is proportional to ∂sR, so each right-hand side must be linear and homogeneous in the `dsR` entries: `Assemble.linearInDsR` tests that doubling every `dsR` doubles the result, and that setting them to 0 gives 0.

@region Kernel/Assemble.wl | eqs = MapThread[Function[{r, rw} | "Seconds" -> AbsoluteTime[] - t1|>]]], {closure["Equations"], raw}];

`MapThread[f, {list1, list2}]` ([[Map at a level]] for the cousin `Map`) runs over the equations and their raw right-hand sides in parallel.

## 4. Choosing a compact form

@quote Kernel/Assemble.wl collectNum compactFrac compactify atomLeaves

:::mma Collect and Factor | `Collect[e, vars, f]` groups e by powers of vars and applies f to each coefficient; `Factor` factors a polynomial.
`Collect[a, dsList, Together]` writes the right-hand side as Σ (coefficient) × `dsR[...]` with each coefficient over a common denominator. `compactify` tries several such forms (per `dsR` coefficient or as one fraction, atoms collected in two orders) and keeps the one with the smallest explicit `LeafCount`.
```wl
{Collect[a^2 x^2 + 2 a b x^2 + b^2 x^2 + c x, x], Collect[a^2 x^2 + 2 a b x^2 + b^2 x^2 + c x, x, Factor], Factor[x^2 - 1]}
```
:::

:::mma TimeConstrained | `TimeConstrained[expr, t, alt]` evaluates expr for at most t seconds, and returns alt if time runs out.
For small results only (explicit `LeafCount` below 15000), `compactify` also tries `Simplify` on each coefficient, limited to 10 seconds, falling back to the coefficient unchanged.
```wl
TimeConstrained[Pause[3]; "finished", 0.5, "timed out"]
```
:::

:::why The one place where Simplify is allowed
CLAUDE.md forbids `Simplify` on large expressions and as a zero test (lesson 00, section 5). Here it is used for presentation only, on small expressions, with a time limit, as one candidate among several, and the winner is checked exactly against the canonical form (`Assemble.compactExact`, with `ExactlyZeroQ`). If `Simplify` is slow, or makes things worse, nothing is lost. This is how the library uses heuristics in general: never trusted, always checked, never needed for correctness.
:::

`atomLeaves` counts each `ScaleAtom` and `dsR` as a single leaf, so the compact `LeafCount` on the card is comparable with the script's hand-made compact forms.

## 5. The explicit form and the final checks

The explicit integrand is `cpt /. atomVals /. dsrRules`: atom values, then the derived ∂sR of lesson 04, substituted only now. `Assemble.noIntermediates` requires that nothing internal survives: no `GammaAtom`, `ScaleAtom`, `Δ`, `DeltaCoefficient`, `dsR`, `LoopMomentumSquared`, and no symbol from a `Private` context (`StringContainsQ[Context[x], "Private`"]`).

@region Kernel/Assemble.wl | validity = prop["ValidityConditions"] | PrintCard[stage];

The stage records everything a reader of the equations needs: the atoms, the ∂sR rules, the loop momentum and scale symbols, the shapes, the couplings, the dimension and the validity conditions.

@quote Kernel/Assemble.wl EquationOf ShowEquation

## Exercises

:::exercise The size of the result
For both equations, print `#LeafCount` and `#CompactForm` (which candidate won). How much smaller is the compact form than the explicit integrand?
:::

:::solution
```wl
{#Flowing, #CompactForm, #LeafCount} & /@ equations["Equations"]
```
:::

:::exercise Linearity in dsCoupling
The explicit integrand is linear in `dsCoupling[nuk]` and `dsCoupling[Dk]`. Extract the coefficient of `dsCoupling[Dk]` in the D equation.
:::

:::solution
```wl
Together[Coefficient[EquationOf[equations, Dk]["Integrand"], dsCoupling[Dk]]]
```
Downstream, `dsCoupling[Dk]` becomes −η_D D, and this coefficient multiplies η_D in the linear system for the anomalous dimensions (lesson 14).
:::

:::exercise Canet's KPZ flow
With the Litim cutoff r(y) = (1/y − 1) θ(1 − y) in both sectors, and ∫dᵈq/(2π)ᵈ f = v_d k^d/2 ∫₀¹ y^{d/2−1} f dy, integrate the ν and D equations, solve for η_ν and η_D, and check that d − 2 − η_D + 3η_ν equals the ratio of Canet's Eq. (20) (arXiv:2509.01472), with g = v_d k^{d−2} λ² D/ν³. This is the test that fixes the loop normalisation of lesson 07.
:::

:::solution
`litimIntegral` replaces q by k√y, inserts the Litim shape (valid for y < 1, where θ = 1), and integrates each power of y from 0 to 1.
```wl
litim = {rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};
litimIntegral[f_] := Module[{gg = Expand[Together[f /. q -> k Sqrt[y] /. litim] y^(d/2 - 1)], terms},
  terms = If[Head[gg] === Plus, List @@ gg, {gg}];
  vd k^d/2 Total[(# /. y -> 1)/(Exponent[#, y] + 1) & /@ terms]];
flowOf[c_] := litimIntegral[EquationOf[equations, c]["Integrand"] /. Dim -> d];
```
```wl
etas = First[Solve[{-etanu nuk == (flowOf[nuk] /. {dsCoupling[nuk] -> -etanu nuk, dsCoupling[Dk] -> -etaD Dk}),
     -etaD Dk == (flowOf[Dk] /. {dsCoupling[nuk] -> -etanu nuk, dsCoupling[Dk] -> -etaD Dk})}, {etanu, etaD}]];
bracket = Together[(d - 2 - etaD + 3 etanu) /. etas /. lamk -> Sqrt[g nuk^3/(vd k^(d - 2) Dk)]]
```
```wl
canet = (2 d^6 + 4 d^5 - d^4 (7 g + 8) - 2 d^3 (3 g + 8) + d^2 g (3 g + 20) - d (g - 8) g - 2 g^2)/
   (2 d^5 + 8 d^4 + d^3 (8 - 3 g) - 2 d^2 g + d g (g + 8) - g^2);
Together[bracket - canet]
```
0, for every d and g. The same computation is `tests/acceptance/KPZ.wl` and the last section of `examples/KPZ.nb`.
:::
