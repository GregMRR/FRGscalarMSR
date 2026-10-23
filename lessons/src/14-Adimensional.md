# Lesson 14 - Adimensional

**Physics step.** Turn the dimensionful equations ∂s(comb) = ∫dᵈq/(2π)ᵈ I(q) into dimensionless β functions. Each state coupling c becomes ĉ = c/[c], with

∂s ĉ = BetaL + ∫daq LoopMeasure(d, aq) BetaNL,  LoopMeasure = S_d aq^{d−1}/(2π)^d,

where BetaL is the local part (canonical scaling and Leibniz terms, not integrated) and BetaNL the loop integrand in the dimensionless momentum aq = q/k. The reference couplings do not get a β: their running gives a linear system for the anomalous dimensions η. The module proves that nothing depends on k or on the references any more.

**Stage.** Input: the Equations stage (lesson 11) and the Dimensions stage (lesson 13), from the same ansatz. Output: the `"Adimensional"` stage, the input of the Julia export.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
equations = lessonUpTo["Equations"];
dims = Block[{Print}, DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}, "Scale" -> k]];
```

## 1. Using it

```wl
adim = AdimensionalEquations[equations, dims];
```

KPZ has one state coupling, λ (ν and D are references). Its flow is purely canonical, because λ is not renormalised (lesson 06):

```wl
ShowBeta[adim, lamk]
```

The anomalous dimensions solve (EtaML + ∫ Measure EtaMNL) η = EtaBL + ∫ Measure EtaBNL, one row per running reference:

```wl
ShowEtaSystem[adim]
```

```wl
ShowCausalPoles[adim]
```

The propagator's validity condition of lesson 05 has become 1 + r_n(aq²) > 0, in minimal dimensionless form: this is what FlowRG will check during the integration (`causal_poles`, lesson 15).

## 2. Units out, by Leibniz

@region Kernel/Adimensional.wl | AdimensionalEquations[equations_, dims_ | gamma2Dim[{a_, b_}] := Expand

The substitution `toUnits` writes each non-reference coupling as ĉ times its dimension, q = k aq, and ∂s ref = −η ref. `Adimensional.scale` compares the regulator's `k` with the `"Scale"` of `DeriveDimensions`: this is where a wrong default `"Scale"` (lesson 13) would be caught.

@region Kernel/Adimensional.wl | (* Leibniz: d_s(comb) = sum_c | dsOf[c_] := Lookup[sol, c, 0];

A flowing combination such as ν·A gives ∂s(νA) = A ∂sν + ν ∂sA: the equations are linear in the unknown ∂s c, and the system must be square and regular. `Array[loopOf, n]` makes the placeholder symbols `loopOf[1], ..., loopOf[n]` for the loop integrals, so that the linear algebra runs on short symbols (lesson 05's abstract atoms, once more). `Variables[poly]` lists the variables of a polynomial.

```wl
{Variables[nuk Dk + lamk^2], Array[loopOf, 3]}
```

## 3. Proving homogeneity

@quote Kernel/Adimensional.wl homVec dimensionVector homClean

`homVec` walks an expression ([[Several definitions of one function]], one rule per head) and returns `{exponent vector, rest}`: the units it contains and what is left. Every term of a sum must carry the same vector; a unit may not appear in an exponent or inside a function argument (a shape argument must reduce to aq²).

:::mma Catch and Throw | `Throw[value, tag]` jumps out of any depth of evaluation to the enclosing `Catch[expr, tag]`, which returns value.
A recursive walk that finds a problem deep inside a sum has to stop and report it. `Throw` does that without every level of the recursion checking a return value. `dimensionVector` catches the reason with the tag `homTag` and turns it into a named check: the library's convention (lesson 00, section 3) is kept at the boundary, while the inner walk uses the cheaper mechanism.
```wl
ui = <|k -> 1, nuk -> 2, Dk -> 3|>;
{Catch[Throw["jumped out", myTag]; "never reached", myTag],
 FRGscalarMSR`Private`homVec[k^2 nuk x, ui, 3],
 Catch[FRGscalarMSR`Private`homVec[k + nuk, ui, 3], FRGscalarMSR`Private`homTag]}
```
:::

@quote Kernel/Adimensional.wl unitInvariantQ

The independent route (review E8): every loop piece, before the units are stripped, is evaluated with every unit multiplied by a random factor, in d = 2 and 3, and must not change. Structural homogeneity and numerical invariance are two different arguments for the same fact.

## 4. The β functions and the η system

@region Kernel/Adimensional.wl | (* a loop piece  alpha x Integral | "linear parts free of aq (not integrated)"]];

`linearPart` (lesson 13) gives the canonical part, ĉ(−x + Σ y_j η_j). The loop pieces become integrands of ∫daq with k^d from the measure.

@quote Kernel/Adimensional.wl LoopMeasure

```wl
{LoopMeasure[3, aq], LoopMeasure[Dim, aq]}
```

@region Kernel/Adimensional.wl | (* eta system: running ref r | "eta system: " <> ToString[Length[running]]

:::mma CoefficientArrays | `CoefficientArrays[{p1, p2, ...}, {x1, x2, ...}]` writes linear (or polynomial) equations as `{b, m, ...}` with p = b + m.x + ...
The η rows are linear in the η's: `CoefficientArrays` returns the constant vector and the matrix in one call (`Normal` turns its sparse arrays into lists). `Length[ca] == 2` means "nothing beyond linear", the check `Adimensional.etaLinear`.
```wl
Normal[CoefficientArrays[{a e1 + b e2 + c, d e1}, {e1, e2}]]
```
:::

## 5. References: running or fixed

@region Kernel/Adimensional.wl | (* references: running (eta from | running = Select[refs

A reference that appears alone in an equation (KPZ: ν and D) is **running**: its η is read from that equation. A reference that appears only multiplied by dimensionless amplitudes (the shallow-water `Dphik aDphik`) is **fixed**: its η is a parameter given at run time, and its row in the η system is a placeholder η = 0, listed as `MustFix` for FlowRG's `fixed_etas` (`design/ADIMENSIONAL.md`).

:::why Which part of the flow is integrated
BetaL contains no momentum integral and BetaNL contains nothing but the integrand. Keeping them apart lets FlowRG evaluate the local part exactly and spend its quadrature only on the loop, and keeps the η system in the form (M_L + ∫M_NL) η = b_L + ∫b_NL that is solved after the integrals are done. The library could integrate nothing here: the shapes are still abstract (lesson 04), and the integrals depend on the state at which FlowRG evaluates them.
:::

## 6. Causal poles in minimal form

@quote Kernel/Adimensional.wl causalPolesCoefficient causalPolesAqExponent causalPolesStrip causalPolesPositiveQ

Each condition e > 0 is brought to p > 0 with p a polynomial (`NumeratorDenominator`, then numerator × denominator, or the numerator times the sign of a numeric denominator), the units are taken out (references, k and aq are positive: an ASSUMPTION line), and the positive numeric content and lowest power of aq are stripped (`GCD`, `Abs`, `Sign`). `Adimensional.causalPolesMinimal` checks exactly that original = stripped factor × minimal (`design/CAUSALITY.md`).

```wl
{NumeratorDenominator[(a + b)/c], GCD[6, 9, 12], Sign[-3]}
```

## 7. The stage

@region Kernel/Adimensional.wl | stage = stampStage[<|"Stage" -> "Adimensional" | PrintCard[stage];

`"Couplings"` are the dimensionless names (`alamk`), `"Etas"` the η symbols, `"Atoms"` the dimensionless atoms `AdimAtom[...]` and their values, `"LoopMomentum"` the symbol `aq`.

```wl
{adim["State"], adim["Couplings"], adim["Etas"], adim["Atoms"]}
```

## Exercises

:::exercise A NonRenormalised reference
Use λ as a reference instead of D. Which check fails, and why must a reference flow?
:::

:::solution
```wl
dims2 = Block[{Print}, DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", lamk -> "l"}, "Scale" -> k]];
showFailure[AdimensionalEquations[equations, dims2]]
```
`Adimensional.referencesFlow`. A reference defines an η through its own flow; λ is declared not to flow, so η_λ would be fixed to 0 by assumption, not computed.
:::

:::exercise The η system with the Litim shape
Substitute the Litim shape r(y) = 1/y − 1 into the first row of EtaMNL (valid for aq < 1). Which atoms simplify?
:::

:::solution
```wl
litim = {rn -> Function[y, 1/y - 1], rd -> Function[y, 1/y - 1]};
Together[adim["EtaMNL"][[1]] /. adim["Atoms"] /. litim]
```
`AdimAtom[h, hbar, 0, 0]` = 1 + r_n(aq²) becomes 1/aq², and `AdimAtom[h, hbar, 0, 1]` = 1 + r_n + aq² r_n′ becomes 0: the Litim shape makes the second atom vanish for aq < 1.
:::

:::exercise The dimensionless measure
Check that `LoopMeasure[d, aq]` integrates to the volume of the unit ball over (2π)^d: ∫₀¹ LoopMeasure daq = V_d/(2π)^d with V_d = π^{d/2}/Γ(d/2 + 1), for d = 3.
:::

:::solution
```wl
{Integrate[LoopMeasure[3, aq], {aq, 0, 1}], Pi^(3/2)/Gamma[5/2]/(2 Pi)^3}
```
Both are 1/(6π²).
:::
