# Lesson 10 - Closure

**Physics step.** Decide whether the flow closes on the ansatz. Every flow that a symmetry, causality, or a number in the ansatz forbids must vanish; every neglected flow is reported with a sample value; every declared coefficient gets its equation ∂s(combination) = RHS, where the Γ_k coefficient is split into a pure number and the flowing combination of couplings.

**Stage.** Input: the Integrated stage and the Symmetries stage it was projected with. Output: the `"Closure"` stage, with the `"Equations"` and a `"Report"` of every monomial.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
integrated = lessonUpTo["FrequencyIntegration"]; symmetries = kpzSymmetries;
```

## 1. Using it

```wl
closure = CheckClosure[integrated, symmetries];
```

Three PASS lines say that the protected, numeric and causality flows vanish exactly. The WARN lines are the neglected monomials that the flow generates; their values are a sample at one random physical point, not a bound. The report puts it all in one table:

```wl
ShowClosureReport[closure]
```

## 2. The flows that must vanish

@region Kernel/Closure.wl | CheckClosure[integrated_, symmetries_] | {c, {{"Protected", "Closure.protected"}

`! MatchQ[Lookup[#, "Samples"], {0 ..}]` selects the records whose samples are not all exactly 0 ([[MatchQ and repeated patterns]]). The `Do` runs over a list of pairs (class, check name), so the three checks share one piece of code. The hint names each offending monomial with the reason of its tag, because a failure here means the tag is wrong or the ansatz breaks the symmetry, and the reason is what you need to decide which.

Tagging ω² of Γ_{h̄h̄} as Protected, which it is not:

```wl
symmetries2 = Block[{Print}, DeclareSymmetries[kpzAnsatz,
     Protected[{hbar, h}, "Frequency" -> 1, "Momentum" -> 0, "Reason" -> "Galilean Ward identity"],
     Protected[{h, hbar}, "Frequency" -> 0, "Momentum" -> 0, "Reason" -> "shift"],
     Protected[{hbar, hbar}, "Frequency" -> 2, "Momentum" -> 0, "Reason" -> "wrong on purpose"],
     Neglected[_, "Reason" -> "beyond"], "TruncationOrders" -> <|"Frequency" -> 2, "Momentum" -> 1|>]];
integrated2 = Block[{Print}, IntegrateFrequencies[ProjectOnMonomials[kpzLoops, symmetries2]]];
showFailure[CheckClosure[integrated2, symmetries2]]
```

:::why A false symmetry is caught, a true one is only assumed
The library cannot derive Ward identities; you declare them. What it can do is test the consequence: a declared symmetry forbids certain flows, and those flows are computed and required to be exactly 0. A wrong tag therefore aborts here, with its reason. A correct tag passes, which supports it but does not prove the symmetry (the flows might vanish for another reason, or only at this truncation), and the card keeps the ASSUMPTION line "Ward identities are assumed, not derived".
:::

## 3. The neglected flows

@region Kernel/Closure.wl | (* neglected flows: reported, not asserted *) | not a bound (" <> #Reason <> ")"]] & /@ neg];

@quote Kernel/Closure.wl monomialString numString recordString

:::mma N and CForm | `N[x]` gives a machine-precision approximation of an exact number; `CForm[x]` writes an expression as C code, here a compact number string.
The samples are exact rationals; the card shows them as decimals. `N[x, 4]` asks for 4 digits.
```wl
{N[2/3], N[2/3, 4], ToString[CForm[N[2/3]]]}
```
:::

## 4. Declared monomials become equations

@quote Kernel/Closure.wl splitCoefficient

:::mma FactorTermsList | `FactorTermsList[poly]` pulls out the numeric factors: `{number, ..., rest}`.
`splitCoefficient` multiplies all numeric parts together (`Times @@ Most[f]`, [[Apply]]) and keeps the rest: −2 Dk gives normalisation −2 and combination Dk, −λ/2 gives −1/2 and λ.
```wl
{FactorTermsList[-2 Dk], FactorTermsList[-lamk/2], FactorTermsList[6 a b + 4 a]}
```
:::

@region Kernel/Closure.wl | (* declared monomials: d_s(comb) = Value / norm *) | ToString[Length[eqs]] <> " equations d_s(comb)

The equation for a coefficient written as −2 Dk is ∂s(Dk) = flow/(−2): the number is fixed, only the couplings flow (decision Q10 of `design/FRGDESIGN.md`).

```wl
KeyTake[#, {"GammaK", "Normalisation", "Flowing"}] & /@ closure["Equations"]
```

`KeyTake[assoc, keys]` keeps only the listed keys, in that order: the opposite of `KeyDrop` ([[Association]]).

## 5. The report

@quote Kernel/Closure.wl outcome ShowClosureReport

:::mma Switch | `Switch[expr, form1, value1, form2, value2, ...]` returns the value of the first form (a pattern) that expr matches.
Like [[Which]], but the tests are patterns matched against one expression: `_` at the end is the default.
```wl
{Switch[2, 1, "one", 2, "two", _, "other"], Switch[x, _Integer, "integer", _, "other"]}
```
:::

## Exercises

:::exercise The wrong Symmetries stage
Pass `integrated` (projected with `symmetries`) together with `symmetries2` to `CheckClosure`. Which check fails, and why is it needed even though both stages come from the same ansatz?
:::

:::solution
```wl
showFailure[CheckClosure[integrated, symmetries2]]
```
`Closure.sameSymmetries`. The classes decide which flows are asserted to vanish; checking against a classification other than the one used for the projection would assert the wrong things. The fingerprint makes the mismatch impossible to miss.
:::

:::exercise Read the equations
Print the right-hand side of the ν equation with `LeafCount`, and list the atoms it contains.
:::

:::solution
```wl
nuEq = SelectFirst[closure["Equations"], #Flowing === nuk &];
{LeafCount[nuEq["RHS"]], Union[Cases[nuEq["RHS"], _GammaAtom[_] | _dsR | DeltaCoefficient[_][_] | Derivative[_][DeltaCoefficient[_]][_], Infinity]]}
```
Still in atoms: `GammaAtom`, `dsR`, `DeltaCoefficient` and its q²-derivatives. They get their values in the next lesson.
:::

:::exercise How big is a neglected term?
Run `CheckClosure` again: do the WARN values change? Why is that the right behaviour?
:::

:::solution
```wl
Block[{Print}, CheckClosure[integrated, symmetries]]["Report"][[{10, 14}, {"Pair", "Samples"}]]
```
No: the samples were drawn in lesson 09 with fixed seeds, and Closure only reads them. Reproducible random points make every card the same from one run to the next.
:::
