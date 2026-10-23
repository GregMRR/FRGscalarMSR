# Lesson 12 - Oracle and Stages

Two tools that sit beside the pipeline rather than in it. **Oracle** checks the assembled equations numerically, by a route that shares nothing with Propagator, Projection, FrequencyIntegration and Assemble. **Stages** saves any stage to disk and loads it back, refusing a file that is stale.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]];
equations = lessonUpTo["Equations"];
```

## 1. The oracle

**What it does.** At one random point (couplings, `dsCoupling[c]`, k random; each shape a concrete function such as 1/(1 + y²) or e^{−y}; d = 2 and 3; an explicit loop momentum vector), it builds Γ2 + R as a numerical matrix, inverts it, evaluates the loop trace of lesson 07, integrates it over ω with `NIntegrate`, and takes the frequency and momentum projections by finite differences. The result must agree with the symbolic integrand at the same point.

```wl
oracle = RunOracle[equations, "ZeroFlows" -> True];
```

20 comparisons: 2 equations and 8 flows that must vanish, in d = 2 and d = 3. One of them in detail:

```wl
oracle["Results"][[1]]
```

:::why A numerical oracle beside an exact pipeline
Every symbolic module checks itself, but a check written by the same hand as the code can share its blind spots. The oracle answers the question "are these equations right?" by brute force: no atoms, no residues, no direction averages, no compact forms, just numbers. Agreement to 10⁻²⁰ at a random point of a rational function is as strong as evidence gets. Its blind spot is stated on the card (ASSUMPTION `Oracle.scope`): it uses the same vertices and the same loop formula, so it cannot test the normalisation of lesson 07. That is what the Canet comparison (lesson 11, last exercise) is for.
:::

## 2. How the oracle computes

@quote Kernel/Oracle.wl oracleShapes oraclePoint

`MapIndexed[f, list]` calls `f[element, {position}]`, so the i-th shape of your regulator gets the i-th concrete function (`Mod` wraps around if there are more shapes than functions).

@region Kernel/Oracle.wl | (* numeric Gamma2 + R as a polynomial in w | {NIntegrate::izero, NIntegrate::precw}]]];

:::mma Arbitrary precision | `N[x, n]` gives n digits, in software arithmetic that tracks its own precision; `Chop[x, eps]` replaces numbers smaller than eps by 0.
Machine numbers have about 16 digits. The finite differences of the oracle, with h = 10⁻⁵, lose about 10 digits per second derivative, so every number is computed with 40 or more digits: `N[..., wp + 10]` for the matrices, `WorkingPrecision -> wp` for the integral. `Inverse[m]` inverts a numerical (or symbolic) matrix.
```wl
{N[Pi, 30], Chop[10.^-20 + 1. - 1.], Inverse[{{1, 2}, {3, 4}}]}
```
:::

:::mma NIntegrate | `NIntegrate[f, {x, a, b}]` integrates numerically; `PrecisionGoal`, `AccuracyGoal` set the relative and absolute target, `WorkingPrecision` the digits used.
`{wO, -Infinity, -10, 0, 10, Infinity}` splits the real axis at the given points, which helps the integrator. `AccuracyGoal -> wp - 8` lets it stop on integrands that vanish, where a relative goal could never be met. `Quiet[..., {NIntegrate::izero, NIntegrate::precw}]` ([[Check and Quiet]]) silences two expected messages and no others.
:::

The memo keys of `mat`, `g3`, `dsr` and `fLoop` include the working precision ([[Memoization]]): a value computed with 40 digits is never reused where 60 are needed (review E15).

@quote Kernel/Oracle.wl stencil oracleLapNorm oracleLapStencil oracleLaplacianN

The derivatives are central finite differences: `stencil[m, h]` lists the shifts and weights of d^m/dx^m at 0. For the momentum, (∇²_p)ⁿ is applied by nested second differences along every tuple of unit vectors (`Tuples`, lesson 03), and the O(h²) error is removed by Richardson extrapolation: (4 f(h/2) − f(h))/3. `oracleLapNorm` is N_n written again, independently of `lapNorm` in Projection.

`wedgeN[a_?VectorQ, b_?VectorQ] /; Length[a] == Length[b] == 2` evaluates ε only once both vectors are explicit lists ([[PatternTest]]): before that, `wedge[...]` stays symbolic.

## 3. Options and checks of the oracle

@region Kernel/Oracle.wl | Options[RunOracle] | "\"Dimensions\" must be a non-empty subset of {2, 3}"]];

```wl
showFailure[RunOracle[equations, "Dimensions" -> {4}]]
```

The validity conditions are also respected: the point is redrawn (up to 100 times, `Break[]` leaves the `Do` loop) until the poles are causal, because the symbolic residues assume they are.

## 4. Saving and loading stages

@quote Kernel/Stages.wl StageFingerprints SaveStage LoadStage

:::mma Export and Import | `Export[file, expr, "WDX"]` writes any expression in Mathematica's binary format; `Import[file, "WDX"]` reads it back identically.
WDX stores the expression exactly, with symbols and their contexts, so a stage comes back as the same association. `FileExistsQ` tests for a file; `$TemporaryDirectory` is the system's temporary folder.
:::

```wl
file = FileNameJoin[{$TemporaryDirectory, "lessonKPZEquations.wdx"}];
SaveStage[equations, file];
reloaded = LoadStage[file, "Expect" -> {kpzAnsatz, kpzRegulator, kpzSymmetries}];
reloaded === equations
```

`StageFingerprints` walks through every stage embedded in the saved one (recursively: a stage's values that are themselves stages) and collects their fingerprints and upstream lists. `"Expect"` lists the stages you are working with now; each must be among them.

```wl
ansatzWithoutLambda = Block[{Print}, DefineAnsatz[kpzTheory, "Density" -> hbar (dt[h] - nuk lap[h]) - Dk hbar^2, "Couplings" -> {nuk, Dk}]];
showFailure[LoadStage[file, "Expect" -> {ansatzWithoutLambda}]]
```

:::why Persistence with a guard
Heavy models take minutes (shallow water) to hours (large p⁴ truncations), so saving intermediate stages matters. A saved file is also the classic source of silent errors: equations derived from last week's ansatz, loaded into today's notebook. `LoadStage` checks two things: the file was not modified (its own fingerprint), and it was built from the stages you are using now (`"Expect"`). Because every stage embeds its ancestors, the whole history is in the file and can be checked (`design/FRGDESIGN.md`, Q13).
:::

## Exercises

:::exercise A tighter tolerance
Run the oracle with `"Tolerance" -> 10^-15` and only d = 2. Does it pass? What sets the size of the error?
:::

:::solution
```wl
Block[{Print}, RunOracle[equations, "Tolerance" -> 10^-15, "Dimensions" -> {2}]]["Results"][[All, {"Label", "Error"}]]
```
It passes: the errors are of order 10⁻²⁰, the size of the finite-difference error after Richardson extrapolation at the precision used.
:::

:::exercise A modified file
Load the saved stage, change one key, save it under a new name with `Export` (not `SaveStage`), and load it with `LoadStage`. Which check fails?
:::

:::solution
```wl
file2 = FileNameJoin[{$TemporaryDirectory, "lessonKPZEquationsEdited.wdx"}];
Export[file2, Append[reloaded, "Couplings" -> {nuk}], "WDX"];
showFailure[LoadStage[file2]]
```
`Stage.fingerprint`, from `validateStage` (lesson 00): the content no longer matches the fingerprint stored with it. `SaveStage` would have refused to write it in the first place.
:::

:::exercise How much history is in a file
How many fingerprints does `StageFingerprints[equations]` find? Which stages do they belong to?
:::

:::solution
```wl
{Length[StageFingerprints[equations]], Union[Cases[equations, KeyValuePattern["Stage" -> s_String] :> s, {0, Infinity}]]}
```
Every stage of the chain, from Theory to Equations. `Cases` with `KeyValuePattern` ([[Cases and level specifications]]) finds every embedded association with a `"Stage"` key.
:::
