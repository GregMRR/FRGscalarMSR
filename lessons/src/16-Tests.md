# Lesson 16 - Tests

The library's checks run inside the pipeline, on your model, every time. The test suite checks the library itself: that each module gives known results on known models, and that each check fires when it should. This lesson reads the harness `tests/runTests.wl` and the way tests are written, and points to the benchmarks.

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"];
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]]
```

## 1. Running the suite

From a shell: `math -script tests/runTests.wl` runs every file (exit code 0 means all pass); `--only TheoryTest` runs the files whose name contains the string; `--verbose` shows the cards and FAIL blocks; `--require-julia` turns the SKIPs of the Julia tests into failures. From a notebook, a separate kernel can run it:

:::mma RunProcess | `RunProcess[{"cmd", "arg", ...}, "StandardOutput"]` runs an external program and returns what it printed.
The tests run in a fresh kernel, as a script, so they cannot disturb this notebook's kernel, and the `Exit[...]` at the end of the harness ends that process only.
```wl
RunProcess[{"math", "-script", FileNameJoin[{$FRGRoot, "tests", "runTests.wl"}], "--only", "TheoryTest"}, "StandardOutput"]
```
:::

## 2. The harness

@region tests/runTests.wl | $HistoryLength = 0; | $FRGExitOnFail = False;

`$FRGExitOnFail = False` turns every failed check into an `Abort[]` instead of `Exit[2]` ([[Abort, CheckAbort, Exit]]), so the harness can catch failures and go on. `$HistoryLength = 0` stops the kernel from remembering every output (`Out[n]`), which matters with million-leaf expressions.

@quote tests/runTests.wl quiet record skip test expectFail

- `quiet[e]` hides the cards while e runs, unless `--verbose` (`Block[{Print}, e]`, as `lessonUpTo` does).
- `test[name, actual, expected]` evaluates `actual` under `CheckAbort` and compares it with `expected` by `===` ([[SameQ and Equal]]): exact structural identity, no tolerance.
- `expectFail[name, expr, check]` passes only if `expr` aborts, *and* the check that stopped it is the named one (`$FRGLastFailure`). An abort from another check is a failure: a test of `Theory.N` that fails on `Theory.scalars.reserved` is not testing what it says.
- `skip[name, why]` records a test that cannot run here (no Julia, no FlowRG).

Both take their arguments held (`HoldRest`), so the expression is evaluated inside `CheckAbort`, not before the call ([[Hold attributes]]).

@region tests/runTests.wl | testFiles = Join[ | Exit[If[nFail == 0

The test files are plain `.wl` files read with `Get` ([[Get]]) one after the other, in alphabetical order, in one kernel: they share global variables. Heavy setups (`tests/reference/ShallowWaterSetup.wl`, about 85 s) are therefore built once, behind `If[! ValueQ[...]]`, and reused by every file that needs them (CLAUDE.md, Gotchas).

## 3. A unit test file

@file tests/unit/TheoryTest.wl

Each check of lesson 01 has an `expectFail`. That is the rule of the library (CLAUDE.md, Where to change X): a new check is a named `FRGAssert`/`assertCheck`, and it comes with an `expectFail[..., "Module.check"]` test. A check that is never seen failing might never fail.

:::why Tests that name the failure
`expectFail` compares the *name* of the check that aborted, which is possible only because every check has a unique name and records it in `$FRGLastFailure` (lesson 00, section 3). The same names appear on the cards, in the FAIL blocks, in `docs/MODULES.md` and in the tests: one vocabulary for the user, the documentation and the test suite. In the notebooks of these lessons, `showFailure` plays the same role as `expectFail`: it returns the name, so a cell can be checked by eye.
:::

## 4. Kinds of tests

- **Unit tests** (`tests/unit/<Module>Test.wl`): one file per module, small inputs, every check made to fail at least once.
- **Acceptance tests** (`tests/acceptance/`): whole models against independent results. Shallow water against the original hand-derived script and the old notebooks (`tests/reference/legacy/`); KPZ against Canet's flow (`KPZ.wl`, the computation of lesson 11's last exercise); the Julia round trip (`tests/reference/JuliaRoundTrip.wl`, lesson 15).
- **Regression tests**: `ProjectionRegression.wl` keeps the old projection operator (`tests/reference/OldProjection.wl`) only to require that the new one agrees with it (lesson 08).
- **Benchmarks** (`benchmarks/`): published results reproduced with the library, one folder each, with the paper's equations (`Reference.wl`), the pipeline (`Pipeline.wl`), the comparisons (`Checks.wl`), a headless driver and a notebook. Each has a suite hook in `tests/acceptance/` that only wraps its checks in `test[...]`:

@file tests/acceptance/KuramotoSivashinsky.wl

The last test is worth copying in any comparison: it perturbs the reference by 1% and requires the comparison to *fail*. A comparison that cannot fail tests nothing.

## Exercises

:::exercise Your own expectFail
Without the harness, write the equivalent of `expectFail["N = 3 aborts", DefineTheory["Scalars" -> {f1, f2, f3}, "Dimension" -> Dim], "Theory.N"]` with `showFailure`, returning `True` or `False`.
:::

:::solution
```wl
showFailure[DefineTheory["Scalars" -> {f1, f2, f3}, "Dimension" -> Dim]] === "Theory.N"
```
:::

:::exercise A test that tests the wrong thing
You want to test the reserved-name check, and write the scalars as `{k, f2, f3}`. Compare the result of `showFailure` with `"Theory.scalars.reserved"`. Why is it `False`, and why is that the right answer?
:::

:::solution
```wl
showFailure[DefineTheory["Scalars" -> {k, f2, f3}, "Dimension" -> Dim]] === "Theory.scalars.reserved"
```
`False`: `DefineTheory` checks the number of scalars before their names, so `Theory.N` fires first. An `expectFail` written this way would fail and reveal that the test does not exercise the check it is named after; a test that only asked "did it abort?" would pass and test nothing about reserved names.
:::

:::exercise The Julia round trip
Run the KPZ acceptance files with `--only KPZ` in a separate kernel. Which tests are about Julia, and do they run or SKIP on this machine?
:::

:::solution
```wl
StringCases[RunProcess[{"math", "-script", FileNameJoin[{$FRGRoot, "tests", "runTests.wl"}], "--only", "KPZ"}, "StandardOutput"], ("SKIP" | "===") ~~ Except["\n"] ..]
```
With Julia and FlowRG installed, the round-trip tests run and count among the passed ones; otherwise they print SKIP with the reason.
:::
