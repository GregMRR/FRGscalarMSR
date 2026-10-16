# FRGscalarMSR

A Mathematica library that turns an ansatz for Γ_k of N MSR scalars into a closed, checked set of Wetterich flow equations,
using the derivative-expansion pipeline. The design is in `FRGDESIGN.md`, the modules are described in `docs/MODULES.md`, and
the conventions in `docs/CONVENTIONS.md`.

**Status: milestone 5 of 6.** The whole pipeline is implemented, from `DefineTheory` to `AssembleEquations`, together with the
numerical oracle (`RunOracle`) and `SaveStage`/`LoadStage`. The 9 shallow-water DE4plus equations reproduce `flowsFAST`
exactly, and the oracle agrees with them in d = 2 and 3. Milestone 6 adds KPZ, the example notebooks and the documentation.

```mathematica
Get["~/sources/FRGscalarMSR/FRGscalarMSR.wl"]
theory = DefineTheory["Scalars" -> {H, Psi}, "Dimension" -> Dim]        (* Hbar, Psibar are created *)
ansatz = DefineAnsatz[theory, "Density" -> Hbar (dt[H] - nk Ak lap[H]) + ..., "Couplings" -> {nk, Ak, ...}]
```

Tests: `math -script tests/runTests.wl [--verbose]` (exit code 0 means every test passed).

Assumptions to keep in mind:
- Ward identities are assumed (declared by the user), not derived.
- The loop normalisation (prefactor, sign, 2π, Γ⁽⁴⁾ tadpole) will be hand-derived in `docs/MODULES.md`; it is not tested
  independently.
