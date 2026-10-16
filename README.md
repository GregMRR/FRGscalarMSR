# FRGscalarMSR

A Mathematica library that turns an ansatz for Γ_k of N MSR scalars into a closed, checked set of Wetterich flow equations,
using the derivative-expansion pipeline. The design is in `FRGDESIGN.md`, the modules are described in `docs/MODULES.md`, and
the conventions in `docs/CONVENTIONS.md`.

**Status: milestone 4 of 6.** Core, Theory, Ansatz, Vertices, Regulator, Propagator, Symmetries, Loops, Projection and FrequencyIntegration are implemented; the 9 shallow-water DE4plus equations reproduce `flowsFAST` exactly. Closure, assembly and oracle follow in milestone 5; KPZ, notebooks and docs in milestone 6. Not yet:
closure, assembly, oracle.

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
