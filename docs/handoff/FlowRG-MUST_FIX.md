# Handoff — two FlowRG interface extensions requested by FRGscalarMSR

For a fresh agent working in `~/Physics/Sources/FlowRG.jl`. Read that repo's `CLAUDE.md` first (golden rules, traps) and
`docs/equations-interface.md` (the contract). Nothing here changes behaviour for existing equations files.

## Why

FRGscalarMSR (`~/Physics/Sources/FRGscalarMSR`, Mathematica) is getting a `JuliaExport` module that writes FlowRG
equations files automatically from a derived set of flow equations. Two of its outputs do not fit the current interface:

1. **Etas fixed by hand.** Some anomalous dimensions are not determined by an equation: the physicist fixes them at run
   time (e.g. η_φ chosen so that the energy injection of a large-scale forcing stays finite). The generated file keeps that
   eta in `ETAS` with a **placeholder row** η = 0 (identity row of `etaMl`, nothing in `etaMnl`, `etabl`, `etabnl`), and
   declares it in a new optional constant:

   ```julia
   const MUST_FIX = (:phi,)       # subset of ETAS; each must be given in fixed_etas
   ```

   Running such a file without fixing the eta would silently use η = 0. FlowRG must refuse that.

2. **Second derivatives of the shapes.** Flows of momentum-derivative couplings (g2, D2φ in the shallow-water DE4plus set)
   contain r″(y). The generated file writes them as `D2armh(aq^2)`, next to `armh(aq^2)` and `Darmh(aq^2)`.

## Task 1 — `MUST_FIX`

- `src/load.jl`, `build_model`: if `isdefined(E, :MUST_FIX)`,
  - assert `Tuple(E.MUST_FIX) ⊆ ETAS`, else `ArgumentError` naming the unknown entries and the file's ETAS;
  - assert every entry is a key of `fixed_etas`, else `ArgumentError` naming the missing etas and the call to write, e.g.
    `load_model(path; fixed_etas = (phi = <value>,))`.
  - A file without `MUST_FIX` behaves exactly as today.
- `docs/equations-interface.md`: add `MUST_FIX` to "What the file must define" as optional, with its meaning (placeholder
  row, value supplied at load time) and the rule that its names must be in `ETAS`. Keep `README.md` in sync (golden rule 3).
- Tests: a toy file `test/toy_must_fix.jl` (copy the style of `toy_equations.jl`) with one eta in `MUST_FIX`:
  - `load_model` without `fixed_etas` throws `ArgumentError` mentioning the eta;
  - with `fixed_etas = (phi = 0.3,)` it loads and `etas(m, u)` returns 0.3 for it;
  - a `MUST_FIX` entry not in `ETAS` throws.

## Task 2 — `D2<shape>`

- `src/regulators.jl`: add the second derivative to `Regulator` (a `D2shape` field, y → r″(y)), for every shipped
  regulator. Keep the existing 4-argument constructor working (e.g. `D2shape = nothing`, and binding a `D2` name to a
  regulator without it throws a clear error when the file uses it). Wetterich: r″(y) = e^y (e^y + 1)/(e^y − 1)³, written
  overflow-safe like the existing `wetterich_Dshape`.
- `src/load.jl`, `bind_regulators!`: bind `D2<name>` for every name in `SHAPES`, next to `D<name>`, as `const` (trap:
  regulator bindings must stay `const`).
- `docs/equations-interface.md` rule 4: shapes, first and second derivatives (`D2armh(aq^2)`), all with respect to y.
- Tests: a toy file using `D2` checked against a finite difference of `D` for Wetterich; a regulator without `D2shape`
  used by a file that calls `D2…` gives a clear error.

**Open physics question — ask the user before implementing it for Litim.** Litim's r′(y) = −1/y² θ(1 − y) jumps from −1
to 0 at y = 1, so r″(y) = 2/y³ θ(1 − y) + δ(y − 1). A pointwise `D2` cannot represent the δ, and the quadrature would
silently drop its contribution. Options to put to the user: (a) Litim gets no `D2shape`, so files with r″ refuse Litim;
(b) a smoothed Litim; (c) FlowRG adds the boundary term explicitly (requires knowing the integrand's coefficient of r″ at
aq = 1, i.e. an interface change). Do not choose silently.

## Constraints

- Follow FlowRG's `CLAUDE.md` traps: `Base.invokelatest` around module bindings, fixed-eta rows masked with `ifelse`,
  static sizes through `Val`.
- Do not edit generated equations files or any Mathematica notebook.
- Do not loosen test tolerances. Run the full test suite (`julia --project -e 'using Pkg; Pkg.test()'`) at the end.
