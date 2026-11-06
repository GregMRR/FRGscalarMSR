# Dimensionless flows and the FlowRG file

Part 2 of the derivative-expansion pipeline: dimensions derived from the ansatz, dimensionless β functions and anomalous
dimensions, and a FlowRG equations file. It starts from the Equations stage of [part 1](flow-equations.md).

## Steps

1. **Dimensions.** `DeriveDimensions[ansatz, "References" -> {…}, "Dimensionless" -> {…}]` derives the dimension of every
   coupling, field and of time from the ansatz (Γ_k is dimensionless, [x] = 1/k), written as k^x ref₁^y₁ … in the basis of k and
   **N+1 reference couplings** that you choose (N for the fields, one for time; for shallow water `nk`, `Dhk`, `Dphik`). Each
   reference names its anomalous dimension: `nk -> "nu"` gives η_ν (symbol `etanu`), with ∂s nk = −η_ν nk.
   - Declare in `"Dimensionless"` the amplitudes that share a monomial with another coupling (`Ak` in `nk Ak`, `aDphik` in
     `Dphik aDphik`): one monomial cannot split the dimension between them. If you forget, the error names them.
   - Checked: exactly N+1 free units, independent references, every monomial consistent, and the action invariant under the
     derived rescaling.
   - `ShowDimensions[dims]` prints coupling | monomials | dimension | dimensionless name | linear part of its flow, then the field
     dimensions, [ω] and z (e.g. [ω] = k² nk, z = 2 − η_ν). `CheckDimensions[dims, <|gk -> k Dphik^(1/2) Dhk^(-1/2) nk, …|>]`
     compares your hand values exactly.
2. **Dimensionless flows.** `AdimensionalEquations[equations, dims, "State" -> {gk, lk, …}]`: the state couplings (dimensionful
   names, in the order you want in Julia) become `agk`, `alk`, … (declared dimensionless couplings keep their name), with
   ∂s ĉ = BetaL + ∫daq LoopMeasure BetaNL: BetaL is the canonical part ĉ(−x + Σ y_j η_j) plus local terms, **not integrated**;
   BetaNL is the loop integrand, with LoopMeasure = S_d aq^{d−1}/(2π)^d and q = k aq. The η's solve a linear system with the same
   split. Every β is checked to be free of k and of the references.
   - A reference that appears **alone** in an equation (`nk` in `- nk lap[Psi]`) is *running*: its η comes from that equation.
   - A reference that appears **only times an amplitude** (`Dphik aDphik`) is *fixed*: you fix its η at run time (e.g. a forcing
     whose energy injection must stay finite); the equation becomes the flow of the amplitude, and the file lists the η in
     `MUST_FIX`. `"References" -> {Dphik -> {"phi", "Mode" -> "Running"}}` instead takes η from the amplitude equation with the
     amplitude counted as 1 (the convention of the old hand-built exports; add `"Equation" -> nk ank` when several qualify).
   - `NonRenormalised` couplings keep their canonical part only.
   - `ShowBeta[adim, gk]` and `ShowEtaSystem[adim]` print them.
   - **Causal poles.** The propagator's validity conditions (all roots of Δ(ω) in the upper half plane) are rewritten in
     minimal dimensionless form, in `adim["CausalPoles"]`: denominators cleared, positive factors (numbers, powers of k, of the
     references, of aq) stripped, assuming the references are > 0. Shallow water gives `Ak + ank + arnh(aq²) + arnphi(aq²) > 0`
     and `agk aHk alk + aq² (ag2k aHk alk + (Ak + arnh)(ank + arnphi)) > 0`; KPZ gives `1 + rn(aq²) > 0` (shape positivity is
     not used: the shapes are chosen in Julia). Each must hold at every aq. `ShowCausalPoles[adim]` prints them
     ([design/CAUSALITY.md](../../design/CAUSALITY.md)).
3. **Export.** `ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{dir, "eqs", "Equations2D.jl"}]]` writes one file per
   space dimension, following FlowRG's `docs/equations-interface.md` (COUPLINGS, ETAS, SHAPES, `betal`, `betanl`, `etaMl`, …). Load
   it with `load_model(path; regulator = Litim, fixed_etas = (phi = 0.0,))`, giving every `MUST_FIX` eta.
   - The file also defines `causal_poles(couplings…, aq)`, the causal-pole conditions: each entry must be > 0 at every aq (only
     the sign means something). FlowRG does not call it: evaluate it yourself, e.g. `M.causal_poles(u..., aq)` on a grid of aq,
     with the shapes bound as in FlowRG. Omitted when no condition is left.

## From a saved stage

From a stage saved by part 1 (e.g. in another notebook):

```mathematica
equations = LoadStage[FileNameJoin[{Dir, "swLPAplusequations.wdx"}]];
dims = DeriveDimensions[AnsatzOf[equations], "References" -> {nk -> "nu", Dhk -> "h", Dphik -> "phi"}, "Dimensionless" -> {Ak, aDphik}];
adim = AdimensionalEquations[equations, dims, "State" -> {gk, lk, Hk, Ak, mk, Dmk, aDphik}];
ExportJulia[adim, "Dimension" -> 2, "File" -> FileNameJoin[{Dir, "eqs", "Equations2D.jl"}]]
```

## Limits

- Regulator prefactors may contain reference couplings only (`nk q^2 arnh`, not `nk (1 + Ak) q^2 arnh`): ∂sR would otherwise
  contain a β itself. Lifting this is planned ([status](../STATUS.md)).
- Flows of momentum-derivative couplings (DE4plus g2, D2φ) need r″, written `D2armh(aq^2)`. FlowRG binds it for every
  shipped regulator except Litim ([docs/handoff/FlowRG-MUST_FIX.md](../handoff/FlowRG-MUST_FIX.md), which also specifies
  the `MUST_FIX` guard): use a smooth regulator for these files.
- Large files (DE4plus, ~1 MB) take minutes for Julia to compile on first use.

Checks: KPZ reproduces Canet's Eq. (20) exactly after adimensionalisation (`tests/acceptance/KPZAdimensional.wl`); DE4plus agrees
with the old hand-built `Equations2D.jl` up to its two known errors (`tests/acceptance/ShallowWaterDE4Adimensional.wl`). Design
and decisions: [design/ADIMENSIONAL.md](../../design/ADIMENSIONAL.md).
