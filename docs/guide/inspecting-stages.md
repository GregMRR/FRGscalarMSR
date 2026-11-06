# Inspecting a loaded stage

A stage is an Association: `Keys[stage]` lists its contents, `stage["Stage"]` its name, `PrintCard[stage]` its card. A stage embeds
the stages it was built from, so one loaded `Equations` file gives access to the whole chain:

```
Equations ⊃ Closure ⊃ Integrated ⊃ Projection ⊃ Loops ⊃ { Vertices, Propagator, Regulator } ⊃ Ansatz ⊃ Theory     (Symmetries beside Closure)
Adimensional ⊃ { Equations, Dimensions }
```

`EmbeddedStage[stage, "Name"]` returns the embedded stage (or `Missing`); then use that module's accessors from the
[API reference](../reference/README.md).

| I want | Call (`eqs` = a loaded Equations stage) |
|---|---|
| the ansatz, its monomials | `ShowMonomials[AnsatzOf[eqs]]` |
| vertices, Γ2 | `v = EmbeddedStage[eqs, "Vertices"]`, `ShowVertex[v, {…}]`, `ShowGamma2[v]` |
| regulator | `ShowRegulator[EmbeddedStage[eqs, "Regulator"]]` |
| propagator, Δ | `prop = EmbeddedStage[eqs, "Propagator"]`, see below |
| symmetry classes | `ShowClassification[EmbeddedStage[eqs, "Symmetries"]]` |
| closure report | `ShowClosureReport[EmbeddedStage[eqs, "Closure"]]` |
| the flow of a coupling | `ShowEquation[eqs, comb]`, `EquationOf[eqs, comb]` |
| dimensions, β function (from an Adimensional stage) | `ShowDimensions[adim["Dimensions"]]`, `ShowBeta[adim, c]` |
| causal-pole conditions, dimensionless | `ShowCausalPoles[adim]` |

## Δ of the propagator matrix

Δ(ω) = det of the (φ, φbar) block of Γ2 + R, whose roots are the poles of the retarded propagator. It is stored in the
Propagator stage as a function of the frequency and q², written in `GammaAtom` atoms:

```mathematica
eqs  = LoadStage["kpz.wdx", "Expect" -> {theory, ansatz, regulator, symmetries}];
prop = EmbeddedStage[eqs, "Propagator"];

prop["Delta"][w, s]                  (* I w + GammaAtom[h, hbar, 0][s]   for KPZ *)
prop["DeltaCoefficients"][s]         (* {c0, c1, c2}: Δ = c0 + c1 (I w) + c2 (I w)^2 *)
prop["DeltaDegree"]                  (* 1 or 2 *)
ValidityConditions[prop][s]          (* conditions for all roots in the upper half plane *)

(* with the couplings and regulator shapes written out *)
prop["Delta"][w, s] /. (g : GammaAtom[__])[x_] :> prop["AtomValues"][g][x]
                                     (* nuk s + I w + nuk s rn[s/k^2]   for KPZ *)
ShowPropagator[prop]                 (* prints atoms, Δ, coefficients, validity *)
```

## The propagator in the couplings

`prop["G"][w, s]` is written in atoms and `Δ`; two substitutions make it explicit:

```mathematica
atoms = (g : GammaAtom[__])[x_] :> prop["AtomValues"][g][x];
G = prop["G"][w, s] /. DeltaRules[prop] /. atoms;      (* matrix in field order (h, hbar), no GammaAtom or Δ left *)
G[[1, 2]]                                              (* G_{h hbar}: 1/(nuk s - I w + nuk s rn[s/k^2]) for KPZ *)
G . (prop["M"][w, s] /. atoms) // Simplify             (* identity matrix: check *)
```

The same stage holds `prop["G"][w, s]` (propagator matrix, same as `PropagatorMatrix[prop, w, s]`) and `prop["M"][w, s]`
(Γ2 + R, `FullMatrix`). From a Loops stage, `loops["Propagator"]` gives it directly.
