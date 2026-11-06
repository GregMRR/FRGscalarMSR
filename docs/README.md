<!-- verified: c5fc51d 2026-11-06 -->
# Documentation

Start with the [project README](../README.md) (install, quickstart), then [status](STATUS.md) for what is implemented and
planned.

## Guides

- [Flow equations, step by step](guide/flow-equations.md): part 1 of the derivative-expansion pipeline, and the 2D `jac` term.
- [Dimensionless flows and the FlowRG file](guide/dimensionless-and-export.md): part 2.
- [LO pipeline](guide/lo-pipeline.md): couplings that are functions of p², the BMWflow module.
- [Inspecting a loaded stage](guide/inspecting-stages.md): what a stage holds and how to reach it.
- [Development](guide/development.md): tests, examples, benchmarks, lessons, the documentation tools.

## Reference

- [API reference](reference/README.md): every exported symbol, generated from the `::usage` strings.
- [Modules](MODULES.md): objective, assumptions, manipulations and checks of each module.
- [Conventions](CONVENTIONS.md): Fourier transform, vertex legs, the density language.
- [Status](STATUS.md): implemented features, checks, limitations, future directions.

## Design and handoffs

- [Design records](../design/): one decision record per feature, with its adversarial reviews in `design/reviews/`.
- [Handoffs](handoff/): specifications for work in other repositories (FlowRG, BMWflow).
