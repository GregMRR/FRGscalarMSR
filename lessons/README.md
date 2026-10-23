# Lessons

One notebook per module, in load order, explaining how each module works inside: the Mathematica tools it uses (blue
"Mathematica:" boxes, each explained once, at first use), the architectural reasons (grey "Why it's built this way" boxes), and
exercises with folded solutions. Every lesson is self-contained: it loads the library and rebuilds the KPZ running example
up to the stage it needs (`lessonUpTo`, in `LessonSetup.wl`). `Glossary.nb` lists every Mathematica concept with the lesson
that explains it.

| # | Lesson | # | Lesson |
|---|---|---|---|
| 00 | Loader and Core | 09 | FrequencyIntegration |
| 01 | Theory | 10 | Closure |
| 02 | Ansatz | 11 | Assemble (ends with the Canet KPZ check) |
| 03 | Vertices | 12 | Oracle and Stages |
| 04 | Regulator | 13 | Dimensions |
| 05 | Propagator | 14 | Adimensional |
| 06 | Symmetries | 15 | JuliaExport |
| 07 | Loops | 16 | Tests (the harness of `tests/runTests.wl`, benchmarks) |
| 08 | Projection | 17 | Momentum functions (flows at generic p, `Related`) |

## Building and checking

- Build: `math -script lessons/makeLessons.wl` (writes `lessons/NN-Name.nb` and `lessons/Glossary.nb`).
- Check one lesson headless: `math -script examples/runNotebook.wl lessons/NN-Name.nb` (exit 0 = every cell evaluated).
- Edit the sources `lessons/src/NN-Name.md`, never the `.nb` files.

The build stops when a quoted symbol no longer exists in the code, when a glossary reference comes before the concept is
explained, or when a concept is explained twice. Code quotes are read from the source files at build time, so they follow
the code; only the prose can go stale.

## Source format (`lessons/src/NN-Name.md`)

- `#`, `##`, `###`, `####`: title, section, subsection, subsubsection.
- Paragraphs are separated by blank lines; `- ` starts a bullet. Inline: `` `code` `` (code containing a backtick: ```` `` Global` `` ````), `**bold**`, `*word*` (italic), `[[Concept]]` (a reference
  to a concept already explained; rendered with the lesson number).
- ```` ```wl ```` fences: Input cells (evaluated by the headless check). ```` ```text ````: displayed only.
- `:::mma Concept | one-line summary` … `:::`: a Mathematica box; registers the concept in the Glossary.
- `:::why Title` … `:::`: an architecture box.
  Boxes hold text and bullets; their code fences become Input cells right after the box.
- `:::exercise Title` … `:::` then `:::solution` … `:::` (folded). A solution must be followed by a heading, another
  exercise, or the end of the file.
- `@quote path sym1 sym2 …`: the definitions of the symbols (with `Options`, `SetAttributes` and the comment just above).
  `@region path | from | to`: the lines from the first one containing `from` to the next containing `to`.
  `@file path`: a whole file. Paths are relative to the repository root.
- A cell that should show a failing check wraps it in `showFailure[…]`, which prints the FAIL block and returns the check's name
  instead of aborting (or exiting, in the headless run).
