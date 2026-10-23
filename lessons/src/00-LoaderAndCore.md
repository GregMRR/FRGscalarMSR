# Lesson 00 - The loader and Core

This is the first of a series of lessons, one per module, in the order the library loads them. Each lesson explains how its module works inside: which Mathematica tools it uses, how, and why the code is organised that way.

Three kinds of boxes recur. A blue **Mathematica** box explains a command or practice the first time it appears; later lessons point back to it, and `Glossary.nb` lists them all. A grey **Why it's built this way** box explains an architectural choice and names the alternative that was rejected. **Exercises** close each lesson; their solutions are folded (click the cell bracket of "Solution" to open them).

Core does no physics. It is the machinery every stage relies on: failing checks, stages and their fingerprints, exact zero tests, the user-symbol checks, and the scalar product `scal`. Every later lesson uses it, which is why it comes first.

## 1. Loading the library

```wl
Get["~/Physics/Sources/FRGscalarMSR/FRGscalarMSR.wl"]
```

:::mma Get | `Get[file]` (or `<< file`) reads a file and evaluates every expression in it, in order.
A `.wl` file is just a sequence of Mathematica expressions. `Get` parses them one by one and evaluates each, exactly as if you had typed them into cells. Definitions made in the file stay defined afterwards. `Get` returns the value of the last expression (here `Null`, because the last line is a `Print`).
:::

The lessons also need two small helpers of their own, `lessonUpTo` and `showFailure`. They live in `lessons/LessonSetup.wl`; we will read them when we use them.

```wl
Get[FileNameJoin[{$FRGRoot, "lessons", "LessonSetup.wl"}]]
```

Here is the loader in full:

@file FRGscalarMSR.wl

The first line, `(* ::Package:: *)`, is a comment: the kernel ignores it, but the front end reads it and opens the file in its package editor rather than as a notebook. Everything between `(*` and `*)` is a comment.

:::mma Symbols and contexts | Every symbol has a full name `` context`name ``; short names are looked up in `$Context`, then in `$ContextPath`.
Each name you type, like `DefineTheory` or `h`, is a **symbol**. Its full name has a context prefix ending in a backquote: `` FRGscalarMSR`DefineTheory ``, `` Global`h ``. When you type a short name, Mathematica looks for it first in the current context `$Context` (`` Global` `` in a notebook), then in each context of the list `$ContextPath`. If it finds nothing, it creates a new symbol in `$Context`. Contexts are Mathematica's namespaces: two packages can both have a symbol called `rank` without clashing, as long as they live in different contexts.
```wl
{$Context, $ContextPath, Context[DefineTheory]}
```
:::

In these lesson notebooks `$Context` is not `` Global` `` but a context private to each notebook (the builder sets the notebook option `CellContext -> Notebook`), so the variables of one lesson do not leak into another. Run headless, as `examples/runNotebook.wl` does, they live in `` Global` ``.

:::mma Set and SetDelayed | `x = rhs` evaluates rhs now and stores the result; `f[x_] := rhs` stores rhs unevaluated and evaluates it at every use.
`` FRGscalarMSR`$FRGRoot = DirectoryName[$InputFileName] `` is a `Set`: the folder of the file is computed once, at load time, and stored. A `SetDelayed` (`:=`) stores the right-hand side as it is written; it is what you use to define functions, whose body must be evaluated afresh for each argument. Note that the loader writes the full name `` FRGscalarMSR`$FRGRoot ``: at that point no package context has been opened yet, so a short name would create `` Global`$FRGRoot ``.
:::

`$InputFileName` is the path of the file currently being read by `Get`. Using it, the library finds its own `Kernel/` folder wherever it is installed.

The reload guard at the top runs only if the library was loaded before. `Contexts[]` lists every context that exists, and `MemberQ[list, x]` tests whether `x` is in the list.

:::mma Attributes and Protected | Attributes are flags on a symbol that change how it evaluates; `Protected` forbids new definitions.
`SetAttributes[f, attr]` attaches a flag. `Protected` makes a symbol read-only: any `f[...] := ...` fails with a message. The library protects the symbols the user writes in input (`dt`, `grad`, `FixedDimension`, ...) so that nobody accidentally gives them a meaning. `Unprotect` removes the flag. `` ClearAll["FRGscalarMSR`*"] `` removes the values, definitions, attributes and options of every symbol whose full name matches the string pattern.
:::

:::why Reloading starts from a clean slate
Definitions in Mathematica accumulate. If a function was defined as `f[x_Integer] := ...` and the new version of the file says `f[x_] := ...`, after a plain reload both rules exist, and the old one still wins for integers. Clearing every library symbol before reading the files again guarantees that what runs is exactly what is in the files. The `Unprotect` is needed because `ClearAll` cannot touch protected symbols.
:::

`$FRGModules` lists the module files. `Scan` reads them in that order.

:::mma Pure functions | `body &` is an anonymous function; `#` is its argument (`#1`, `#2`, ... for several).
`` Get[FileNameJoin[{FRGscalarMSR`$FRGRoot, "Kernel", # <> ".wl"}]] & `` is a function of one argument `#`: given `"Core"`, it reads `.../Kernel/Core.wl`. `<>` joins strings (`StringJoin`). `Function[x, body]` is the same thing with a named argument.
```wl
(# <> ".wl" &) /@ {"Core", "Theory"}
```
:::

:::mma Map and Scan | `Map[f, list]` (`f /@ list`) returns the list of results; `Scan[f, list]` applies f only for its side effects and returns nothing.
The loader uses `Scan` because it wants the effect of each `Get` (the definitions), not the return values. `/@` is the short form of `Map`.
:::

:::why One file per stage, loaded in a fixed order
Each module may use only the modules listed above it in `$FRGModules`. The dependencies therefore form a line, not a web: Core uses nothing, Theory uses Core, Ansatz uses Theory and Core, and so on. This has three consequences. There are no circular dependencies. You can understand the modules one at a time in load order, which is the order of these lessons. And a module can be tested on its own with the stages built by the modules before it (`tests/unit/`). The rejected alternative is a single large file, as in the older `DE4plus_fast.wl`, where any helper can depend on any other (`design/FRGDESIGN.md`, decision Q8).
:::

## 2. Packages: public, internal and private names

Every module starts and ends the same way. This is the top of Core:

@region Kernel/Core.wl | BeginPackage["FRGscalarMSR`Internal`"]; | Begin["`Private`"];

:::mma BeginPackage, Begin, End, EndPackage | The four commands that give a package its public and private namespaces.
`` BeginPackage["FRGscalarMSR`", {"FRGscalarMSR`Internal`"}] `` makes `` FRGscalarMSR` `` the current context and sets `$ContextPath` to just that context, the listed ones, and `` System` ``. Every symbol *mentioned* from now on is created in `` FRGscalarMSR` ``: these are the public names. `` Begin["`Private`"] `` switches to the sub-context `` FRGscalarMSR`Private` `` (the leading backquote means "relative to the current context"); symbols first mentioned after it, such as helper functions, are created there. `End[]` goes back, and `EndPackage[]` restores the previous `$Context` and puts `` FRGscalarMSR` `` at the front of `$ContextPath`, so its public symbols are visible to you by their short names.
:::

:::mma Usage messages | `f::usage = "text"` attaches help text to f; `?f` shows it.
In a package, the usage line has a second job: it is the first mention of the symbol, so it creates the symbol in the public context, before `` Begin["`Private`"] ``. The list of usage messages at the top of each module *is* its public interface.
```wl
?FRGAssert
```
:::

Core actually opens two packages. The first, `` FRGscalarMSR`Internal` ``, holds helpers that the modules share with each other but that are not meant for you (`assertCheck`, `stampStage`, `validateStage`, ...). The second, `` FRGscalarMSR` ``, holds the public names. Later modules list `` FRGscalarMSR`Internal` `` among the contexts they need, so they can use the internal helpers by their short names.

```wl
{Context[FRGAssert], Context[assertCheck], MemberQ[$ContextPath, "FRGscalarMSR`Internal`"]}
```

Note the last answer: `True`. The comment in `Kernel/Core.wl` says the internal helpers are "not on the user's $ContextPath", but they are, because the first `EndPackage[]` puts `` FRGscalarMSR`Internal` `` on `$ContextPath` like any package. In practice this means you can call `assertCheck` or `stampStage` from a notebook (we will, in section 4), and that a symbol of yours called `inputForm` would refer to the library's helper.

The private parts of all modules share one context. Every module calls `` Begin["`Private`"] `` inside `` BeginPackage["FRGscalarMSR`"] ``, so all helpers land in `` FRGscalarMSR`Private` ``:

```wl
Length[Names["FRGscalarMSR`Private`*"]]
```

`Names["pattern"]` returns the names of all symbols matching a string pattern.

:::why Three layers of names
The user's symbols (`h`, `nuk`, `Dim`) live in `` Global` ``, the API in `` FRGscalarMSR` ``, the shared plumbing in `` FRGscalarMSR`Internal` ``, and the implementation in `` FRGscalarMSR`Private` ``. A helper called `rank` never clashes with a coupling you call `rank`. The library made one compromise: a single private context for all modules instead of one per module (`` FRGscalarMSR`Ansatz`Private` ``, ...). It is simpler, but a helper defined in one module is visible to, and can be overwritten by, every other. That is why CLAUDE.md asks for distinctive helper names (the Regulator's `sq` is the example given there).
:::

## 3. Failing loudly: FRGAssert

Every check in the library goes through one function:

@quote Kernel/Core.wl FRGAssert failBlock

:::mma Patterns | `x_` matches any expression and names it x; `x_h` only expressions with head h; `x_ : d` is optional with default d; `patt /; cond` adds a condition.
The left-hand side of a definition is a **pattern**. `FRGAssert[name_, cond_, details_, hint_ : ""]` matches a call with three or four arguments; inside the body, `name`, `cond`, ... stand for whatever was passed. `hint_ : ""` makes the fourth argument optional. The **head** of an expression is what stands before the brackets: the head of `f[1, 2]` is `f`, of `a + b` is `Plus`, of `{1, 2}` is `List`. So `e_Plus` matches only sums. A definition is really a rule "rewrite anything that matches this pattern into this", which Mathematica stores and applies automatically.
```wl
f[x_Integer] := "integer"; f[x_] := "something else"; {f[3], f[3.5], f[a + b]}
```
:::

:::mma Hold attributes | `HoldAll`, `HoldFirst`, `HoldRest` stop Mathematica from evaluating some arguments before the function sees them.
Normally Mathematica evaluates every argument first, then calls the function. `SetAttributes[FRGAssert, HoldRest]` keeps all arguments but the first unevaluated. Inside the body, `cond` is evaluated where it is used, inside `TrueQ`. `details` and `hint` are evaluated only when `failBlock` is called, that is, only when the check fails.
```wl
g[x_] := "done"; SetAttributes[g2, HoldAll]; g2[x_] := "done";
{g[Print["evaluated"]], g2[Print["not evaluated"]]}
```
:::

The library uses `HoldRest` for speed: `details` is often something like `Select[monomials, ...]` over a large list, and on success it should cost nothing.

:::mma If and TrueQ | `If[c, a, b]` gives a if c is True, b if c is False, and stays unevaluated otherwise; `TrueQ[c]` is True only if c is literally True.
In Mathematica a condition can be neither true nor false: `nuk > 0` stays as it is when `nuk` has no value, and `If[nuk > 0, a, b]` then returns itself unevaluated. `TrueQ` turns "undecided" into `False`. For a check this is the safe choice: a condition the library cannot decide is a failed check, not a passed one.
```wl
{If[nuk > 0, "yes", "no"], TrueQ[nuk > 0]}
```
:::

A check that holds returns its name:

```wl
FRGAssert["Demo.positive", 3 > 0, "never computed"]
```

A failing check prints its FAIL block and aborts. In this notebook we wrap it in `showFailure`, which shows the block and returns the name of the check instead of stopping:

```wl
showFailure[FRGAssert["Demo.positive", -3 > 0, {-3}, "the number must be positive"]]
```

:::mma Abort, CheckAbort, Exit | `Abort[]` stops the whole evaluation; `CheckAbort[expr, alt]` catches an abort inside expr and returns alt; `Exit[n]` quits the kernel with exit code n.
`Abort[]` unwinds every function call back to the top level; the cell prints `$Aborted`. Nothing downstream runs. `CheckAbort` is the only way to catch it. `Exit[2]` ends the kernel process itself; the library uses it when it runs as a script (`math -script`), so the shell sees a non-zero exit code. `$FRGExitOnFail` chooses between the two: `$FrontEnd === Null` is True when no front end is attached, and `$CommandLine` holds the command-line arguments.
:::

Here is `showFailure`:

@quote lessons/LessonSetup.wl showFailure

:::mma Block and Module | `Module[{x}, body]` creates a fresh local x seen only by the code written inside; `Block[{x = v}, body]` temporarily gives the existing x the value v, seen by every function called from body.
`Module` is **lexical** scoping: `Module[{x}, ...]` renames `x` to a unique `x$123`, so only the text inside the `Module` can refer to it. `Block` is **dynamic** scoping: during the evaluation of the body, the global symbol itself has the new value, for every function that runs, however deep, and the old value comes back afterwards. `showFailure` uses `Block` on `$FRGExitOnFail` because `FRGAssert`, called deep inside some constructor, must see the changed value. Every constructor does the same with `Block[{$FRGModule = "Theory"}, ...]`, which is how `failBlock` knows which module to print.
```wl
yy = "global"; show[] := yy;
{Module[{yy = "module"}, show[]], Block[{yy = "block"}, show[]], yy}
```
:::

`failBlock` prints plain lines in a script and a red framed column in a notebook. Constructors do not call `FRGAssert` directly for the checks they report; they call `assertCheck`, which asserts and then returns a small record for the stage's card:

@quote Kernel/Core.wl assertCheck infoCheck

:::why Named checks that abort
The library is meant to be read as a set of grey boxes: you read the cards and the FAIL blocks, not the code. A failure therefore has a name (`Ansatz.MSR`, `Stage.fingerprint`) that you can look up in `docs/MODULES.md`, the offending part of the input, and a hint. It then aborts. Mathematica's own convention is different: a function that fails prints a message and returns unevaluated or `$Failed`, and evaluation continues. In a pipeline whose expressions reach a million leaves, that convention lets a wrong intermediate result travel silently into the flow equations. Aborting stops at the first wrong step. `Exit[2]` makes the same failure visible to shell scripts and to the test runner.
:::

## 4. Stages: immutable, fingerprinted Associations

Every module returns a **stage**: an `Association` holding everything it computed, plus bookkeeping.

:::mma Association | `<|key -> value, ...|>` is a table of key-value pairs; read with `a[key]` or `Lookup[a, key, default]`.
Keys are usually strings. `a["Stage"]` reads a key, `Lookup[a, "Warnings", {}]` reads it with a default when it is absent, `Keys[a]` lists the keys. `Append[a, k -> v]` and `KeyDrop[a, {k1, ...}]` return a **new** association and leave `a` unchanged: an association is a value, not a reference. Inside a pure function, `#Status` is short for `#["Status"]`. `AssociationThread[{k1, k2} -> {v1, v2}]` builds one from two lists.
```wl
toyA = <|"Stage" -> "Toy", "Data" -> {1, 2}|>;
{toyA["Data"], Lookup[toyA, "Missing", "default"], Keys[Append[toyA, "New" -> 3]], Keys[toyA]}
```
:::

The bookkeeping code is short:

@quote Kernel/Core.wl volatileKeys contentHash stampStage validateStage

:::mma Hash | `Hash[expr, "SHA256", "HexString"]` gives a 64-character fingerprint of expr; any change to expr changes it.
A hash is a short digest of an arbitrarily large expression. Two equal expressions have the same hash, and in practice two different ones never do. Comparing hashes is a cheap way to ask "is this still exactly the expression that was built?".
:::

Let us build a toy stage. `stampStage` is an internal helper, but as we saw it is reachable from the notebook.

```wl
toy = stampStage[<|"Stage" -> "Toy", "Data" -> {a, b}, "Checks" -> {}|>];
validateStage[toy, "Toy"]["Fingerprint"]
```

A stage passed to the wrong constructor, or modified by hand, is refused:

```wl
showFailure[validateStage[toy, "Theory"]]
```

```wl
showFailure[validateStage[Append[toy, "Data" -> {a}], "Toy"]]
```

The subtle case is a symbol of the stage that gets a value after the stage was built. `contentHash` is private, so we call it by its full name:

```wl
FRGscalarMSR`Private`contentHash[toy] === toy["Fingerprint"]
```

```wl
a = 5;
{toy, toy["Data"], FRGscalarMSR`Private`contentHash[toy] === toy["Fingerprint"]}
```

The association still displays `{a, b}`: it does not re-evaluate its contents when you use it. The part you extract, `toy["Data"]`, is evaluated and gives `{5, b}`. `Normal` turns the association back into a list of rules, which does get re-evaluated; that is why `contentHash` hashes `Normal[...]`, and why the fingerprint no longer matches:

```wl
showFailure[validateStage[toy, "Toy"]]
```

```wl
Clear[a]
```

`PrintCard` prints the summary card that every constructor shows:

@quote Kernel/Core.wl PrintCard

```wl
PrintCard[Append[toy, "Checks" -> {assertCheck["Demo.check", True, Null, "", "a passing check"]}]]
```

`Count[list, KeyValuePattern["Status" -> "PASS"]]` counts the records whose `"Status"` key is `"PASS"`; `StringPadRight` pads a string to a fixed width.

:::why Immutable stages instead of global state
A classical Mathematica notebook keeps its state in global variables that later cells read and overwrite. Re-evaluating cells out of order then gives results that depend on history. Here each constructor takes stages as arguments and returns a new one; nothing is stored anywhere else (`design/FRGDESIGN.md`, Q8: "immutable stage Associations passed explicitly, no global state"). Re-running a cell is always safe, two models can live in one kernel, and any stage can be saved to disk and loaded back (lesson 12). The fingerprint closes the remaining hole: a stage that was edited, or whose couplings acquired values, would silently give wrong results downstream, so every constructor recomputes the hash of its inputs and aborts with `Stage.fingerprint`. The keys in `volatileKeys` are excluded because they describe how the stage was built (timing, card lines), not what it contains.
:::

## 5. Exact zero tests

The library constantly has to decide whether a large expression is exactly zero: whether two forms of a vertex agree, whether a projection cancels. `ExactlyZeroQ` does it numerically, but exactly:

@quote Kernel/Core.wl RandomRationalPoint systemSymbolQ userSymbols valueAt ExactlyZeroQ

Several new tools appear in these few lines.

:::mma Options | `Options[f] = {"Name" -> default}` declares named options; `f[x_, OptionsPattern[]]` accepts them; `OptionValue["Name"]` reads them.
Options are optional arguments passed by name, `ExactlyZeroQ[e, "Seeds" -> {5, 6}]`, with defaults stored in `Options[f]`. They keep a function's signature readable when it has many settings.
:::

:::mma Cases and level specifications | `Cases[e, patt, {0, Infinity}]` returns every subexpression of e, at any depth, that matches patt.
`Cases` walks the expression tree. The level specification `{0, Infinity}` means "the whole expression and everything inside it". `Heads -> False` skips heads, so in `r[y]` it sees `y` but not `r`. `userSymbols` thus collects the plain symbols (couplings, `q`, `k`, `Dim`), and `Union` sorts them and removes duplicates. `FreeQ[e, patt]`, used everywhere later, is the opposite question: True if nothing in e matches.
```wl
Cases[nuk q^2 r[q^2/k^2] + Dim, s_Symbol, {0, Infinity}, Heads -> False] // Union
```
:::

:::mma Rules and replacement | `lhs -> rhs` and `lhs :> rhs` are rules; `e /. rules` applies them to every part of e.
`e /. x -> 2` replaces `x` by 2 everywhere in `e`. With `->` the right-hand side is evaluated once, when the rule is made; with `:>` (`RuleDelayed`) it is evaluated afresh at each match, using the matched parts. `Thread[{a, b} -> {1, 2}]` makes `{a -> 1, b -> 2}`. `/.` scans the expression from the outside in and does not look inside a part it has just replaced.
```wl
{(x + x^2) /. x -> 2, {r[1], r[2]} /. r[n_] :> n^2, Thread[{a, b} -> {1, 2}]}
```
:::

:::mma Alternatives and named patterns | `p1 | p2` matches either pattern; `x : p` matches p and names the whole match x.
`atomRule = x : (Derivative[__][_][___] | (_[___])[___] | ...) :> memo[x]` matches any of three shapes of "function atom" (`r'[y]`, `f[a][y]`, `r[y]`) and calls the whole match `x`. `__` matches one or more arguments, `___` zero or more.
:::

:::mma Memoization | `f[x_] := f[x] = value` computes value the first time and stores it as a specific definition for that x.
The first call `memo[r[y]]` evaluates the right-hand side, which is itself an assignment `memo[r[y]] = 37/61`: a new, specific definition. Every later `memo[r[y]]` uses it. Here it guarantees that the same atom gets the same random number wherever it appears.
:::

Why the atoms are replaced *whole* is the subtle point. `/.` goes from the outside in, so `r'[y]` (which is `Derivative[1][r][y]`) is matched as one atom before `r` alone could be. Had the symbol `r` been replaced by a number first, `Derivative[1][37/61]` would be zero, and every derivative term would vanish without notice (`design/MOMENTUMPROJECTION.md`, step 0).

:::mma Check and Quiet | `Check[expr, alt]` returns alt if evaluating expr produced a message; `Quiet[expr, {msgs}]` evaluates expr without printing those messages.
A random point can hit a pole (`1/0`). `valueAt` then returns `$Failed` without printing anything, and `ExactlyZeroQ` retries with another seed.
:::

:::mma Exact numbers | Integers and rationals are exact (`Precision` is `Infinity`); `Together` puts a rational expression over a common denominator.
`3/7` stays a rational number forever; `0.4286` is a machine number with rounding errors. After all symbols are replaced by rationals, `Together` reduces the expression to a single exact rational, and `PossibleZeroQ` of an exact number is an exact test. `ExactlyZeroQ` asserts that what is left is an exact number: if a symbol or pattern survived the substitution, the test would be meaningless.
:::

```wl
{ExactlyZeroQ[(a + b)^2 - a^2 - 2 a b - b^2], ExactlyZeroQ[D[r[y]^2, y] - 2 r[y] r'[y]], ExactlyZeroQ[r'[y] - r''[y]]}
```

:::why Random rationals instead of Simplify
`Simplify[expr] == 0` is the obvious test, and the library forbids it (CLAUDE.md, Conventions). On expressions with 10^5 to 10^6 leaves, `Simplify` can run for hours, and when it returns something non-zero that proves nothing: it may simply have failed to find the simplification. Evaluating at random rational points is fast, has no rounding error, and is reliable: a non-zero rational function of degree D vanishes at a random point drawn from M values with probability at most D/M (the Schwartz-Zippel bound), so three independent points make an accidental pass negligible. Giving each function atom its own value tests the identity for arbitrary shape functions, not for one particular regulator.
:::

## 6. User symbols

Fields and couplings are symbols that you choose. Before using them, the library checks that they are usable:

@quote Kernel/Core.wl reservedNames assertFreshSymbols

:::mma OwnValues and DownValues | A symbol's definitions are stored as lists of rules: `OwnValues[x]` for `x = ...`, `DownValues[f]` for `f[...] := ...`, `UpValues[f]` for definitions attached to f as an argument.
A symbol with no definitions has all three empty. A field that has a value, say `h = 3`, would be replaced by 3 inside the library, so `assertFreshSymbols` requires them empty.
```wl
v = 3; w[x_] := x^2; {OwnValues[v], DownValues[w], OwnValues[nuk]}
```
:::

`Select[list, test]` keeps the elements that pass `test`; `AllTrue`, `NoneTrue` and `DuplicateFreeQ` test a whole list; `SymbolName[s]` is the short name of a symbol, as a string.

:::why Reserved names
The library writes some symbols itself into its results: the RG scale `k`, the loop momentum `q`, the external momentum `p`, the shape argument `y = q^2/k^2`. If one of your couplings were called `k`, it would become indistinguishable from the RG scale in the flow equations. Checking the names once, at the entrance, is cheaper than debugging the confusion at the exit.
:::

## 7. scal and wedge: algebra by rewriting rules

Vectors (momenta) are bare symbols, and their scalar product is `scal[a, b]`. Its algebra is entirely defined by rewriting rules:

@quote Kernel/Core.wl scal

:::mma Orderless | The attribute `Orderless` sorts a function's arguments into a canonical order before anything else.
With `Orderless`, `scal[b, a]` is turned into `scal[a, b]`: the symmetry of the scalar product is built into how the expression is stored. `Plus` and `Times` are Orderless too, which is why `b + a` prints as `a + b`.
:::

:::mma Evaluation as rewriting | Mathematica applies matching definitions repeatedly until no rule matches any more.
`scal[a + 2 b, c]` matches `scal[(a_) + (b_), c_]` and becomes `scal[a, c] + scal[2 b, c]`; the second term matches the number rule and becomes `2 scal[b, c]`; then nothing matches, and that is the result. You never call these rules: writing the expression is enough.
```wl
{scal[a + 2 b, c], scal[b, a], scal[nuk a, b], scal[-a, b]}
```
:::

The third result, `scal[a nuk, b]`, stays as it is. `scal` is linear over explicit numbers only. Vectors and couplings are both bare symbols, so `scal` cannot tell whether `nuk` in `scal[nuk a, b]` is a coupling or a vector; only numbers can be pulled out safely. The modules that build `scal` terms keep couplings outside.

`wedge[a, b]` is the antisymmetric product of d = 2 (it appears with `jac`, the Levi-Civita operator of `design/LEVICIVITA.md`). It cannot be `Orderless`; instead a conditional rule puts its arguments in canonical order and flips the sign:

@quote Kernel/Core.wl wedge wedgePairRules reduceWedge2D

:::mma FixedPoint | `FixedPoint[f, x]` applies f repeatedly, starting from x, until the result stops changing.
`reduceWedge2D` replaces products of two wedges by scalar products, expands, and repeats until no product of wedges is left.
```wl
{wedge[b, a], wedge[a, a], reduceWedge2D[wedge[a, b]^2]}
```
:::

:::why Let the evaluator do the algebra
Because `scal` and `wedge` rewrite themselves into canonical forms, two equal expressions built along different paths end up looking the same. Terms combine automatically, `Collect` and `Together` see common factors, and `ExactlyZeroQ` works on smaller expressions. The alternative, a vector package with explicit components, would fix d and multiply the size of every expression by the number of components; here d stays the symbol `Dim` to the end.
:::

## Exercises

:::exercise Read a FAIL block
Write a check named `"Me.even"` that `7` is an even number (`EvenQ`), with the offending value and a hint, and run it through `showFailure`. Then evaluate `$FRGLastFailure` and `$FRGModule`. Why does `$FRGModule` print `FRGscalarMSR` here, not a module name?
:::

:::solution
```wl
showFailure[FRGAssert["Me.even", EvenQ[7], 7, "pass an even number"]]
```
```wl
{$FRGLastFailure, $FRGModule}
```
No constructor is running, so no `Block[{$FRGModule = ...}]` is active and the global default set in `Kernel/Core.wl` is printed.
:::

:::exercise HoldRest at work
`FRGAssert` holds its `details` argument. Verify it: call it once with a true condition and once with a false one, passing `Print["details computed"]` as details. How many times is the message printed?
:::

:::solution
```wl
FRGAssert["Me.hold", True, Print["details computed"]]
```
```wl
showFailure[FRGAssert["Me.hold", False, Print["details computed"]]]
```
Once: only the failing call evaluates `details`, inside `failBlock` (the printed "offending" is then `Null`, the value of `Print`).
:::

:::exercise An identity with a shape function
Use `ExactlyZeroQ` to check that the derivative of `y r[y]` is `r[y] + y r'[y]`, and that it is not `r[y] + r'[y]`. Then ask yourself what would go wrong if `valueAt` replaced the symbol `r` by a number before the atoms.
:::

:::solution
```wl
{ExactlyZeroQ[D[y r[y], y] - (r[y] + y r'[y])], ExactlyZeroQ[D[y r[y], y] - (r[y] + r'[y])]}
```
`r'[y]` is `Derivative[1][r][y]`. With `r` replaced by 37/61 it would become `Derivative[1][37/61][y]`, the derivative of a constant, which is 0. The difference in the wrong identity is `(y - 1) r'[y]`, so it would pass.
:::
