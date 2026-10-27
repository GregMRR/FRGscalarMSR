(* ::Package:: *)

(* DE regression baseline (design/LOPIPELINE.md step 0). Content hashes of the DE stages and of the exported .jl text, recorded
   before the LO pipeline was added; tests/acceptance/DEBaseline.wl recomputes them. Reproducible across kernels and histories:
   the volatile keys (Timing, Checks, Fingerprint, Warnings, Upstream, Seconds) are dropped at every nesting level (the
   stage fingerprint drops them at the top level only), and the .jl text loses its two header lines (date, fingerprint).
   Needs testRoot and quiet (the harness, or a recording script that defines them). *)

deBaselineVolatile = {"Timing", "Checks", "Fingerprint", "Warnings", "Upstream", "Seconds"};
deBaselineStrip[a_Association] := deBaselineStrip /@ KeyDrop[a, deBaselineVolatile];
deBaselineStrip[l_List] := deBaselineStrip /@ l;
deBaselineStrip[e_] := e;
deBaselineHash[stage_Association] := Hash[deBaselineStrip[stage], "SHA256", "HexString"];   (* the form of the r5/r6 probes *)
deBaselineJlHash[file_String] := Hash[StringRiffle[Drop[StringSplit[Import[file, "Text"], "\n", All], 2], "\n"], "SHA256", "HexString"];

(* every hash of the baseline, rebuilt from the shared setups (no reliance on other test files' globals beyond these setups) *)
deBaselineCompute[] := Module[{dir = FileNameJoin[{$TemporaryDirectory, "FRGscalarMSR-baseline-" <> ToString[$ProcessID]}], kpzD, kpzA, swD, swA, jl},
   Get[FileNameJoin[{testRoot, "reference", "KPZSetup.wl"}]];
   Get[FileNameJoin[{testRoot, "reference", "ShallowWaterSetup.wl"}]];
   quiet[Get[FileNameJoin[{ParentDirectory[testRoot], "benchmarks", "KuramotoSivashinsky", "Pipeline.wl"}]]];
   kpzD = quiet[DeriveDimensions[kpzAnsatz, "References" -> {nuk -> "nu", Dk -> "D"}]];
   kpzA = quiet[AdimensionalEquations[kpzEquations, kpzD]];
   swD = quiet[DeriveDimensions[swaAnsatz, "References" -> {Dphik -> {"phi", "Mode" -> "Running"}, Dhk -> {"h", "Mode" -> "Running"},
        nk -> {"nu", "Mode" -> "Running", "Equation" -> nk ank}}, "Dimensionless" -> {Ak, ank, aDhk, aDphik}]];
   swA = quiet[AdimensionalEquations[swaEquations, swD, "State" -> {gk, g2k, lk, Hk, aDphik, D2phik, Kk, Ak, ank, aDhk}]];
   jl[name_, adim_, d_] := deBaselineJlHash[quiet[ExportJulia[adim, "Dimension" -> d, "File" -> FileNameJoin[{dir, name <> ToString[d] <> "D.jl"}]]]];
   <|"KPZ Equations" -> deBaselineHash[kpzEquations], "ShallowWater DE4 Equations" -> deBaselineHash[swaEquations],
    "KS Equations" -> deBaselineHash[ksEquations],
    "KPZ Adimensional" -> deBaselineHash[kpzA], "ShallowWater DE4 Adimensional" -> deBaselineHash[swA],
    "KPZ .jl d = 2" -> jl["KPZ", kpzA, 2], "KPZ .jl d = 3" -> jl["KPZ", kpzA, 3], "ShallowWater DE4 .jl d = 2" -> jl["SW", swA, 2]|>];

(* recorded at commit bfe883a (branch BMW-LO), identical in two fresh kernels *)
deBaselineExpected = <|
   "KPZ Equations" -> "2306608ad3ff5b1deb909aa468d941d5d55562c1c92e4ed4c41904744e69bf8c",
   "ShallowWater DE4 Equations" -> "d3e8af22e17a7effbd9d1574cd36112dea2738b2ee169aac18265f6dddbc4e6d",
   "KS Equations" -> "3a053cfbad1f651720132cc78a1ffe8ab7a7fce2314ffcd7f40a63ea726f7fcc",
   "KPZ Adimensional" -> "597b779e7005c56acc6f8adcfc2b73c354fd41bb3ccdada8491da98b262d5922",
   "ShallowWater DE4 Adimensional" -> "d26d67b8ea33c61457b450f20173a601cf60a3b80e675ace6a5c941e99b1b805",
   "KPZ .jl d = 2" -> "ad530621b6c819f422fa49393da7b79d0da0e5f10beee3027d69ba1c88f6a9bb",
   "KPZ .jl d = 3" -> "b9c318222d7e9ed285891f956b8abd5caa563629e7978c2e640f5bfd1813057e",
   "ShallowWater DE4 .jl d = 2" -> "2d8ba5f28dc0d15f84a84f8a2ff2850a22bf6808a5f2db6f8c81182a19e7e3ce"|>;
