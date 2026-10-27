(* ::Package:: *)

(* DE regression guard (design/LOPIPELINE.md step 0): the derivative-expansion stages and the exported .jl text are unchanged
   since the LO pipeline was added. Content hashes (volatile keys dropped at every level, .jl without its header) of KPZ,
   shallow-water DE4 and Kuramoto-Sivashinsky, against the values recorded at commit bfe883a (tests/reference/DEBaseline.wl). *)

Get[FileNameJoin[{testRoot, "reference", "DEBaseline.wl"}]];
deBaselineNow = deBaselineCompute[];
KeyValueMap[test["DE baseline unchanged: " <> #1, deBaselineNow[#1], #2] &, deBaselineExpected];
