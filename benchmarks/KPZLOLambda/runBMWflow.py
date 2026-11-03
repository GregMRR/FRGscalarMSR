"""Acceptance 2 of design/LOFLAMBDA.md: KPZ at LO with and without a flowing f_lambda, integrated by BMWflow in d = 1, 2, 3.

Prints z = 2 - eta_nu and chi = (2 - d + eta_D - eta_nu)/2 at the end of each run (the fixed point, if reached), and saves
PNGs of eta_D(s), eta_nu(s) and of the functions f_X(p) at the end in benchmarks/KPZLOLambda/out/ (git-ignored). The two
runs are compared by the user. Needs BMWflow with contract v3 (docs/handoff/BMWflow_LO.md) and its .venv:

    ~/sources/BMWflow/.venv/bin/python benchmarks/KPZLOLambda/runBMWflow.py [--dims 1 2 3] [--s-fin -20] [--g0 1:1,2:50,3:300]

The modules are exported first by  math -script benchmarks/KPZLOLambda/export.wl  (writes out/KPZ_LO.py, out/KPZLambda_LO.py).
Grids and regulator: those of BMWflow's usage_examples/KPZ_1D_twogrids (tests/test_model_generated.py). In d = 1 the
f_lambda model breaks the KPZ time-reversal symmetry (design/LOFLAMBDA.md section 5): its exponents are not protected there.
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
BMW = Path(os.environ.get("BMWFLOW_PATH", Path.home() / "sources" / "BMWflow"))
sys.path.insert(0, str(BMW))

from flow import flow_dataclasses, regulator  # noqa: E402
from flow.evolution import Evolution  # noqa: E402
from flow.models import ModelGenerated  # noqa: E402
from flow.models.model_generated import load_generated  # noqa: E402

ALPHA, BETA = 2, 1
r = lambda x: regulator.wett_alpha_beta(x, ALPHA, BETA)
r_ = lambda x: regulator.wett_alpha_beta_(x, ALPHA, BETA)
GRID_EXT = flow_dataclasses.ParamsGridExternal(N_p=100, p_min=0.01, p_max=100, N_w=100, w_min=0.01, w_max=100)
GRID_INT = flow_dataclasses.ParamsGridInternal(deg_q=16, q_max=4, deg_theta=10)
MODULES = {"KPZ (f_lambda = 1)": OUT / "KPZ_LO.py", "KPZ + f_lambda": OUT / "KPZLambda_LO.py"}


def run(path, dim, g0, s_fin, ds):
    """Euler steps of BMWflow's Evolution; returns the history (s, eta, g) and the final model state.
    Evolution's own files go to a fresh temporary directory (it asks for input if its path exists)."""
    tmp = tempfile.TemporaryDirectory(prefix="kpzlolambda-")
    mod = load_generated(str(path))
    n_eta = sum(e is not None for e in mod.ETAS)
    model = ModelGenerated(str(path), "LO", dim, GRID_EXT, GRID_INT, r, r_, "const")
    ic = flow_dataclasses.IC_NLO(g_in=g0, etas_in=np.zeros(n_eta), fs_in=np.ones((len(mod.FUNCTIONS), GRID_EXT.N_p, 1)))
    evo = Evolution(model, ic, ds, str(Path(tmp.name) / "run"))
    hist, s, n = [], 0.0, 0
    while s > s_fin:
        evo.step_update()
        if n % 100 == 0:
            hist.append((s, *np.atleast_1d(evo.eta), evo.g))
        n += 1
        s -= ds
    evo.close_files()
    tmp.cleanup()
    names = [e for e in mod.ETAS if e is not None]
    return np.array(hist), names, mod.FUNCTIONS, model.p, np.array(evo.f)


def exponents(eta, names, dim):
    e = dict(zip(names, eta))
    eD, enu = e["etaD"], e["etanu"]
    return 2 - enu, (2 - dim + eD - enu) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dims", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--s-fin", type=float, default=-20.0)
    ap.add_argument("--ds", type=float, default=1e-3)
    ap.add_argument("--g0", default="1:1,2:50,3:300", help="initial g per dimension, d:g,...")
    a = ap.parse_args()
    g0 = {int(k): float(v) for k, v in (x.split(":") for x in a.g0.split(","))}
    for p in MODULES.values():
        if not p.exists():
            sys.exit(f"{p} missing: run  math -script benchmarks/KPZLOLambda/export.wl  first")
    OUT.mkdir(exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
        print("matplotlib not found: no plots, table only")

    rows = []
    for dim in a.dims:
        res = {}
        for label, path in MODULES.items():
            hist, names, funcs, p, f = run(path, dim, g0[dim], a.s_fin, a.ds)
            eta = hist[-1, 1:1 + len(names)]
            drift = np.max(np.abs(hist[-1, 1:1 + len(names)] - hist[max(0, len(hist) - 11), 1:1 + len(names)]))
            z, chi = exponents(eta, names, dim)
            rows.append((dim, label, z, chi, hist[-1, -1], drift))
            res[label] = (hist, names, funcs, p, f)
        if plt is not None:
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            for label, (hist, names, funcs, p, f) in res.items():
                for i, nm in enumerate(names):
                    ax[0].plot(hist[:, 0], hist[:, 1 + i], label=f"{nm}, {label}")
                for i, fn in enumerate(funcs):
                    ax[1].semilogx(p, f[i, :, 0] if f.ndim == 3 else f[i], label=f"{fn}, {label}")
            ax[0].set_xlabel("s"); ax[0].set_ylabel("eta"); ax[0].legend(fontsize=7)
            ax[1].set_xlabel("p"); ax[1].set_ylabel("f(p) at s_fin"); ax[1].legend(fontsize=7)
            fig.suptitle(f"KPZ LO, d = {dim}")
            fig.tight_layout(); fig.savefig(OUT / f"kpz_flambda_d{dim}.png", dpi=130); plt.close(fig)

    print(f"\n{'d':>2}  {'model':<20} {'z':>10} {'chi':>10} {'g':>10}  {'|d eta| last unit of s':>22}")
    for dim, label, z, chi, g, drift in rows:
        print(f"{dim:>2}  {label:<20} {z:>10.5f} {chi:>10.5f} {g:>10.4f}  {drift:>22.2e}")
    print(f"\nPlots: {OUT}/kpz_flambda_d<d>.png. A large |d eta| means the run has not reached its fixed point (raise --s-fin or change --g0).")


if __name__ == "__main__":
    main()
