#!/usr/bin/env python3
"""
view_run_split.py - save THREE separate diagnostic images from one run's .feh.

The UMAT writes three record types to unit 16:

  NEWTON  l_iter_fp  l_iter_nw  ||R||  tol_nw  n_active   (every inner Newton step)
  FIXPT   l_iter_fp  ||R||      n_active                  (every fixed-point pass)
  CONVRG  l_iter_fp  ||R||      n_active                  (each converged increment)

Produces:
  <out>_newton_residual.png   - ||R|| at every Newton iteration   (NEWTON)
  <out>_fixpt_residual.png    - ||R|| at every fixed-point pass    (FIXPT)
  <out>_active_systems.png    - active slip systems per increment  (CONVRG)

Usage:
  python view_run_split.py cp_uniaxial_1.feh
  python view_run_split.py cp_uniaxial_1.feh --label phase1 --out phase1
  python view_run_split.py cp_uniaxial_1.feh --tail 500
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read(path, tag, res_col, act_col):
    res, act = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            p = line.split()
            if len(p) > max(res_col, act_col) and p[0] == tag:
                try:
                    res.append(float(p[res_col]))
                    act.append(int(p[act_col]))
                except (ValueError, IndexError):
                    continue
    return np.array(res, dtype=float), np.array(act, dtype=int)


def safelog(x):
    x = np.array(x, dtype=float)
    x[x <= 0] = np.nan
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh")
    ap.add_argument("--label", default="run")
    ap.add_argument("--tail", type=int, default=0)
    ap.add_argument("--out", default="view_run")
    args = ap.parse_args()

    nw_res, nw_act = read(args.feh, "NEWTON", 3, 5)
    fp_res, fp_act = read(args.feh, "FIXPT", 2, 3)
    cv_res, cv_act = read(args.feh, "CONVRG", 2, 3)

    if nw_res.size == 0 and fp_res.size == 0 and cv_res.size == 0:
        print(f"No NEWTON/FIXPT/CONVRG records found in {args.feh}", file=sys.stderr)
        sys.exit(1)

    def tail(x):
        return x[-args.tail:] if args.tail > 0 and x.size else x

    nw_res = tail(nw_res)
    fp_res = tail(fp_res)
    cv_res, cv_act = tail(cv_res), tail(cv_act)

    def summ(name, r):
        if r.size:
            print(f"  {name:7s}: n={r.size:7d}  min={np.nanmin(safelog(r)):.2e}"
                  f"  median={np.nanmedian(safelog(r)):.2e}  max={np.nanmax(safelog(r)):.2e}")
        else:
            print(f"  {name:7s}: (no records)")
    print(f"run: {args.label}")
    summ("NEWTON", nw_res); summ("FIXPT", fp_res); summ("CONVRG", cv_res)
    if cv_act.size:
        vals, counts = np.unique(cv_act, return_counts=True)
        print("active-system distribution (CONVRG): " +
              "  ".join(f"{v}:{c}" for v, c in zip(vals, counts)))

    if nw_res.size:
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.semilogy(np.arange(nw_res.size), safelog(nw_res), "-", lw=0.6, color="tab:blue")
        ax.set_title(f"[{args.label}] Residual at every Newton iteration")
        ax.set_xlabel("Newton iteration # (cumulative, in order)")
        ax.set_ylabel(r"$\|R\|_2$")
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        o = f"{args.out}_newton_residual.png"; fig.savefig(o, dpi=150); plt.close(fig); print(f"saved: {o}")

    if fp_res.size:
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.semilogy(np.arange(fp_res.size), safelog(fp_res), "-", lw=0.8, color="tab:orange")
        ax.set_title(f"[{args.label}] Residual at every fixed-point iteration")
        ax.set_xlabel("fixed-point iteration # (cumulative, in order)")
        ax.set_ylabel(r"$\|R\|_2$")
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        o = f"{args.out}_fixpt_residual.png"; fig.savefig(o, dpi=150); plt.close(fig); print(f"saved: {o}")

    if cv_act.size:
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(np.arange(cv_act.size), cv_act, "s-", ms=3, lw=0.8, color="tab:green")
        ax.set_title(f"[{args.label}] Active slip-system count at convergence")
        ax.set_xlabel("converged increment # (in order)")
        ax.set_ylabel("active systems")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.5, cv_act.max() + 0.5)
        fig.tight_layout()
        o = f"{args.out}_active_systems.png"; fig.savefig(o, dpi=150); plt.close(fig); print(f"saved: {o}")


if __name__ == "__main__":
    main()
