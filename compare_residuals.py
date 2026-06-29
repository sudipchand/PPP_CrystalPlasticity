#!/usr/bin/env python3
"""
compare_residuals.py - overlay the AL residual from two runs (e.g. base code
vs. Phase-1 code) at two identical measurement points.

Both UMATs write these records to unit 16 (.feh):

  NWEND  l_iter_fp  ||residual_nw||  n_active     (after the inner Newton loop ends,
                                                   once per fixed-point iteration)
  FPEND  l_iter_fp  ||residual_nw||  n_active     (after the fixed-point loop ends,
                                                   once per converged increment)

||residual_nw|| is the augmented-Lagrangian residual norm (dimensionless), defined
identically in both codes.

Usage:
  python compare_residuals.py base.feh phase1.feh
  python compare_residuals.py base.feh phase1.feh --labels base phase1
  python compare_residuals.py base.feh phase1.feh --out compare.png
  python compare_residuals.py base.feh phase1.feh --tail 300   # zoom last N records
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read(path, tag):
    """Return arrays (residual, n_active) for records matching `tag`."""
    res, nact = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 4 and p[0] == tag:
                try:
                    res.append(float(p[2]))
                    nact.append(int(p[3]))
                except (ValueError, IndexError):
                    continue
    return np.array(res, dtype=float), np.array(nact, dtype=int)


def safelog(x):
    x = np.array(x, dtype=float)
    x[x <= 0] = np.nan
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh_a", help="first .feh (e.g. base code)")
    ap.add_argument("feh_b", help="second .feh (e.g. phase-1 code)")
    ap.add_argument("--labels", nargs=2, default=["run A", "run B"])
    ap.add_argument("--tail", type=int, default=0,
                    help="only plot the last N records of each series (0 = all)")
    ap.add_argument("--out", default="compare_residuals.png")
    args = ap.parse_args()

    la, lb = args.labels

    # read both measurement points from both files
    nw_a, nwact_a = read(args.feh_a, "NWEND")
    nw_b, nwact_b = read(args.feh_b, "NWEND")
    fp_a, fpact_a = read(args.feh_a, "FPEND")
    fp_b, fpact_b = read(args.feh_b, "FPEND")

    if nw_a.size == 0 and fp_a.size == 0:
        print(f"No NWEND/FPEND records in {args.feh_a}", file=sys.stderr)
        sys.exit(1)
    if nw_b.size == 0 and fp_b.size == 0:
        print(f"No NWEND/FPEND records in {args.feh_b}", file=sys.stderr)
        sys.exit(1)

    def tail(x):
        return x[-args.tail:] if args.tail > 0 and x.size else x

    nw_a, nwact_a = tail(nw_a), tail(nwact_a)
    nw_b, nwact_b = tail(nw_b), tail(nwact_b)
    fp_a, fpact_a = tail(fp_a), tail(fpact_a)
    fp_b, fpact_b = tail(fp_b), tail(fpact_b)

    # textual summary
    def summ(name, r):
        if r.size:
            print(f"  {name:14s}: n={r.size:6d}  "
                  f"min={np.nanmin(safelog(r)):.2e}  "
                  f"median={np.nanmedian(safelog(r)):.2e}  "
                  f"max={np.nanmax(safelog(r)):.2e}")
        else:
            print(f"  {name:14s}: (no records)")
    print("residual after Newton loop (NWEND):")
    summ(la, nw_a); summ(lb, nw_b)
    print("residual after fixed-point loop (FPEND):")
    summ(la, fp_a); summ(lb, fp_b)
    print("final active-system count (last FPEND record):")
    print(f"  {la}: {fpact_a[-1] if fpact_a.size else 'NA'}    "
          f"{lb}: {fpact_b[-1] if fpact_b.size else 'NA'}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 13))

    # Panel 1: residual after Newton loop ends
    ax = axes[0]
    ax.semilogy(np.arange(nw_a.size), safelog(nw_a), "-", lw=0.8, label=la)
    ax.semilogy(np.arange(nw_b.size), safelog(nw_b), "-", lw=0.8, label=lb)
    ax.set_title("Residual after inner Newton loop ends (per fixed-point iteration)")
    ax.set_xlabel("NWEND record # (in order)")
    ax.set_ylabel(r"$\|R\|_2$")
    ax.grid(True, which="both", alpha=0.3); ax.legend()

    # Panel 2: residual after fixed-point loop ends
    ax = axes[1]
    ax.semilogy(np.arange(fp_a.size), safelog(fp_a), "o-", ms=3, lw=0.8, label=la)
    ax.semilogy(np.arange(fp_b.size), safelog(fp_b), "o-", ms=3, lw=0.8, label=lb)
    ax.set_title("Residual after fixed-point loop ends (per converged increment)")
    ax.set_xlabel("FPEND record # (in order)")
    ax.set_ylabel(r"$\|R\|_2$")
    ax.grid(True, which="both", alpha=0.3); ax.legend()

    # Panel 3: active-system count at convergence (the key comparison)
    ax = axes[2]
    ax.plot(np.arange(fpact_a.size), fpact_a, "s-", ms=3, lw=0.8, label=la)
    ax.plot(np.arange(fpact_b.size), fpact_b, "^-", ms=3, lw=0.8, label=lb)
    ax.set_title("Active slip-system count at fixed-point convergence")
    ax.set_xlabel("FPEND record # (in order)")
    ax.set_ylabel("active systems")
    ax.grid(True, alpha=0.3); ax.legend()
    allact = np.concatenate([a for a in (fpact_a, fpact_b) if a.size])
    if allact.size:
        ax.set_ylim(-0.5, allact.max() + 0.5)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nsaved: {args.out}")


if __name__ == "__main__":
    main()
