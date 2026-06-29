#!/usr/bin/env python3
"""
view_run.py - plot the AL residual from ONE run at the two measurement points.

Reads the records written by the UMAT to unit 16 (.feh):

  NWEND  l_iter_fp  ||residual_nw||  n_active   (after inner Newton loop ends,
                                                 once per fixed-point iteration)
  FPEND  l_iter_fp  ||residual_nw||  n_active   (after fixed-point loop ends,
                                                 once per converged increment)

||residual_nw|| is the augmented-Lagrangian residual norm (dimensionless).

Usage:
  python view_run.py cp_uniaxial_1.feh
  python view_run.py cp_uniaxial_1.feh --label phase1
  python view_run.py cp_uniaxial_1.feh --out run.png
  python view_run.py cp_uniaxial_1.feh --tail 300     # zoom last N records
  python view_run.py cp_uniaxial_1.feh --print        # also dump table to stdout
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read(path, tag):
    """Return (residual, n_active) arrays for records matching `tag`."""
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
    ap.add_argument("feh", help="path to the .feh file for one run")
    ap.add_argument("--label", default="run")
    ap.add_argument("--tail", type=int, default=0,
                    help="only plot the last N records (0 = all)")
    ap.add_argument("--print", dest="do_print", action="store_true",
                    help="print the FPEND table to the terminal")
    ap.add_argument("--out", default="view_run.png")
    args = ap.parse_args()

    nw, nwact = read(args.feh, "NWEND")
    fp, fpact = read(args.feh, "FPEND")

    if nw.size == 0 and fp.size == 0:
        print(f"No NWEND/FPEND records found in {args.feh}", file=sys.stderr)
        sys.exit(1)

    def tail(x):
        return x[-args.tail:] if args.tail > 0 and x.size else x

    nw, nwact = tail(nw), tail(nwact)
    fp, fpact = tail(fp), tail(fpact)

    # textual summary
    def summ(name, r):
        if r.size:
            print(f"  {name:6s}: n={r.size:6d}  "
                  f"min={np.nanmin(safelog(r)):.2e}  "
                  f"median={np.nanmedian(safelog(r)):.2e}  "
                  f"max={np.nanmax(safelog(r)):.2e}")
        else:
            print(f"  {name:6s}: (no records)")
    print(f"run: {args.label}")
    print("residual after Newton loop (NWEND):"); summ("NWEND", nw)
    print("residual after fixed-point loop (FPEND):"); summ("FPEND", fp)
    if fpact.size:
        vals, counts = np.unique(fpact, return_counts=True)
        dist = "  ".join(f"{v}:{c}" for v, c in zip(vals, counts))
        print(f"active-system count distribution at convergence: {dist}")
        print(f"final active-system count: {fpact[-1]}")

    if args.do_print and fp.size:
        print(f"\n{'#':>6}  {'||R||':>14}  {'n_active':>8}")
        for i in range(fp.size):
            print(f"{i:>6}  {fp[i]:>14.6e}  {fpact[i]:>8}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 12))

    ax = axes[0]
    ax.semilogy(np.arange(nw.size), safelog(nw), "-", lw=0.8, color="tab:blue")
    ax.set_title(f"[{args.label}] Residual after inner Newton loop "
                 f"(per fixed-point iteration)")
    ax.set_xlabel("NWEND record # (in order)")
    ax.set_ylabel(r"$\|R\|_2$")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[1]
    ax.semilogy(np.arange(fp.size), safelog(fp), "o-", ms=3, lw=0.8,
                color="tab:orange")
    ax.set_title(f"[{args.label}] Residual after fixed-point loop "
                 f"(per converged increment)")
    ax.set_xlabel("FPEND record # (in order)")
    ax.set_ylabel(r"$\|R\|_2$")
    ax.grid(True, which="both", alpha=0.3)

    ax = axes[2]
    ax.plot(np.arange(fpact.size), fpact, "s-", ms=3, lw=0.8, color="tab:green")
    ax.set_title(f"[{args.label}] Active slip-system count at convergence")
    ax.set_xlabel("FPEND record # (in order)")
    ax.set_ylabel("active systems")
    ax.grid(True, alpha=0.3)
    if fpact.size:
        ax.set_ylim(-0.5, fpact.max() + 0.5)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nsaved: {args.out}")


if __name__ == "__main__":
    main()
