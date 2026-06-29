#!/usr/bin/env python3
"""
plot_tangent_condition.py - plot the condition number of the consistent tangent
DDSDDE vs increment, to evidence the active-set / tangent instability.

Mirrors Niehueser & Mosler Fig. 2(b) (condition number of the Jacobian).

The UMAT writes, each accepted increment:
  DDSDDE  KINC  NTENS  d11 d12 ... d16 d21 ... d66   (row-major NTENS*NTENS)

For each record we form the NTENS x NTENS matrix and compute its 2-norm
condition number (ratio of largest to smallest singular value). A large /
spiking condition number indicates a near-singular or ill-conditioned tangent,
which is what an unstable active set produces.

Usage:
  python plot_tangent_condition.py phase1.feh --label phase1
  python plot_tangent_condition.py base.feh phase1.feh --labels base phase1
  python plot_tangent_condition.py phase1.feh --out cond.png
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read_cond(path):
    """Return (kinc_array, condnum_array) from DDSDDE records."""
    kincs, conds = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            p = line.split()
            if not p or p[0] != "DDSDDE":
                continue
            try:
                kinc = int(p[1])
                ntens = int(p[2])
                vals = [float(x) for x in p[3:3 + ntens * ntens]]
                if len(vals) != ntens * ntens:
                    continue
                M = np.array(vals, dtype=float).reshape(ntens, ntens)
                # 2-norm condition number via singular values
                sv = np.linalg.svd(M, compute_uv=False)
                smin = sv[-1]
                cond = sv[0] / smin if smin > 0 else np.inf
                kincs.append(kinc)
                conds.append(cond)
            except (ValueError, IndexError, np.linalg.LinAlgError):
                continue
    return np.array(kincs), np.array(conds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh", nargs="+", help="one or two .feh files")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--out", default="tangent_condition.png")
    args = ap.parse_args()

    files = args.feh
    labels = args.labels if args.labels else [f"run {i+1}" for i in range(len(files))]
    if len(labels) < len(files):
        labels += [f"run {i+1}" for i in range(len(labels), len(files))]

    fig, ax = plt.subplots(figsize=(11, 5))
    colors = ["tab:blue", "tab:red", "tab:green"]
    any_data = False
    for i, (fpath, lab) in enumerate(zip(files, labels)):
        kinc, cond = read_cond(fpath)
        if kinc.size == 0:
            print(f"  {lab}: no DDSDDE records in {fpath}", file=sys.stderr)
            continue
        any_data = True
        # cond can be plotted vs record order (multiple GP per KINC -> use index)
        x = np.arange(cond.size)
        ax.semilogy(x, cond, "-", lw=0.8, color=colors[i % len(colors)], label=lab)
        finite = cond[np.isfinite(cond)]
        print(f"  {lab}: {cond.size} tangents  "
              f"median cond={np.median(finite):.2e}  max cond={finite.max():.2e}")

    if not any_data:
        print("No DDSDDE records found in any file.", file=sys.stderr)
        sys.exit(1)

    ax.set_xlabel("DDSDDE record # (in order)")
    ax.set_ylabel("condition number of DDSDDE")
    ax.set_title("Consistent-tangent condition number")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
