#!/usr/bin/env python3
"""
plot_newton_per_increment.py - Newton (or fixed-point) iterations per converged
increment, overlaying any number of runs on one figure.

Each CONVRG record closes one converged solve; the NEWTON / FIXPT records since
the previous CONVRG are that solve's counts.

Usage:
  # two runs (back-compatible)
  python plot_newton_per_increment.py base.feh phase1.feh --labels base phase1

  # eta sweep: four (or more) runs overlaid
  python plot_newton_per_increment.py v1.feh v2.feh v3.feh v4.feh \
        --labels 10/C44 smallG 10x 100x

  # fixed-point counts instead of Newton
  python plot_newton_per_increment.py v1.feh v2.feh --metric fixpt

  # also print a summary table
  python plot_newton_per_increment.py v1.feh v2.feh v3.feh --labels a b c --table
"""
import argparse, sys
import numpy as np
import matplotlib.pyplot as plt


def per_solve_counts(path):
    """Return (newton_per_solve, fixpt_per_solve) arrays, one entry per CONVRG."""
    nw, fp = [], []
    nwc = fpc = 0
    seen_convrg = False
    for line in open(path, errors="ignore"):
        p = line.split()
        if not p:
            continue
        if p[0] == "NEWTON":
            nwc += 1
        elif p[0] == "FIXPT":
            fpc += 1
        elif p[0] == "CONVRG":
            nw.append(nwc); fp.append(fpc)
            nwc = fpc = 0
            seen_convrg = True
    return np.array(nw), np.array(fp), seen_convrg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh", nargs="+", help="one or more .feh files to overlay")
    ap.add_argument("--labels", nargs="*", default=None)
    ap.add_argument("--metric", choices=["newton", "fixpt"], default="newton")
    ap.add_argument("--out", default="newton_per_increment.png")
    ap.add_argument("--table", action="store_true",
                    help="print a summary table (mean/median/max + completion)")
    args = ap.parse_args()

    files = args.feh
    labels = args.labels if args.labels else [f"run {k+1}" for k in range(len(files))]
    if len(labels) < len(files):
        labels += [f"run {k+1}" for k in range(len(labels), len(files))]

    ylab = ("Newton iterations per converged increment" if args.metric == "newton"
            else "fixed-point iterations per converged increment")

    # a distinct color cycle that stays readable for up to ~6 runs
    colors = ["tab:blue", "tab:red", "tab:green", "tab:purple",
              "tab:orange", "tab:brown"]

    fig, ax = plt.subplots(figsize=(12, 5.5))
    summary = []
    any_data = False
    for k, (fpath, lab) in enumerate(zip(files, labels)):
        nw, fp, ok = per_solve_counts(fpath)
        series = nw if args.metric == "newton" else fp
        if series.size == 0:
            print(f"  {lab}: no CONVRG records (new-tag run needed)", file=sys.stderr)
            summary.append((lab, 0, np.nan, np.nan, np.nan))
            continue
        any_data = True
        ax.plot(np.arange(series.size), series, "-", lw=0.8,
                color=colors[k % len(colors)],
                label=f"{lab} (mean {series.mean():.1f}, n={series.size})")
        summary.append((lab, series.size, series.mean(),
                        float(np.median(series)), int(series.max())))

    if not any_data:
        print("No CONVRG records found in any file.", file=sys.stderr)
        sys.exit(1)

    ax.set_xlabel("converged increment # (in order)")
    ax.set_ylabel(ylab)
    ax.set_title(f"{ylab}: comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"saved: {args.out}")

    if args.table:
        print()
        print(f"{'run':<14}{'solves':>8}{'mean':>9}{'median':>9}{'max':>8}")
        print("-" * 48)
        for lab, n, mean, med, mx in summary:
            if n == 0:
                print(f"{lab:<14}{'-':>8}{'(no data)':>9}")
            else:
                print(f"{lab:<14}{n:>8}{mean:>9.1f}{med:>9.0f}{mx:>8}")


if __name__ == "__main__":
    main()
