#!/usr/bin/env python3
"""
plot_newton_per_increment.py
Newton (or fixed-point) iterations per converged increment.

File format expected:
  NWEND  <iter>  <residual>  <active_systems>   <- one Newton iteration
  FPEND  <iter>  <residual>  <active_systems>   <- end of converged increment

Usage:
  # compare two runs
  python plot_newton_per_increment.py base.feh phase1.feh --labels base phase1

  # multiple runs (eta sweep etc.)
  python plot_newton_per_increment.py v1.feh v2.feh v3.feh --labels 10/C44 10x 100x

  # fixed-point pass count instead of Newton iters
  python plot_newton_per_increment.py base.feh phase1.feh --metric fixpt

  # also print summary table
  python plot_newton_per_increment.py base.feh phase1.feh --labels base phase1 --table
"""
import argparse, sys
import numpy as np
import matplotlib.pyplot as plt


def per_increment_counts(path):
    """
    Parse a .feh file and return:
      newton_per_inc  : np.array, Newton iterations per converged increment
      fixpt_per_inc   : np.array, fixed-point passes per converged increment
    Tags used:
      NWEND -> one Newton iteration (replaces NEWTON / FIXPT in old format)
      FPEND -> end of a converged increment (replaces CONVRG in old format)
    Also handles old-style tags (NEWTON / FIXPT / CONVRG) as fallback.
    """
    nw_counts, fp_counts = [], []
    nwc = fpc = 0
    new_fp_started = False   # track whether we're inside a new fixed-point pass

    for raw in open(path, errors="ignore"):
        p = raw.split()
        if not p:
            continue
        tag = p[0]

        # ── new-style tags ──────────────────────────────────────────
        if tag == "NWEND":
            nwc += 1
            # iter index resets to 0 at the start of each new fixed-point pass
            try:
                if int(p[1]) == 0:
                    fpc += 1
            except (IndexError, ValueError):
                pass

        elif tag == "FPEND":
            nw_counts.append(nwc)
            fp_counts.append(fpc)
            nwc = fpc = 0

        # ── old-style tags (fallback) ───────────────────────────────
        elif tag == "NEWTON":
            nwc += 1
        elif tag == "FIXPT":
            fpc += 1
        elif tag == "CONVRG":
            nw_counts.append(nwc)
            fp_counts.append(fpc)
            nwc = fpc = 0

    return np.array(nw_counts), np.array(fp_counts)


def main():
    ap = argparse.ArgumentParser(
        description="Plot Newton/fixed-point iterations per converged increment.")
    ap.add_argument("feh", nargs="+", help=".feh file(s) to overlay")
    ap.add_argument("--labels", nargs="*", default=None,
                    help="legend labels, one per file")
    ap.add_argument("--metric", choices=["newton", "fixpt"], default="newton",
                    help="newton = NWEND count, fixpt = fixed-point pass count")
    ap.add_argument("--out", default="newton_per_increment.png",
                    help="output image filename")
    ap.add_argument("--table", action="store_true",
                    help="print a summary table to stdout")
    args = ap.parse_args()

    files  = args.feh
    labels = args.labels or [f"run {k+1}" for k in range(len(files))]
    if len(labels) < len(files):
        labels += [f"run {k+1}" for k in range(len(labels), len(files))]

    ylab = ("Newton iterations per converged increment"
            if args.metric == "newton"
            else "Fixed-point passes per converged increment")

    colors = ["#2166AC", "#D6604D", "#1A9641", "#762A83",
              "#E08214", "#4D9221", "#8073AC"]

    fig, ax = plt.subplots(figsize=(13, 5.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    summary = []
    any_data = False

    for k, (fpath, lab) in enumerate(zip(files, labels)):
        nw, fp = per_increment_counts(fpath)
        series = nw if args.metric == "newton" else fp

        if series.size == 0:
            print(f"  WARNING: '{lab}' has no converged records — check tags.",
                  file=sys.stderr)
            summary.append((lab, 0, float("nan"), float("nan"), 0))
            continue

        any_data = True
        col = colors[k % len(colors)]
        ax.plot(np.arange(series.size), series, "-",
                lw=0.9, color=col, alpha=0.85,
                label=f"{lab}")
        #ax.axhline(series.mean(), color=col, lw=1.4, ls="--", alpha=0.5)

        summary.append((lab, series.size, series.mean(),
                        float(np.median(series)), int(series.max())))

    if not any_data:
        print("ERROR: No converged records found in any file. "
              "Check that your .feh files contain NWEND/FPEND or "
              "NEWTON/FIXPT/CONVRG tags.", file=sys.stderr)
        sys.exit(1)

    ax.set_xlabel("converged increment # (in order)", fontsize=12)
    ax.set_ylabel(ylab, fontsize=12)
    ax.set_title(f"{ylab}", fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(args.out, dpi=160, bbox_inches="tight")
    print(f"Saved: {args.out}")

    if args.table:
        print()
        print(f"{'run':<20} {'n':>7} {'mean':>8} {'median':>8} {'max':>6}")
        print("-" * 54)
        for lab, n, mean, med, mx in summary:
            if n == 0:
                print(f"{lab:<20} {'—':>7}   (no data)")
            else:
                print(f"{lab:<20} {n:>7} {mean:>8.1f} {med:>8.0f} {mx:>6}")


if __name__ == "__main__":
    main()
