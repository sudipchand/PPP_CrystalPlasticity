#!/usr/bin/env python3
"""
plot_eta_vs_active.py - overlay the penalty eta and the active slip-system count
per converged increment, to test why the active set fluctuates.

Requires the eta-instrumented UMAT, whose CONVRG record is:
  CONVRG  l_iter_fp  ||R||  n_active  eta

Reads CONVRG records and plots, on a shared x-axis (converged-increment order):
  - active slip-system count (left axis)
  - penalty eta (right axis, log scale)

If eta spikes exactly where the active count starts fluctuating, the fluctuation
is penalty-driven (Candidate 2). If eta is steady while the count fluctuates,
it points to Taylor-ambiguity / reorientation instead.

Usage:
  python plot_eta_vs_active.py phase1.feh --label phase1 --out eta_active.png
  python plot_eta_vs_active.py phase1.feh --tail 600
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read_convrg(path):
    act, eta = [], []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 5 and p[0] == "CONVRG":
                try:
                    act.append(int(p[3]))
                    eta.append(float(p[4]))
                except (ValueError, IndexError):
                    continue
    return np.array(act), np.array(eta, dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh")
    ap.add_argument("--label", default="run")
    ap.add_argument("--tail", type=int, default=0)
    ap.add_argument("--out", default="eta_active.png")
    args = ap.parse_args()

    act, eta = read_convrg(args.feh)
    if act.size == 0:
        print("No CONVRG records with eta found. Is this the eta-instrumented run?",
              file=sys.stderr)
        sys.exit(1)
    if args.tail > 0:
        act, eta = act[-args.tail:], eta[-args.tail:]

    x = np.arange(act.size)
    eta_plot = eta.copy(); eta_plot[eta_plot <= 0] = np.nan

    fig, ax1 = plt.subplots(figsize=(12, 5))
    ax1.plot(x, act, "-", color="tab:green", lw=0.8, label="active systems")
    ax1.set_xlabel("converged increment # (in order)")
    ax1.set_ylabel("active slip systems", color="tab:green")
    ax1.tick_params(axis="y", labelcolor="tab:green")
    ax1.set_ylim(-0.5, act.max() + 0.5)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.semilogy(x, eta_plot, "-", color="tab:purple", lw=0.8, label="eta")
    ax2.set_ylabel("penalty eta (log)", color="tab:purple")
    ax2.tick_params(axis="y", labelcolor="tab:purple")

    plt.title(f"[{args.label}] active systems vs penalty eta per converged increment")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"saved: {args.out}")
    print(f"eta: median={np.nanmedian(eta_plot):.3e}  min={np.nanmin(eta_plot):.3e}  max={np.nanmax(eta_plot):.3e}")
    print(f"active: at 8 = {(act==8).sum()} of {act.size}")


if __name__ == "__main__":
    main()
