#!/usr/bin/env python3
"""
show_residual.py  -  quick look at the AL residual recorded in the .feh dump.

The UMAT writes two record types:
  NWDUMP  fp_iter  nw_iter  norm_R  tol_nw  norm_step  n_active
  FPDUMP  fp_iter  norm_R   phitilde_norm  kkt(T/F)  n_active  max_absPhi

"Residual" = norm_R = ||R||_2 of the augmented-Lagrangian residual (dimensionless).

Usage:
  python show_residual.py cp_uniaxial_1.feh                 # FP-level (default)
  python show_residual.py cp_uniaxial_1.feh --level newton   # every Newton iter
  python show_residual.py cp_uniaxial_1.feh --level fp       # each fixed-point iter
  python show_residual.py cp_uniaxial_1.feh --print          # also print the table
  python show_residual.py cp_uniaxial_1.feh --out resid.png  # save plot here
  python show_residual.py cp_uniaxial_1.feh --tail 200       # only last 200 records
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def read_feh(path, level):
    """Return a list of (norm_R, tol, n_active, kkt) for the chosen level."""
    rows = []
    tag = "FPDUMP" if level == "fp" else "NWDUMP"
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            p = line.split()
            if not p or p[0] != tag:
                continue
            try:
                if tag == "NWDUMP" and len(p) >= 7:
                    # fp_iter nw_iter norm_R tol_nw norm_step n_active
                    rows.append((float(p[3]), float(p[4]), int(p[6]), None))
                elif tag == "FPDUMP" and len(p) >= 7:
                    # FPDUMP fp_iter norm_R phitilde kkt n_active max_absPhi
                    #   p[0]   p[1]    p[2]   p[3]    p[4] p[5]     p[6]
                    kkt = p[4].strip().upper().startswith("T")
                    rows.append((float(p[2]), None, int(p[5]), kkt))
            except (ValueError, IndexError):
                continue
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh")
    ap.add_argument("--level", choices=["fp", "newton"], default="fp",
                    help="'fp' = residual at each fixed-point iteration (default), "
                         "'newton' = residual at every inner Newton iteration")
    ap.add_argument("--tail", type=int, default=0,
                    help="only show the last N records (0 = all)")
    ap.add_argument("--print", dest="do_print", action="store_true",
                    help="print the residual table to the terminal")
    ap.add_argument("--out", default="residual.png")
    args = ap.parse_args()

    rows = read_feh(args.feh, args.level)
    if not rows:
        print("No matching records found.", file=sys.stderr)
        sys.exit(1)
    if args.tail > 0:
        rows = rows[-args.tail:]

    idx = np.arange(len(rows))
    res = np.array([r[0] for r in rows], dtype=float)
    res_plot = res.copy()
    res_plot[res_plot <= 0] = np.nan  # for log axis

    if args.do_print:
        print(f"{'#':>6}  {'norm_R':>14}  {'n_active':>8}  {'kkt':>4}")
        for i, r in enumerate(rows):
            kkt = "" if r[3] is None else ("T" if r[3] else "F")
            print(f"{i:>6}  {r[0]:>14.6e}  {r[2]:>8}  {kkt:>4}")

    print(f"\nrecords: {len(rows)}   "
          f"min={np.nanmin(res_plot):.3e}   "
          f"median={np.nanmedian(res_plot):.3e}   "
          f"max={np.nanmax(res_plot):.3e}")

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.semilogy(idx, res_plot, "o-", ms=3, lw=0.8, label=r"$\|R\|_2$")
    # mark KKT-satisfied points if available (fp level)
    if args.level == "fp":
        kk = np.array([1.0 if r[3] else 0.0 for r in rows])
        if kk.any():
            ax.semilogy(idx[kk > 0.5], res_plot[kk > 0.5], "g*", ms=8,
                        label="KKT satisfied")
    lvl = "fixed-point iteration" if args.level == "fp" else "inner Newton iteration"
    ax.set_xlabel(f"record # (one per {lvl}, in order)")
    ax.set_ylabel(r"residual $\|R\|_2$")
    ax.set_title(f"AL residual per {lvl}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
