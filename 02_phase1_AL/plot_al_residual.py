#!/usr/bin/env python3
"""
plot_al_residual.py

Parse the augmented-Lagrangian diagnostic dump written to the UMAT's unit-16
file (cp_uniaxial_1.feh) and plot the convergence behaviour of the inner
Newton loop and the outer fixed-point loop.

The UMAT writes two record types:

  NWDUMP  l_iter_fp  l_iter_nw  norm_R       tol_nw       norm_step     n_active
  FPDUMP  l_iter_fp  norm_R     phi_tilde_n  kkt_ok(L)    n_active      max_absPhi

where
  l_iter_fp   = outer fixed-point iteration counter
  l_iter_nw   = inner Newton iteration counter
  norm_R      = ||R||_2 of the AL residual (dimensionless)
  tol_nw      = current adaptive inner-Newton tolerance
  norm_step   = ||delta(delta_lambda)||_2 (the Newton step size)
  n_active    = number of active slip systems
  phi_tilde_n = ||Phi_tilde_+||_inf normalised by tau0 (dimensionless)
  kkt_ok      = T/F, whether the KKT exit condition was satisfied
  max_absPhi  = max_alpha |Phi^alpha| in stress units (MPa)

Because the UMAT has no access to the Abaqus increment number, the dump is a
flat stream of all calls in order. A new UMAT call (material point evaluation)
is detected whenever the fixed-point counter resets to a low value after having
been higher, OR whenever an NWDUMP with l_iter_fp smaller than the previous
record appears. The LAST detected "call block" corresponds to the final
(failing) attempt of the last increment, which is usually what we want to see.

Usage:
    python plot_al_residual.py cp_uniaxial_1.feh
    python plot_al_residual.py cp_uniaxial_1.feh --block last
    python plot_al_residual.py cp_uniaxial_1.feh --block all
    python plot_al_residual.py cp_uniaxial_1.feh --out diag.png
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt


def parse_feh(path):
    """Read the dump file and return lists of NWDUMP and FPDUMP records."""
    nw = []  # (fp, nw, norm_R, tol_nw, norm_step, n_active)
    fp = []  # (fp, norm_R, phi_tilde_n, kkt_ok, n_active, max_absPhi)
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            tag = parts[0]
            try:
                if tag == "NWDUMP" and len(parts) >= 7:
                    nw.append((
                        int(parts[1]), int(parts[2]),
                        float(parts[3]), float(parts[4]),
                        float(parts[5]), int(parts[6]),
                    ))
                elif tag == "FPDUMP" and len(parts) >= 7:
                    kkt = parts[4].strip().upper().startswith("T")
                    fp.append((
                        int(parts[1]), float(parts[2]),
                        float(parts[3]), kkt,
                        int(parts[5]), float(parts[6]),
                    ))
            except (ValueError, IndexError):
                # skip any non-conforming line (e.g. real warning text)
                continue
    return nw, fp


def split_blocks(nw):
    """
    Split the flat NWDUMP stream into separate UMAT-call blocks.
    A new block starts whenever the fixed-point counter does not increase
    monotonically (i.e. it resets to a value <= the running max for that block
    while the Newton counter is also 0), which marks a fresh material-point
    evaluation / attempt.
    """
    blocks = []
    cur = []
    prev_fp = None
    for rec in nw:
        fp_i, nw_i = rec[0], rec[1]
        if prev_fp is not None and fp_i < prev_fp and nw_i == 0:
            # fixed-point counter went backwards at the start of a Newton
            # sweep -> new call block
            if cur:
                blocks.append(cur)
            cur = []
        cur.append(rec)
        prev_fp = fp_i
    if cur:
        blocks.append(cur)
    return blocks


def plot_block(nw_block, fp_all, title_suffix, out_path):
    """Produce the 3-panel diagnostic figure for one call block."""
    nw_arr = np.array([(i, r[2], r[3], r[4], r[5])
                       for i, r in enumerate(nw_block)], dtype=float)
    # columns: seq_index, norm_R, tol_nw, norm_step, n_active
    seq = nw_arr[:, 0]
    norm_R = nw_arr[:, 1]
    tol_nw = nw_arr[:, 2]
    norm_step = nw_arr[:, 3]
    n_active = nw_arr[:, 4]

    # guard against zeros for log plots
    def safe(x):
        x = np.array(x, dtype=float)
        x[x <= 0] = np.nan
        return x

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=False)

    # Panel 1: residual vs tolerance over the whole call block
    ax = axes[0]
    ax.semilogy(seq, safe(norm_R), "o-", ms=3, label=r"$\|R\|_2$ (residual)")
    ax.semilogy(seq, safe(tol_nw), "--", label=r"$tol_{nw}$ (adaptive tol.)")
    ax.semilogy(seq, safe(norm_step), ":", label=r"$\|\delta(\Delta\lambda)\|$ (step)")
    ax.set_xlabel("diagnostic record # within call block")
    ax.set_ylabel("magnitude")
    ax.set_title("Inner-Newton residual, tolerance and step size " + title_suffix)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    # Panel 2: number of active slip systems
    ax = axes[1]
    ax.plot(seq, n_active, "s-", ms=3, color="tab:green")
    ax.set_xlabel("diagnostic record # within call block")
    ax.set_ylabel("active slip systems")
    ax.set_title("Active slip-system count")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.5, max(n_active.max(), 1) + 0.5)

    # Panel 3: outer fixed-point convergence (from FPDUMP, full history)
    ax = axes[2]
    if fp_all:
        fp_arr = np.array([(i, r[1], r[2], 1.0 if r[3] else 0.0, r[5])
                           for i, r in enumerate(fp_all)], dtype=float)
        fseq = fp_arr[:, 0]
        f_normR = fp_arr[:, 1]
        f_phitilde = fp_arr[:, 2]
        f_kkt = fp_arr[:, 3]
        ax.semilogy(fseq, safe(f_normR), "o-", ms=2,
                    label=r"$\|R\|_2$ at FP exit")
        ax.semilogy(fseq, safe(f_phitilde), ".-", ms=2,
                    label=r"$\|\tilde\Phi_+\|_\infty / \tau_0$ (norm.)")
        # mark KKT-pass points
        passes = fseq[f_kkt > 0.5]
        if passes.size:
            ax.semilogy(passes, safe(f_phitilde)[f_kkt > 0.5], "g*",
                        ms=10, label="KKT satisfied")
        ax.set_xlabel("FPDUMP record # (all calls, in order)")
        ax.set_ylabel("magnitude")
        ax.set_title("Outer fixed-point convergence (whole run)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    else:
        ax.text(0.5, 0.5, "no FPDUMP records found",
                ha="center", va="center", transform=ax.transAxes)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"saved figure: {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("feh", help="path to the .feh diagnostic dump")
    ap.add_argument("--block", default="last",
                    help="'last' (default), 'all', or an integer block index")
    ap.add_argument("--out", default="al_diagnostic.png",
                    help="output PNG path")
    args = ap.parse_args()

    nw, fp = parse_feh(args.feh)
    if not nw:
        print("No NWDUMP records parsed - is this the right file?",
              file=sys.stderr)
        sys.exit(1)

    blocks = split_blocks(nw)
    print(f"parsed {len(nw)} NWDUMP and {len(fp)} FPDUMP records "
          f"in {len(blocks)} call blocks")

    # quick textual summary of the final block (the failing attempt)
    last = blocks[-1]
    fp_last = max(r[0] for r in last)
    nw_max = max(r[1] for r in last)
    minR = min(r[2] for r in last)
    print(f"final call block: {len(last)} records, "
          f"max fixed-point iter = {fp_last}, max inner Newton iter = {nw_max}, "
          f"min ||R|| = {minR:.3e}")
    # how often did the inner Newton hit its cap?
    capped = sum(1 for r in last if r[1] >= 20)
    print(f"  records at inner-Newton iter >= 20 (cap): {capped}")

    if args.block == "all":
        for bi, blk in enumerate(blocks):
            plot_block(blk, fp, f"(block {bi})",
                       args.out.replace(".png", f"_block{bi}.png"))
    elif args.block == "last":
        plot_block(last, fp, "(final / failing call block)", args.out)
    else:
        bi = int(args.block)
        plot_block(blocks[bi], fp, f"(block {bi})", args.out)


if __name__ == "__main__":
    main()
