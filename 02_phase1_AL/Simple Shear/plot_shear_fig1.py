# plot_shear_fig1.py
# Run locally with Python 3 after copying CSVs from cluster:
#   python plot_shear_fig1.py
#
# Reads:  shear_stress_F12.csv
#         shear_slip_SDV.csv
# Writes: Fig1_simple_shear.png

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

STRESS_CSV = 'shear_stress_F12.csv'
SLIP_CSV   = 'shear_slip_SDV.csv'
OUT_FIG    = 'Fig1_simple_shear.png'

for f in [STRESS_CSV, SLIP_CSV]:
    if not os.path.isfile(f):
        print('ERROR: missing ' + f)
        raise SystemExit(1)

stress = np.loadtxt(STRESS_CSV, delimiter=',', skiprows=1)
slip   = np.loadtxt(SLIP_CSV,   delimiter=',', skiprows=1)

F12      = stress[:, 0]
tau12    = stress[:, 1]
F12_sdv  = slip[:, 0]
gammas   = slip[:, 1:]   # shape (n_frames, 12)

# colours matching Mosler paper Fig.1b legend order
COLORS = ['#e41a1c','#377eb8','#4daf4a','#984ea3',
          '#ff7f00','#a65628','#f781bf','#999999',
          '#1b9e77','#d95f02','#7570b3','#66a61e']
LS     = ['-','-','-','-','-','-','--','--','--','--','--','--']

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# --- panel (a): tau_12 vs F_12 ----------------------------------------------
ax = axes[0]
ax.plot(F12, tau12, color='#1f77b4', lw=2, label='Phase 1 (Algo 2)')
ax.set_xlabel('Shear deformation $F_{12}$', fontsize=13)
ax.set_ylabel('Kirchhoff stress $\\tau_{12}$ in GPa', fontsize=13)
ax.set_xlim(0, max(F12.max(), 0.5))
ax.set_ylim(bottom=0)
ax.legend(fontsize=11)
ax.set_title('(a)', fontsize=13, loc='left')
ax.grid(True, alpha=0.3)

# --- panel (b): accumulated slips -------------------------------------------
ax = axes[1]
if gammas.max() > 0:
    for a in range(min(12, gammas.shape[1])):
        ax.plot(F12_sdv, gammas[:, a],
                color=COLORS[a], ls=LS[a], lw=1.5,
                label='$\\gamma^{' + str(a+1) + '}$')
    ax.legend(fontsize=8, ncol=2, loc='upper left')
else:
    ax.text(0.5, 0.5, 'No slip data (SDV all zero).\nCheck SDV field output in inp.',
            ha='center', va='center', transform=ax.transAxes,
            fontsize=11, color='red')

ax.set_xlabel('Shear deformation $F_{12}$', fontsize=13)
ax.set_ylabel('Accumulated slips $\\gamma^{\\alpha}$', fontsize=13)
ax.set_xlim(0, max(F12_sdv.max(), 0.5))
ax.set_ylim(bottom=0)
ax.set_title('(b)', fontsize=13, loc='left')
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
print('Saved: ' + OUT_FIG)
