# plot_shear_fig1.py
# Run locally: python plot_shear_fig1.py
# Produces two separate figures matching Mosler & Niehuser Fig.1a and Fig.1b

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

STRESS_CSV = 'shear_stress_F12.csv'
SLIP_CSV   = 'shear_slip_SDV.csv'

for f in [STRESS_CSV, SLIP_CSV]:
    if not os.path.isfile(f):
        print('ERROR: missing ' + f)
        raise SystemExit(1)

stress = np.loadtxt(STRESS_CSV, delimiter=',', skiprows=1)
slip   = np.loadtxt(SLIP_CSV,   delimiter=',', skiprows=1)

F12     = stress[:, 0]
tau12   = stress[:, 1]
F12_sdv = slip[:, 0]
gammas  = slip[:, 1:]   # shape (n_frames, 12)

COLORS = ['#e41a1c','#377eb8','#4daf4a','#984ea3',
          '#ff7f00','#a65628','#f781bf','#999999',
          '#1b9e77','#d95f02','#7570b3','#66a61e']
LS = ['-','-','-','-','-','-','--','--','--','--','--','--']

# ── Fig 1(a): Kirchhoff shear stress ─────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.5, 5))

ax.plot(F12, tau12, color='#1f77b4', lw=2, label='Phase 1')

ax.set_xlabel('Shear deformation $F_{12}$', fontsize=13)
ax.set_ylabel('Kirchhoff stress $\\tau_{12}$ in GPa', fontsize=13)
ax.set_xlim(0, 6)
ax.set_ylim(0, 0.30)
ax.set_yticks(np.arange(0, 0.31, 0.05))
ax.legend(fontsize=11, loc='upper right')
ax.set_title('(a)', fontsize=13, loc='left', fontweight='bold')
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('Fig1a_shear_stress.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: Fig1a_shear_stress.png')

# ── Fig 1(b): Accumulated slips ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.5, 5))

if gammas.max() > 0:
    for a in range(min(12, gammas.shape[1])):
        ax.plot(F12_sdv, gammas[:, a],
                color=COLORS[a], ls=LS[a], lw=1.5,
                label='$\\gamma^{' + str(a+1) + '}$')
    # legend in two columns matching paper layout
    ax.legend(fontsize=8, ncol=2, loc='upper left',
              handlelength=2.5, columnspacing=1.0)

ax.set_xlabel('Shear deformation $F_{12}$', fontsize=13)
ax.set_ylabel('Accumulated slips $\\gamma^{\\alpha}$', fontsize=13)
ax.set_xlim(0, 6)
ax.set_ylim(0, 1.2)
ax.set_yticks(np.arange(0, 1.21, 0.20))
ax.set_title('(b)', fontsize=13, loc='left', fontweight='bold')
ax.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('Fig1b_accumulated_slips.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved: Fig1b_accumulated_slips.png')
