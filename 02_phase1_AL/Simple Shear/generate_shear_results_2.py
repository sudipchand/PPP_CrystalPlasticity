# -*- coding: ascii -*-
# generate_shear_results.py
# Run with: abaqus python generate_shear_results.py
#
# Reads:  Simple_Shear_three_dim.odb
# Writes: shear_stress_F12.csv   -- tau_12 vs F_12  (Mosler Fig.1a)
#         shear_slip_SDV.csv     -- 12 combined slips vs F_12  (Mosler Fig.1b)
#         Fig1_simple_shear.png  -- figure (if matplotlib available)
#
# SDV layout (crystal_plasticity_augm_lagr):
#   SDV 1-9  : Fp_inv tensor components
#   SDV 10-33: plastic_slip(1..24)  accumulated slip on 24 FCC systems
#   SDV 34   : sum(plastic_slip)
#   SDV 35-37: Euler angles of Re
#   SDV 38   : aver_adj_stp
#
# Mosler Fig.1b combined slips:
#   gamma_combined[a] = SDV[9+a] + SDV[9+a+12]  for a=1..12
#
# Geometry:
#   H = 20 mm  (y: -15 to +5),  A = 400 mm^2  (20x20)
#   F_12    = U1_RP / H
#   tau_12  = RF1_RP / A  [MPa] / 1000  [GPa]

import os
import sys

ODB_NAME  = 'Simple_Shear_three_dim.odb'
STEP_NAME = 'displ_contr_test'
RP_NSET   = 'RP'
H         = 20.0    # element height mm
A         = 400.0   # cross-section mm^2
SDV_START = 10      # 1-based index of first slip SDV
N_SYS     = 12      # combined slip systems to plot
OUT_STRESS = 'shear_stress_F12.csv'
OUT_SLIP   = 'shear_slip_SDV.csv'
OUT_FIG    = 'Fig1_simple_shear.png'

# --- open ODB ---------------------------------------------------------------
try:
    from odbAccess import openOdb
except ImportError:
    print 'ERROR: odbAccess not available. Run with: abaqus python'
    sys.exit(1)

if not os.path.isfile(ODB_NAME):
    print 'ERROR: ODB not found: ' + ODB_NAME
    sys.exit(1)

print 'Opening ' + ODB_NAME
odb = openOdb(ODB_NAME, readOnly=True)

if STEP_NAME not in odb.steps:
    print 'ERROR: step not found. Available: ' + str(list(odb.steps.keys()))
    odb.close()
    sys.exit(1)

step   = odb.steps[STEP_NAME]
frames = step.frames
print 'Frames: ' + str(len(frames))

# --- find RP node label -----------------------------------------------------
rootAsm = odb.rootAssembly
inst_name = rootAsm.instances.keys()[0]
inst = rootAsm.instances[inst_name]
print 'Instance: ' + inst_name

rp_node_label = None
if RP_NSET in rootAsm.nodeSets:
    rp_node_label = rootAsm.nodeSets[RP_NSET].nodes[0][0].label
elif RP_NSET in inst.nodeSets:
    rp_node_label = inst.nodeSets[RP_NSET].nodes[0].label
else:
    rp_node_label = 9
    print 'WARNING: RP nset not found, defaulting to node label 9'
print 'RP node label: ' + str(rp_node_label)

# --- frame-by-frame extraction ----------------------------------------------
F12_list  = []
tau12_list = []
slip_list  = []

for i, frame in enumerate(frames):

    # U1 at RP
    u1 = 0.0
    if 'U' in frame.fieldOutputs:
        uf = frame.fieldOutputs['U']
        for v in uf.values:
            if v.nodeLabel == rp_node_label:
                u1 = v.data[0]
                break
    F12 = u1 / H

    # RF1 at RP
    rf1 = 0.0
    if 'RF' in frame.fieldOutputs:
        rf = frame.fieldOutputs['RF']
        for v in rf.values:
            if v.nodeLabel == rp_node_label:
                rf1 = v.data[0]
                break
    tau12 = rf1 / A / 1000.0   # GPa

    # SDV: try named fields SDV10..SDV33 first, then bulk array
    gammas = [0.0] * N_SYS
    named_ok = True
    for a in range(N_SYS):
        key_pos = 'SDV' + str(SDV_START + a)        # SDV10..SDV21
        key_neg = 'SDV' + str(SDV_START + a + 12)   # SDV22..SDV33
        if key_pos in frame.fieldOutputs and key_neg in frame.fieldOutputs:
            pos = frame.fieldOutputs[key_pos].values[0].data
            neg = frame.fieldOutputs[key_neg].values[0].data
            gammas[a] = pos + neg
        else:
            named_ok = False
            break

    if not named_ok and 'SDV' in frame.fieldOutputs:
        sdv_vals = frame.fieldOutputs['SDV'].values[0].data
        for a in range(N_SYS):
            idx_pos = (SDV_START - 1) + a        # 0-based: 9..20
            idx_neg = (SDV_START - 1) + a + 12   # 0-based: 21..32
            if idx_neg < len(sdv_vals):
                gammas[a] = sdv_vals[idx_pos] + sdv_vals[idx_neg]

    F12_list.append(F12)
    tau12_list.append(tau12)
    slip_list.append(gammas[:])

    if i % 50 == 0:
        print 'frame ' + str(i) + '  F12=' + str(round(F12,3)) + \
              '  tau12=' + str(round(tau12,5)) + ' GPa'

odb.close()
print 'ODB closed.'

# --- write stress CSV -------------------------------------------------------
with open(OUT_STRESS, 'w') as f:
    f.write('F12,tau12_GPa\n')
    for F12, tau12 in zip(F12_list, tau12_list):
        f.write(str(F12) + ',' + str(tau12) + '\n')
print 'Written: ' + OUT_STRESS

# --- write slip CSV ---------------------------------------------------------
header = 'F12,' + ','.join(['gamma_' + str(a+1) for a in range(N_SYS)])
with open(OUT_SLIP, 'w') as f:
    f.write(header + '\n')
    for F12, gammas in zip(F12_list, slip_list):
        f.write(str(F12) + ',' + ','.join([str(g) for g in gammas]) + '\n')
print 'Written: ' + OUT_SLIP

# --- plot -------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    F12_arr  = np.array(F12_list)
    t12_arr  = np.array(tau12_list)
    slip_arr = np.array(slip_list)

    COLORS = ['#e41a1c','#377eb8','#4daf4a','#984ea3',
              '#ff7f00','#a65628','#f781bf','#999999',
              '#1b9e77','#d95f02','#7570b3','#66a61e']
    LS = ['-','-','-','-','-','-','--','--','--','--','--','--']

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # panel (a)
    ax = axes[0]
    ax.plot(F12_arr, t12_arr, 'b-', lw=2, label='Phase 1 (Algo 2)')
    ax.set_xlabel('Shear deformation F12', fontsize=13)
    ax.set_ylabel('Kirchhoff stress tau_12 [GPa]', fontsize=13)
    ax.set_xlim(0, max(float(F12_arr.max()), 0.1))
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.set_title('(a)', fontsize=13, loc='left')
    ax.grid(True, alpha=0.3)

    # panel (b)
    ax = axes[1]
    if slip_arr.max() > 0:
        for a in range(N_SYS):
            ax.plot(F12_arr, slip_arr[:, a],
                    color=COLORS[a], ls=LS[a], lw=1.5,
                    label='g' + str(a+1))
        ax.legend(fontsize=8, ncol=2, loc='upper left')
    else:
        ax.text(0.5, 0.5, 'No SDV slip data found.\nCheck SDV in Element Output.',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=11, color='red')
    ax.set_xlabel('Shear deformation F12', fontsize=13)
    ax.set_ylabel('Accumulated slips gamma', fontsize=13)
    ax.set_xlim(0, max(float(F12_arr.max()), 0.1))
    ax.set_ylim(bottom=0)
    ax.set_title('(b)', fontsize=13, loc='left')
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        'Simple shear - Phase 1 (Algo 2) | '
        'kappa=49.98 GPa, G=21.1 GPa, Q0=0.06 GPa, DeltaQ=0.049 GPa, p=10 | '
        'Orientation: 0/0/0 deg',
        fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
    print 'Written: ' + OUT_FIG

except ImportError:
    print 'matplotlib not available -- CSVs written, no figure produced.'
    print 'Copy CSVs locally and plot with standard Python 3 + matplotlib.'
