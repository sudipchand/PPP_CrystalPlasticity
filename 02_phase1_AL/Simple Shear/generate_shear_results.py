"""
generate_shear_results.py
=========================
Abaqus Python script -- run with:
    abaqus python generate_shear_results.py

Reads:  Simple_Shear_three_dim.odb
Writes: shear_stress_F12.csv        -- tau_12 vs F_12  (Panel a)
        shear_slip_SDV.csv           -- 12 combined slips vs F_12  (Panel b)
        Fig1_simple_shear.png        -- both panels (requires matplotlib)

SDV layout (from crystal_plasticity_augm_lagr):
  SDV 1-9  : Fp_inv components
  SDV 10-33: plastic_slip(1..24)  -- accumulated slip on 24 FCC systems
  SDV 34   : sum(plastic_slip)
  SDV 35-37: Euler angles of Re
  SDV 38   : aver_adj_stp

Mosler Fig.1(b) combined slips:
  gamma_combined[alpha] = SDV[9+alpha] + SDV[9+alpha+12]  for alpha=1..12
  i.e. positive + negative direction on same slip system

Geometry:
  Element height H = 20 mm  (y: -15 to +5)
  F_12 = U1_RP / H
  tau_12 = RF1_RP / A  (in MPa, then /1000 -> GPa)
  A = 20 mm x 20 mm = 400 mm^2  (x: -10 to +10, z: 0 to 20)
"""

import os
import sys

# ── parameters ────────────────────────────────────────────────────────────────
ODB_NAME    = 'Simple_Shear_three_dim.odb'
STEP_NAME   = 'displ_contr_test'
RP_NSET     = 'RP'                  # nset containing the reference point node
ELEMENT_H   = 20.0                  # mm
ELEMENT_A   = 400.0                 # mm^2  (cross-section area)
SDV_SLIP_START = 10                 # 1-based index of first slip SDV (SDV10)
N_SYSTEMS   = 12                    # Mosler plots 12 combined slip systems

OUT_STRESS  = 'shear_stress_F12.csv'
OUT_SLIP    = 'shear_slip_SDV.csv'
OUT_FIG     = 'Fig1_simple_shear.png'

# ── open ODB ─────────────────────────────────────────────────────────────────
try:
    from odbAccess import openOdb
    odb_available = True
except ImportError:
    odb_available = False

if not odb_available:
    print('ERROR: odbAccess not available. Run with: abaqus python generate_shear_results.py')
    sys.exit(1)

if not os.path.isfile(ODB_NAME):
    print('ERROR: ODB not found: ' + ODB_NAME)
    sys.exit(1)

print('Opening ' + ODB_NAME + ' ...')
odb = openOdb(ODB_NAME, readOnly=True)

if STEP_NAME not in odb.steps:
    print('ERROR: step "' + STEP_NAME + '" not in ODB.')
    print('Available steps: ' + str(list(odb.steps.keys())))
    odb.close()
    sys.exit(1)

step = odb.steps[STEP_NAME]
frames = step.frames
n_frames = len(frames)
print('Number of frames: ' + str(n_frames))

# ── find RP node ──────────────────────────────────────────────────────────────
rootAssembly = odb.rootAssembly
instance_name = list(rootAssembly.instances.keys())[0]
instance = rootAssembly.instances[instance_name]
print('Instance: ' + instance_name)

# Find the RP nset -- try assembly level first, then instance level
rp_node_label = None
if RP_NSET in rootAssembly.nodeSets:
    rp_nset = rootAssembly.nodeSets[RP_NSET]
    rp_node_label = rp_nset.nodes[0][0].label
elif RP_NSET in instance.nodeSets:
    rp_nset = instance.nodeSets[RP_NSET]
    rp_node_label = rp_nset.nodes[0].label
else:
    # RP is node 9 in the inp (the reference point)
    print('WARNING: nset "' + RP_NSET + '" not found, using node label 9')
    rp_node_label = 9

print('RP node label: ' + str(rp_node_label))

# ── element label for SDV extraction ─────────────────────────────────────────
elem_label = instance.elements[0].label
print('Element label: ' + str(elem_label))

# ── extract data frame by frame ───────────────────────────────────────────────
F12_list   = []
tau12_list = []
# slip data: list of lists, one per frame
slip_list  = []   # each entry: [gamma_1..gamma_12]

print('Extracting frames ...')
for i, frame in enumerate(frames):
    frame_time = frame.frameValue

    # ── U1 at RP ────────────────────────────────────────────────────────────
    u1 = 0.0
    if 'U' in frame.fieldOutputs:
        u_field = frame.fieldOutputs['U']
        # get value at RP node
        try:
            u_val = u_field.getSubset(region=rootAssembly.nodeSets[RP_NSET]) \
                           if RP_NSET in rootAssembly.nodeSets \
                           else u_field.getSubset(region=instance.nodeSets[RP_NSET])
            u1 = u_val.values[0].data[0]   # component 1 = U1
        except Exception:
            # fallback: scan all values for RP node label
            for v in u_field.values:
                if v.nodeLabel == rp_node_label:
                    u1 = v.data[0]
                    break

    F12 = u1 / ELEMENT_H

    # ── RF1 at RP ───────────────────────────────────────────────────────────
    rf1 = 0.0
    if 'RF' in frame.fieldOutputs:
        rf_field = frame.fieldOutputs['RF']
        try:
            rf_val = rf_field.getSubset(region=rootAssembly.nodeSets[RP_NSET]) \
                             if RP_NSET in rootAssembly.nodeSets \
                             else rf_field.getSubset(region=instance.nodeSets[RP_NSET])
            rf1 = rf_val.values[0].data[0]
        except Exception:
            for v in rf_field.values:
                if v.nodeLabel == rp_node_label:
                    rf1 = v.data[0]
                    break

    # Kirchhoff tau_12: RF1 / A [MPa] / 1000 -> GPa
    # For simple shear J ≈ 1 so Kirchhoff ≈ Cauchy
    tau12 = rf1 / ELEMENT_A / 1000.0

    F12_list.append(F12)
    tau12_list.append(tau12)

    # ── SDV at integration point 1 ──────────────────────────────────────────
    gammas = [0.0] * N_SYSTEMS
    if 'SDV' in frame.fieldOutputs:
        # SDV is a scalar field -- need individual SDV components
        # Try SDV10..SDV33 directly
        all_sdv_ok = True
        pos_slips = []
        neg_slips = []
        for alpha in range(1, N_SYSTEMS + 1):
            sdv_pos_key = 'SDV' + str(SDV_SLIP_START - 1 + alpha)        # SDV10..SDV21
            sdv_neg_key = 'SDV' + str(SDV_SLIP_START - 1 + alpha + 12)   # SDV22..SDV33
            pos_val = 0.0
            neg_val = 0.0
            if sdv_pos_key in frame.fieldOutputs:
                for v in frame.fieldOutputs[sdv_pos_key].values:
                    pos_val = v.data; break
            else:
                all_sdv_ok = False
            if sdv_neg_key in frame.fieldOutputs:
                for v in frame.fieldOutputs[sdv_neg_key].values:
                    neg_val = v.data; break
            else:
                all_sdv_ok = False
            gammas[alpha - 1] = pos_val + neg_val

        if not all_sdv_ok and i == 0:
            # Try bulk SDV field (older Abaqus: all SDVs in one array field)
            sdv_field = frame.fieldOutputs['SDV']
            if hasattr(sdv_field.values[0], 'dataDouble'):
                sdv_arr = list(sdv_field.values[0].dataDouble)
            else:
                sdv_arr = list(sdv_field.values[0].data)
            if len(sdv_arr) >= 33:
                for alpha in range(N_SYSTEMS):
                    idx_pos = SDV_SLIP_START - 1 + alpha        # 0-based: 9..20
                    idx_neg = SDV_SLIP_START - 1 + alpha + 12   # 0-based: 21..32
                    gammas[alpha] = sdv_arr[idx_pos] + sdv_arr[idx_neg]
                print('Using bulk SDV array, length=' + str(len(sdv_arr)))

    slip_list.append(gammas[:])

    if i % 50 == 0:
        print('  frame ' + str(i) + '/' + str(n_frames) +
              '  F12=' + str(round(F12, 3)) +
              '  tau12=' + str(round(tau12, 5)) + ' GPa')

odb.close()
print('ODB closed.')

# ── write CSV: stress ─────────────────────────────────────────────────────────
with open(OUT_STRESS, 'w') as f:
    f.write('F12,tau12_GPa\n')
    for F12, tau12 in zip(F12_list, tau12_list):
        f.write(str(F12) + ',' + str(tau12) + '\n')
print('Written: ' + OUT_STRESS)

# ── write CSV: slips ──────────────────────────────────────────────────────────
header = 'F12,' + ','.join(['gamma_' + str(a+1) for a in range(N_SYSTEMS)])
with open(OUT_SLIP, 'w') as f:
    f.write(header + '\n')
    for F12, gammas in zip(F12_list, slip_list):
        row = str(F12) + ',' + ','.join([str(g) for g in gammas])
        f.write(row + '\n')
print('Written: ' + OUT_SLIP)

# ── plot ──────────────────────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    F12_arr   = np.array(F12_list)
    tau12_arr = np.array(tau12_list)
    slip_arr  = np.array(slip_list)   # shape (n_frames, 12)

    # colours matching Mosler paper legend order (gamma1..gamma12)
    COLORS = [
        '#e41a1c','#377eb8','#4daf4a','#984ea3',
        '#ff7f00','#a65628','#f781bf','#999999',
        '#1b9e77','#d95f02','#7570b3','#66a61e',
    ]
    LS = ['-','-','-','-','-','-','--','--','--','--','--','--']

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel (a)
    ax = axes[0]
    ax.plot(F12_arr, tau12_arr, 'b-', lw=2, label='Phase 1 (Algo 2)')
    ax.set_xlabel(r'Shear deformation $F_{12}$', fontsize=13)
    ax.set_ylabel(r'Kirchhoff stress $\tau_{12}$ in GPa', fontsize=13)
    ax.set_xlim(0, max(F12_arr.max(), 0.1))
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=11)
    ax.set_title('(a)', fontsize=13, loc='left')
    ax.grid(True, alpha=0.3)

    # Panel (b)
    ax = axes[1]
    if slip_arr.shape[1] == N_SYSTEMS and slip_arr.max() > 0:
        for alpha in range(N_SYSTEMS):
            ax.plot(F12_arr, slip_arr[:, alpha],
                    color=COLORS[alpha], ls=LS[alpha], lw=1.5,
                    label=r'$\gamma^{' + str(alpha+1) + r'}$')
        ax.legend(fontsize=8, ncol=2, loc='upper left')
    else:
        ax.text(0.5, 0.5, 'No SDV slip data extracted.\nCheck SDV field output in ODB.',
                ha='center', va='center', transform=ax.transAxes,
                fontsize=11, color='red')
    ax.set_xlabel(r'Shear deformation $F_{12}$', fontsize=13)
    ax.set_ylabel(r'Accumulated slips $\gamma^\alpha$', fontsize=13)
    ax.set_xlim(0, max(F12_arr.max(), 0.1))
    ax.set_ylim(bottom=0)
    ax.set_title('(b)', fontsize=13, loc='left')
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        'Simple shear – Phase 1 (Algo 2)  |  '
        r'$\kappa=49.98$ GPa, $G=21.1$ GPa, $Q_0=0.06$ GPa, $\Delta Q=0.049$ GPa, $p=10$'
        '\nOrientation: 0°/0°/0°',
        fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
    print('Written: ' + OUT_FIG)

except ImportError:
    print('matplotlib not available -- CSVs written, no figure produced.')
    print('Copy CSVs to local machine and run plot_from_csv.py to plot.')
