# -*- coding: ascii -*-
# generate_shear_results.py
# Run with: abaqus python generate_shear_results.py
#
# Reads:  Simple_Shear_three_dim.odb
# Writes: shear_stress_F12.csv   -- tau_12 vs F_12
#         shear_slip_SDV.csv     -- 12 combined slips vs F_12

import os
import sys

ODB_NAME   = 'Simple_Shear_three_dim.odb'
STEP_NAME  = 'displ_contr_test'
H          = 20.0    # element height mm  (y: -15 to +5)
A          = 400.0   # cross-section mm^2 (20x20)
SDV_START  = 10      # 1-based index of first slip SDV
N_SYS      = 12
OUT_STRESS = 'shear_stress_F12.csv'
OUT_SLIP   = 'shear_slip_SDV.csv'

try:
    from odbAccess import openOdb
except ImportError:
    print 'ERROR: run with abaqus python'
    sys.exit(1)

if not os.path.isfile(ODB_NAME):
    print 'ERROR: ODB not found: ' + ODB_NAME
    sys.exit(1)

print 'Opening ' + ODB_NAME
odb  = openOdb(ODB_NAME, readOnly=True)
step = odb.steps[STEP_NAME]

# --- find RP node label from field output -----------------------------------
rootAsm   = odb.rootAssembly
inst_name = rootAsm.instances.keys()[0]
inst      = rootAsm.instances[inst_name]
print 'Instance: ' + inst_name

# RP is node 9 (the rigid body reference point)
RP_LABEL = 9

# --- read RF1 from history output -------------------------------------------
# History output is the reliable source for RP reaction forces
# Key names in Abaqus history: 'RF1 at Node 9 in NSET RP' or similar
# We scan all history regions to find the one containing RF1

print 'Scanning history output regions...'
hist_regions = step.historyRegions
rf1_hist = None
u1_hist  = None

for key in hist_regions.keys():
    region = hist_regions[key]
    outputs = region.historyOutputs
    # look for RF1
    for okey in outputs.keys():
        if 'RF1' in okey and rf1_hist is None:
            rf1_hist = outputs[okey]
            print 'Found RF1 in region: ' + key + '  key: ' + okey
        if 'U1' in okey and u1_hist is None:
            u1_hist = outputs[okey]
            print 'Found U1 in region: ' + key + '  key: ' + okey
    if rf1_hist is not None and u1_hist is not None:
        break

if rf1_hist is None:
    print 'WARNING: RF1 not found in history output.'
    print 'Available history regions:'
    for key in hist_regions.keys():
        print '  ' + key + ': ' + str(hist_regions[key].historyOutputs.keys())

if u1_hist is None:
    print 'WARNING: U1 not found in history output.'

# history data is list of (time, value) tuples
if u1_hist is not None:
    u1_data  = u1_hist.data    # [(time, u1), ...]
else:
    u1_data = []

if rf1_hist is not None:
    rf1_data = rf1_hist.data   # [(time, rf1), ...]
else:
    rf1_data = []

print 'U1  history points: ' + str(len(u1_data))
print 'RF1 history points: ' + str(len(rf1_data))

# build dicts keyed by time for alignment
u1_dict  = dict(u1_data)
rf1_dict = dict(rf1_data)

# --- read SDV from field output frame by frame ------------------------------
frames = step.frames
print 'Field output frames: ' + str(len(frames))

F12_list   = []
tau12_list = []
slip_list  = []

for i, frame in enumerate(frames):
    t = frame.frameValue

    # U1 from history dict (align by time)
    u1 = u1_dict.get(t, None)
    if u1 is None:
        # fallback: scan field output U
        if 'U' in frame.fieldOutputs:
            for v in frame.fieldOutputs['U'].values:
                if v.nodeLabel == RP_LABEL:
                    u1 = v.data[0]
                    break
    if u1 is None:
        u1 = 0.0
    F12 = u1 / H

    # RF1 from history dict
    rf1 = rf1_dict.get(t, 0.0)
    tau12 = rf1 / A / 1000.0   # MPa -> GPa

    # SDV: try named fields first, then bulk array
    gammas = [0.0] * N_SYS
    named_ok = True
    for a in range(N_SYS):
        kp = 'SDV' + str(SDV_START + a)
        kn = 'SDV' + str(SDV_START + a + 12)
        if kp in frame.fieldOutputs and kn in frame.fieldOutputs:
            pos = frame.fieldOutputs[kp].values[0].data
            neg = frame.fieldOutputs[kn].values[0].data
            gammas[a] = pos + neg
        else:
            named_ok = False
            break

    if not named_ok and 'SDV' in frame.fieldOutputs:
        sdv_vals = frame.fieldOutputs['SDV'].values[0].data
        for a in range(N_SYS):
            ip = (SDV_START - 1) + a
            iq = (SDV_START - 1) + a + 12
            if iq < len(sdv_vals):
                gammas[a] = sdv_vals[ip] + sdv_vals[iq]

    F12_list.append(F12)
    tau12_list.append(tau12)
    slip_list.append(gammas[:])

    if i % 20 == 0:
        print 'frame ' + str(i) + '  F12=' + str(round(F12,3)) + \
              '  tau12=' + str(round(tau12,5)) + ' GPa' + \
              '  gamma_sum=' + str(round(sum(gammas),4))

odb.close()
print 'ODB closed.'

# --- write CSVs -------------------------------------------------------------
with open(OUT_STRESS, 'w') as f:
    f.write('F12,tau12_GPa\n')
    for F12, tau12 in zip(F12_list, tau12_list):
        f.write(str(F12) + ',' + str(tau12) + '\n')
print 'Written: ' + OUT_STRESS

header = 'F12,' + ','.join(['gamma_' + str(a+1) for a in range(N_SYS)])
with open(OUT_SLIP, 'w') as f:
    f.write(header + '\n')
    for F12, gammas in zip(F12_list, slip_list):
        f.write(str(F12) + ',' + ','.join([str(g) for g in gammas]) + '\n')
print 'Written: ' + OUT_SLIP
print 'Done. Copy CSVs locally and run plot_shear_fig1.py to generate figure.'
