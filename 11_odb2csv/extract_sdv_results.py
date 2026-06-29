from odbAccess import *

odb = openOdb(path='cp_uniaxial_1.odb')
step = odb.steps['TENSION']

out = open('sdv_results.csv', 'w')

header = ['frame', 'time', 'S33', 'LE33']
for i in range(1, 39):
    header.append('SDV%d' % i)

out.write(','.join(header) + '\n')

for frame_id, frame in enumerate(step.frames):

    time = frame.frameValue

    S = frame.fieldOutputs['S'].values[0].data
    LE = frame.fieldOutputs['LE'].values[0].data

    S33 = S[2]
    LE33 = LE[2]

    row = [str(frame_id), str(time), str(S33), str(LE33)]

    for i in range(1, 39):
        key = 'SDV%d' % i
        val = frame.fieldOutputs[key].values[0].data
        row.append(str(val))

    out.write(','.join(row) + '\n')

out.close()
odb.close()

print('Created sdv_results.csv')
