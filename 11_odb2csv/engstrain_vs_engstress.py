from odbAccess import *

odb = openOdb(path='cp_uniaxial_1.odb')
step = odb.steps['TENSION']

top_nodes = [
    'Node PART-1-1.5',
    'Node PART-1-1.6',
    'Node PART-1-1.7',
    'Node PART-1-1.8'
]

L0 = 1.0
A0 = 1.0

first_region = step.historyRegions[top_nodes[0]]
u_data = first_region.historyOutputs['U3'].data

csv_name = 'stress_strain_force_displacement.csv'
plot_script = 'plot_main_results.py'

with open(csv_name, 'w') as f:
    f.write('time,U3_avg,RF3_total,engineering_strain,engineering_stress\n')

    for i in range(len(u_data)):
        time = u_data[i][0]

        u_sum = 0.0
        rf_sum = 0.0

        for node in top_nodes:
            region = step.historyRegions[node]
            u_sum += region.historyOutputs['U3'].data[i][1]
            rf_sum += region.historyOutputs['RF3'].data[i][1]

        u_avg = u_sum / len(top_nodes)

        eng_strain = u_avg / L0
        eng_stress = rf_sum / A0

        f.write('%g,%g,%g,%g,%g\n' %
                (time, u_avg, rf_sum, eng_strain, eng_stress))

odb.close()

# Create separate normal Python plotting script
with open(plot_script, 'w') as f:
    f.write("""import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

time = []
u = []
rf = []
strain = []
stress = []

with open('stress_strain_force_displacement.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        time.append(float(row['time']))
        u.append(float(row['U3_avg']))
        rf.append(float(row['RF3_total']))
        strain.append(float(row['engineering_strain']))
        stress.append(float(row['engineering_stress']))

plt.figure()
plt.plot(u, rf, marker='o')
plt.xlabel('Average displacement U3')
plt.ylabel('Total reaction force RF3')
plt.title('Force-Displacement Curve')
plt.grid(True)
plt.tight_layout()
plt.savefig('force_displacement.png', dpi=300)
plt.close()

plt.figure()
plt.plot(strain, stress, marker='o')
plt.xlabel('Engineering strain')
plt.ylabel('Engineering stress')
plt.title('Engineering Stress-Strain Curve')
plt.grid(True)
plt.tight_layout()
plt.savefig('stress_strain.png', dpi=300)
plt.close()

print('Created force_displacement.png')
print('Created stress_strain.png')
""")

print('Created', csv_name)
print('Created', plot_script)
print('Now run: python plot_main_results.py')
