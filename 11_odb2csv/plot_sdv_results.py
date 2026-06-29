import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = {}

with open('sdv_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    for name in reader.fieldnames:
        data[name] = []

    for row in reader:
        for name in reader.fieldnames:
            data[name].append(float(row[name]))

strain = data['LE33']

# Plot 1: stress-strain from integration point
plt.figure()
plt.plot(strain, data['S33'], marker='o')
plt.xlabel('LE33')
plt.ylabel('S33')
plt.title('Integration Point Stress-Strain')
plt.grid(True)
plt.tight_layout()
plt.savefig('ip_stress_strain.png', dpi=300)
plt.close()

# Plot 2: total accumulated slip
plt.figure()
plt.plot(strain, data['SDV34'], marker='o')
plt.xlabel('LE33')
plt.ylabel('SDV34')
plt.title('Total Accumulated Plastic Slip')
plt.grid(True)
plt.tight_layout()
plt.savefig('total_slip_SDV34.png', dpi=300)
plt.close()

# Plot 3: 24 slip systems
plt.figure()
for i in range(10, 34):
    plt.plot(strain, data['SDV%d' % i], label='SDV%d' % i)

plt.xlabel('LE33')
plt.ylabel('Accumulated slip')
plt.title('Slip System Evolution: SDV10-SDV33')
plt.grid(True)
plt.legend(fontsize=6, ncol=3)
plt.tight_layout()
plt.savefig('slip_systems_SDV10_33.png', dpi=300)
plt.close()

# Plot 4: Euler angles
plt.figure()
plt.plot(strain, data['SDV35'], label='SDV35 phi1')
plt.plot(strain, data['SDV36'], label='SDV36 Phi')
plt.plot(strain, data['SDV37'], label='SDV37 phi2')
plt.xlabel('LE33')
plt.ylabel('Euler angles')
plt.title('Crystal Orientation Evolution')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('euler_angles_SDV35_37.png', dpi=300)
plt.close()

print('Created plots:')
print('ip_stress_strain.png')
print('total_slip_SDV34.png')
print('slip_systems_SDV10_33.png')
print('euler_angles_SDV35_37.png')
