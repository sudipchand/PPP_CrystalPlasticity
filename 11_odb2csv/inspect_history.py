from odbAccess import *

odb = openOdb(path='cp_uniaxial_1.odb')
step = odb.steps['TENSION']

for region_name, region in step.historyRegions.items():
    print("\nREGION:", region_name)
    for output_name in region.historyOutputs.keys():
        print("   ", output_name)

odb.close()
