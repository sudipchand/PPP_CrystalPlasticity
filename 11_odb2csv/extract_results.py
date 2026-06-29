from odbAccess import *

odb = openOdb(path='cp_uniaxial_1.odb')

print("\nAVAILABLE STEPS:")
for step_name in odb.steps.keys():
    print(step_name)

step = odb.steps.values()[0]

print("\nAVAILABLE HISTORY REGIONS:")
for key in step.historyRegions.keys():
    print(key)

print("\nAVAILABLE FIELD OUTPUTS:")

lastFrame = step.frames[-1]

for field in lastFrame.fieldOutputs.keys():
    print(field)

odb.close()
