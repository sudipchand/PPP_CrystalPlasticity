from odbAccess import *

odb = openOdb('cp_uniaxial_1.odb')

step = odb.steps['TENSION']
lastFrame = step.frames[-1]

for name, field in lastFrame.fieldOutputs.items():
    print("\nFIELD:", name)
    print("Description:", field.description)
    print("Number of values:", len(field.values))

odb.close()
