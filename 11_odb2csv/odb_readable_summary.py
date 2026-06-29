from odbAccess import *

odb_name = 'cp_uniaxial_1.odb'
out_name = 'odb_original_structure.txt'

odb = openOdb(path=odb_name, readOnly=True)

with open(out_name, 'w') as f:

    f.write(str(odb))
    f.write("\n\nODB KEYS:\n")
    f.write(str(odb.__dict__))

odb.close()

print("Created:", out_name)
