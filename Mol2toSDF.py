import re
from rdkit import Chem
from openbabel import openbabel

## change input_file to your mol2 file with all your ligands docked
input_file = "resultsWithHbonds.mol2"
output_file = "output.sdf"

# Step 1: Split the mol2 file into individual molecule blocks
with open(input_file, 'r') as f:
    content = f.read()

# Split at molecule header
molecules = content.split("@<TRIPOS>MOLECULE")

writer = Chem.SDWriter(output_file)

for block in molecules[1:]:
    # Extract properties from ########## lines
    props = dict(re.findall(r"^#{10}\s*(.+?):\s*(.+)$", block, re.MULTILINE))

    # Recreate mol2 text (must start with @<TRIPOS>MOLECULE)
    mol2_text = "@<TRIPOS>MOLECULE" + block

    # Convert mol2 to RDKit molecule using OpenBabel backend
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "sdf")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, mol2_text)

    # Convert to RDKit molecule object
    molblock = obConversion.WriteString(mol)
    rdmol = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
    if rdmol is None:
        continue

    # Add the extracted properties as SD tags
    for key, value in props.items():
        rdmol.SetProp(key.strip(), value.strip())

    writer.write(rdmol)

writer.close()
print(f"✅ Done! Saved to {output_file}")
