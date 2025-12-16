from rdkit import Chem
import re
import os


def convert_multimol2_to_sdf(mol2_file_path, sdf_file_path):
    """
    Converts a multi-molecule MOL2 file (with properties *before* the
    @<TRIPOS>MOLECULE delimiter) to a multi-molecule SDF file.

    *** v6: Adds SMILES string as a property to the SDF. ***
    """

    if not os.path.exists(mol2_file_path):
        print(f"Error: Input file not found at {mol2_file_path}")
        return

    print(f"Starting conversion from {mol2_file_path} to {sdf_file_path}...")

    try:
        with open(mol2_file_path, 'r') as f:
            mol2_content = f.read()
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Split by the unique property header
    mol_chunks = mol2_content.split('########## Name:')

    if len(mol_chunks) < 2:
        print(f"Error: Could not find '########## Name:' delimiters. Is this the correct file?")
        return

    print(f"Found {len(mol_chunks) - 1} potential molecule blocks.")

    # Open an SDWriter for the output file
    try:
        writer = Chem.SDWriter(sdf_file_path)
    except Exception as e:
        print(f"Error opening output file for writing: {e}")
        return

    mol_count = 0
    fail_count = 0

    # Iterate through the blocks, *skipping the first one (index 0)*
    for chunk in mol_chunks[1:]:

        mol = None  # Ensure mol is defined in this scope
        try:
            # --- 1. Separate Properties from MOL2 Data ---
            mol2_delimiter = '@<TRIPOS>MOLECULE'
            delimiter_pos = chunk.find(mol2_delimiter)

            if delimiter_pos == -1:
                print("Warning: Found '########## Name:' but no '@<TRIPOS>MOLECULE'. Skipping.")
                fail_count += 1
                continue

            property_text = "Name:" + chunk[:delimiter_pos]
            mol2_block = chunk[delimiter_pos:]

            # --- 2. Load the MOL2 Block ---
            mol = Chem.MolFromMol2Block(mol2_block, sanitize=False, removeHs=False)

            if mol is None:
                print("Warning: Failed to parse MOL2 data part. Skipping.")
                fail_count += 1
                continue

            # --- 3. Parse and Add Properties ---
            prop_matches = re.findall(r'#*\s*([^:\n]+):\s*(.*)', property_text)

            for key, value in prop_matches:
                key = key.strip().replace(" ", "_")
                value = value.strip()

                if key == "Name":
                    mol.SetProp("_Name", value)
                else:
                    mol.SetProp(key, value)

            # --- 4. NEW: Calculate SMILES and Add as Property ---
            # This line will generate the canonical SMILES.
            # It may fail and trigger the 'except' block for the 22 bad molecules.
            smiles_string = Chem.MolToSmiles(mol, isomericSmiles=True)
            mol.SetProp("SMILES", smiles_string)

            # --- 5. Write to SDF (Main attempt) ---
            mol.UpdatePropertyCache(strict=False)
            writer.write(mol)
            mol_count += 1

        except Exception as e:
            # --- 6. FALLBACK for "bad" molecules ---
            if "kekulize" in str(e) and mol is not None:
                print(f"Warning: Kekulization failed. Writing de-aromatized.")
                try:
                    # Manually clear all aromatic flags
                    for atom in mol.GetAtoms():
                        atom.SetIsAromatic(False)
                    for bond in mol.GetBonds():
                        bond.SetIsAromatic(False)

                    # --- NEW: Add fallback SMILES ---
                    # Now that it's de-aromatized, this will generate a
                    # non-aromatic SMILES string without failing.
                    fallback_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    mol.SetProp("SMILES", fallback_smiles)

                    # Now try to write the molecule again
                    writer.write(mol)
                    mol_count += 1  # It's a success!
                except Exception as e2:
                    print(f"Error: Still could not write molecule after de-aromatizing: {e2}")
                    fail_count += 1
            else:
                print(f"An unexpected error occurred: {e}")
                fail_count += 1

    # Close the writer
    writer.close()

    print("---")
    print("✅ Conversion Complete.")
    print(f"Successfully processed and wrote: {mol_count} molecules.")
    print(f"Failed to parse: {fail_count} blocks.")


# --- How to use the function ---

# 1. Define your file paths
input_mol2_file = "site2r5.mol2"
output_sdf_file = "site2r5.sdf"

# 2. Run the conversion
# Make sure to uncomment the line below and change the paths:
convert_multimol2_to_sdf(input_mol2_file, output_sdf_file)
