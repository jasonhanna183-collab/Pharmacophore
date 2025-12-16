import os
import zipfile
from openbabel import openbabel, pybel

def sdf_to_pdb_zip(input_sdf, output_zip="converted_pdbs.zip"):
    # Load all molecules from SDF
    mols = list(pybel.readfile("sdf", input_sdf))

    # Make output folder
    outdir = "converted_pdbs"
    os.makedirs(outdir, exist_ok=True)

    pdb_files = []

    for i, mol in enumerate(mols, start=1):
        # Create a safe base filename
        name = mol.title if mol.title.strip() else f"mol_{i}"
        safe_name = "".join(c if c.isalnum() or c in "._-" else "_" for c in name)

        outfile = os.path.join(outdir, f"{safe_name}.pdb")

        # Write PDB equivalent to: obabel input.sdf -O name.pdb -m
        mol.write("pdb", outfile, overwrite=True)
        pdb_files.append(outfile)

    # Create ZIP file containing all PDBs
    with zipfile.ZipFile(output_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
        for pdb in pdb_files:
            zipf.write(pdb, arcname=os.path.basename(pdb))

    print(f"✔ Converted {len(pdb_files)} molecules")
    print(f"✔ PDB files saved to: {outdir}/")
    print(f"✔ ZIP archive created: {output_zip}")


# Example usage:
sdf_to_pdb_zip("combinedpharmithitsfiltered.sdf")
