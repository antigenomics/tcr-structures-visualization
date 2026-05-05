import os
import pymolPy3


def align_folder(folder: str, output_folder: str) -> None:
    """
    Align all PDB files in a folder to the first one using PyMOL.

    Args:
        folder (str): Path to the folder containing PDB files.
        output_folder (str): Path to the output folder for aligned PDB files.
    """
    os.makedirs(output_folder, exist_ok=True)

    pdb_files = [
        os.path.join(folder, file)
        for file in os.listdir(folder)
        if file.endswith('.pdb')
    ]

    if not pdb_files:
        print("No PDB files found in the folder.")
        return

    pm = pymolPy3.pymolPy3(0)  # Initialize in headless mode (no GUI)

    reference_file = pdb_files[0]
    ref_name = os.path.splitext(os.path.basename(reference_file))[0]
    pm(f"load {reference_file}, {ref_name}")

    aligned_ref = f"aligned_{os.path.basename(reference_file)}"
    pm(f"save {os.path.join(output_folder, aligned_ref)}, {ref_name}")
    print(f"Saved: {aligned_ref}")

    for pdb_file in pdb_files[1:]:  # Skip the reference file
        target_name = os.path.splitext(os.path.basename(pdb_file))[0]
        pm(f"load {pdb_file}, {target_name}")

        pm(f"align {target_name}, {ref_name}")

        aligned_file = f"aligned_{os.path.basename(pdb_file)}"
        pm(f"save {os.path.join(output_folder, aligned_file)}, {target_name}")
        print(f"Aligned and saved: {aligned_file}")

    pm("quit")
