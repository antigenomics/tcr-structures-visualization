import pymolPy3
import os


def align_folder(folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)

    pdb_files = []
    for file in os.listdir(folder):
        if file.endswith('.pdb'):
            pdb_files.append(f'{folder}/{file}')

    pm = pymolPy3.pymolPy3(0)  # Initialize in headless mode (no GUI)

    reference_file = pdb_files[0]

    ref_name = reference_file.split('/')[-1][:-4]
    pm(f"load {reference_file}, {ref_name}")

    aligned_ref = 'aligned_' + pdb_files[0].split('/')[-1]
    pm(f"save {output_folder}/{aligned_ref}, {ref_name}")
    print(f"Saved: {aligned_ref}")

    for pdb_file in pdb_files:
        target_name = pdb_file.split('/')[-1][:-4]
        pm(f"load {pdb_file}, {target_name}")

        pm(f"align {target_name}, {ref_name}")

        aligned_file = 'aligned_' + target_name + '.pdb'
        pm(f"save {output_folder}/{aligned_file}, {target_name}")
        print(f"Aligned and saved: {aligned_file}")

    pm("quit")
