import argparse
import os
import pickle

import pandas as pd
from matplotlib import pyplot as plt

from alignment import align_folder
from coordinates import extract_coords_from_pdb_by_seq, apply_pca
from plotting import plot_combined_residue_graph_pca, svg_to_html

# script_path = './get_files.sh'
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Courier New'] + plt.rcParams['font.monospace']

chain_fict = {
    'A': 'TCR_alpha',
    'B': 'TCR_beta',
    'C': 'peptide',
    'D': 'MHC'
}

amino_acid_dict = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_ca_to_dataframe(pdb_file_path: str, save_dir: str) -> pd.DataFrame:
    """
    Extract CA atom coordinates from a PDB file and save a dataframe.

    Args:
        pdb_file_path: Path to the PDB file.
        save_dir: Directory where the TSV should be written.

    Returns:
        The pandas DataFrame with CA atom coordinates.
    """
    structure_hash = os.path.basename(pdb_file_path).split('.')[0].split('aligned_')[-1]
    ca_data = []

    with open(pdb_file_path, 'r') as f:
        for line in f:
            if not line.startswith('ATOM  '):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue

            chain = line[21].strip() or '-'
            resname = line[17:20].strip()
            resnum = line[22:27].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            ca_data.append([chain, resname, resnum, round(x, 3), round(y, 3), round(z, 3)])

    df = pd.DataFrame(ca_data, columns=['Chain', 'Residue', 'ResNum', 'X', 'Y', 'Z'])
    df['Chain'] = df['Chain'].map(chain_fict).fillna(df['Chain'])
    df['Residue'] = df['Residue'].map(amino_acid_dict).fillna(df['Residue'])

    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f'{structure_hash}_aa_coordinates.tsv')
    df.to_csv(output_path, sep='\t', index=False)
    return df


data_generation = pd.read_csv('/projects/structures/clusters/HomoSapiens_MHCI_all_clusters.tsv', sep='\t')
if 'TCR_hash' not in data_generation.columns:
    raise KeyError("Expected 'TCR_hash' column in data_generation TSV")
data_generation = data_generation.set_index('TCR_hash')


def process_folder(folder: str, pca) -> None:
    """
    Process a folder of PDB files, extract coordinates, apply PCA, and generate plots.

    Args:
        folder: Path to the folder containing PDB files.
        pca: PCA object for dimensionality reduction.
    """
    for pdb_file in os.listdir(folder):
        if not pdb_file.endswith('.pdb'):
            continue

        complex_hash = pdb_file[:-4].split('_')[-1]
        metadata = data_generation.loc[complex_hash]

        tra_aa_seq = metadata['cdr3.alpha']
        trb_aa_seq = metadata['cdr3.beta']
        antigen_epitope = metadata['antigen.epitope']

        coords_ca, coords_all, atom_to_ca_map = extract_coords_from_pdb_by_seq(
            os.path.join(folder, pdb_file), tra_aa_seq, trb_aa_seq, antigen_epitope
        )

        coords_ca = apply_pca(coords_ca, pca)

        plot_combined_residue_graph_pca(
            coords_ca,
            coords_all,
            atom_to_ca_map,
            max_distance=5.0,
            pdb_filename=complex_hash,
            save_dir=os.path.join(folder, "contacts_and_skeleton_plots")
        )
        plot_combined_residue_graph_pca(
            coords_ca,
            coords_all,
            atom_to_ca_map,
            max_distance=5.0,
            pdb_filename=complex_hash,
            save_dir=os.path.join(folder, "contacts_and_skeleton_plots"),
            simplified=True
        )

        svg_to_html(
            os.path.join(folder, "contacts_and_skeleton_plots", f"{complex_hash}.svg"),
            os.path.join(folder, "contacts_and_skeleton_plots", f"{complex_hash}.html")
        )
        svg_to_html(
            os.path.join(folder, "contacts_and_skeleton_plots", f"{complex_hash}_simplified.svg"),
            os.path.join(folder, "contacts_and_skeleton_plots", f"{complex_hash}_simplified.html")
        )

        coordinates_db = pd.DataFrame(coords_ca).T
        coordinates_db.rename(columns={0: 'x', 1: 'y', 2: 'z', 3: 'PC1', 4: 'PC2'}, inplace=True)
        coordinates_db.to_csv(
            os.path.join(folder, "contacts_and_skeleton_plots", f"{complex_hash}_coordinates.tsv"),
            sep='\t'
        )

        extract_ca_to_dataframe(
            os.path.join(folder, pdb_file),
            os.path.join(folder, "contacts_and_skeleton_plots")
        )


# Load PCA model
filename = 'pca_all_structures.sav'
pca = pickle.load(open(filename, 'rb'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process PDB files for TCR structure visualization')
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        default='.',
        help='Input folder with .pdb files'
    )

    args = parser.parse_args()
    input_folder = args.input

    input_folder_align_first = input_folder.rstrip('/') + '_align'
    print(f"Aligning to: {input_folder_align_first}")

    align_folder(input_folder, input_folder_align_first)

    input_folder_align_second = input_folder.rstrip('/') + '_align_2'

    # Distribute models by antigen folders
    for structure in os.listdir(input_folder_align_first):
        if not structure.endswith('.pdb'):
            continue

        tcr_hash = structure[8:-4]
        structure_epitope = data_generation.loc[tcr_hash]['antigen.epitope']

        os.makedirs(os.path.join(input_folder_align_first, structure_epitope), exist_ok=True)
        os.rename(
            os.path.join(input_folder_align_first, structure),
            os.path.join(input_folder_align_first, structure_epitope, structure)
        )

    # Align models within antigen folders
    for epitope in os.listdir(input_folder_align_first):
        if '_' in epitope:
            continue
        align_folder(
            os.path.join(input_folder_align_first, epitope) + '/',
            os.path.join(input_folder_align_second, epitope) + '/'
        )

    # Generate plots
    for folder in os.listdir(input_folder_align_second):
        if "_" in folder:
            continue
        process_folder(os.path.join(input_folder_align_second, folder) + '/', pca)
