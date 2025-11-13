from turtle import pd

import matplotplib.pyplot as plt
import subprocess
import pickle

from alingment import *
from coordinates import *
from plotting import *

script_path = './get_files.sh'
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Courier New'] + plt.rcParams['font.monospace']

data_generation = pd.read_csv('/projects/structures/clusters/HomoSapiens_MHCI_all_clusters.tsv', sep='\t', index_col=0)


def process_folder(folder, pca):
    for pdb_file in os.listdir(folder):
        if pdb_file.endswith('.pdb'):
            complex_hash = pdb_file.split('/')[-1][:-4].split('_')[-1]

            metadata_for_generation_i = data_generation[data_generation['TCR_hash'] == complex_hash].iloc[0]

            TRA_aa_seq = metadata_for_generation_i['cdr3.alpha']
            TRB_aa_seq = metadata_for_generation_i['cdr3.beta']
            antigen_epitope = metadata_for_generation_i['antigen.epitope']

            # NLVPMVATV
            # coords_ca, coords_all, atom_to_ca_map = extract_coords_from_pdb(pdb_file,  chains=['C', 'A', 'B'], sequences=['NLVPMVATV', 'CIRDNNNDMRF', 'CASSLAPGATNEKLFF'])
            coords_ca, coords_all, atom_to_ca_map = extract_coords_from_pdb_by_seq(folder + pdb_file, TRA_aa_seq,
                                                                                   TRB_aa_seq, antigen_epitope)

            coords_ca = apply_pca(coords_ca, pca)

            plot_combined_residue_graph_pca(coords_ca,
                                            coords_all,
                                            atom_to_ca_map,
                                            max_distance=5.0,
                                            pdb_filename=complex_hash,
                                            save_dir=f"{folder}/contacts_and_skeleton_plots/")
            plot_combined_residue_graph_pca_simple(coords_ca,
                                                   coords_all,
                                                   atom_to_ca_map,
                                                   max_distance=5.0,
                                                   pdb_filename=complex_hash,
                                                   save_dir=f"{folder}/contacts_and_skeleton_plots/")

            svg_to_html(f'{folder}/contacts_and_skeleton_plots/{complex_hash}.svg',
                        f'{folder}/contacts_and_skeleton_plots/{complex_hash}.html')

            svg_to_html(f'{folder}/contacts_and_skeleton_plots/{complex_hash}_simplified.svg',
                        f'{folder}/contacts_and_skeleton_plots/{complex_hash}_simplified.html')

            coordinates_db = pd.DataFrame(coords_ca).T
            coordinates_db.rename(columns={0: 'x', 1: 'y', 2: 'z', 3: 'PC1', 4: 'PC2'}, inplace=True)
            coordinates_db.to_csv(f'{folder}/contacts_and_skeleton_plots/{complex_hash}_coordinates.tsv', sep='\t')
        else:
            print(pdb_file)
            continue


filename = 'pca_all_structures.sav'
pca = pickle.load(open(filename, 'rb'))

if __name__ == '__main__':

    # get models and aling them all
    args = ['/projects/structures/clusters/HLA-A/', '/projects/structures/clusters/MHCI_data/']
    subprocess.run([script_path] + args)

    args = ['/projects/structures/clusters/HLA-B/', '/projects/structures/clusters/MHCI_data/']
    subprocess.run([script_path] + args)

    align_folder('/projects/structures/clusters/MHCI_data/', '/projects/structures/clusters/MHCI_to_align')

    # destribute models by antigen folders
    for structure in os.listdir('/projects/structures/clusters/MHCI_to_align/'):

        if structure.endswith('.pdb'):
            tcr_hash = structure[8:-4]

            structure_epitope = data_generation.loc[tcr_hash]['antigen.epitope']

            os.makedirs(f'/projects/structures/clusters/MHCI_to_align/{structure_epitope}', exist_ok=True)
            os.rename(f'/projects/structures/clusters/MHCI_to_align/{structure}',
                      f'/projects/structures/clusters/MHCI_to_align/{structure_epitope}/{structure}')

    # aling models within antigen folders
    for epitope in os.listdir('/projects/structures/clusters/MHCI_to_align/'):
        if '_' not in epitope:
            align_folder(f'/projects/structures/clusters/MHCI_to_align/{epitope}/',
                         f'/projects/structures/clusters/MHCI_align/{epitope}/')

    # generate plots
    for folder in os.listdir('/projects/structures/clusters/MHCI_align/'):
        if "_" not in folder:
            process_folder(f'/projects/structures/clusters/MHCI_align/{folder}/', pca)
