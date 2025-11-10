import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Bio.PDB import PDBParser
import numpy as np
from Bio.Data.IUPACData import protein_letters_1to3
from matplotlib.lines import Line2D
from Bio.Data.IUPACData import protein_letters_3to1
from Bio import PDB


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_1to3
from itertools import combinations
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import Polypeptide
import re
import os
from Bio.SeqUtils import seq1

import pandas as pd

from sklearn.decomposition import PCA


def extract_coords_from_pdb_by_seq(pdb_file, CDR3a, CDR3b, peptide):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    
    coords_ca = {}       
    coords_all = {}
    atom_to_ca_map = {}  

    for model in structure:
        for chain in model.get_chains():
            residues = [res for res in chain if PDB.Polypeptide.is_aa(res)]

            pdb_sequence = "".join(seq1(res.get_resname()) for res in residues)
            
            if CDR3a in pdb_sequence:
                curr_seq = 'CDR3_alpha'
                start_idx = pdb_sequence.find(CDR3a)
                sequence = CDR3a
            elif CDR3b in pdb_sequence:
                curr_seq = 'CDR3_beta'
                start_idx = pdb_sequence.find(CDR3b)
                sequence = CDR3b
            elif peptide in pdb_sequence:
                curr_seq = 'peptide'
                start_idx = pdb_sequence.find(peptide)
                sequence = peptide
            else:
                start_idx = -1
                continue
                
            if start_idx == -1:
                print(f"Warning: sequence not found in chain {curr_seq}")
                continue

            for i, res in enumerate(residues[start_idx:start_idx + len(sequence)]):
                res_name = seq1(res.get_resname())
                res_num = res.id[1]

               
                if 'CA' in res:
                    ca_coord = res['CA'].coord.tolist()
                    ca_key = (curr_seq, res_name, res_num)
                    coords_ca[ca_key] = ca_coord
                else:
                    ca_coord = None
                    ca_key = None

                
                for atom in res.get_atoms():
                    atom_name = atom.get_name()
                    atom_key = (curr_seq, res_name, res_num, atom_name)
                    coords_all[atom_key] = atom.coord.tolist()
                    if ca_key:
                        atom_to_ca_map[atom_key] = ca_key

    return coords_ca, coords_all, atom_to_ca_map

def apply_pca(data, pca):
    
    keys = list(data.keys())
    points = np.array([data[key] for key in keys])

    projected = pca.transform(points)

    for i, key in enumerate(keys):
        data[key].extend(projected[i].tolist())
    
    return data

def calculate_3d_distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))