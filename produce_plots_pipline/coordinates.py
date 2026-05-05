from typing import Dict, List, Tuple, Any
from Bio import PDB
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


def extract_coords_from_pdb_by_seq(
    pdb_file: str,
    cdr3a: str,
    cdr3b: str,
    peptide: str
) -> Tuple[Dict[Tuple[str, str, int], List[float]],
           Dict[Tuple[str, str, int, str], List[float]],
           Dict[Tuple[str, str, int, str], Tuple[str, str, int]]]:
    """
    Extract coordinates from a PDB file based on sequences.

    Args:
        pdb_file: Path to the PDB file.
        cdr3a: CDR3 alpha sequence.
        cdr3b: CDR3 beta sequence.
        peptide: Peptide sequence.

    Returns:
        A tuple of (coords_ca, coords_all, atom_to_ca_map).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    coords_ca = {}
    coords_all = {}
    atom_to_ca_map = {}

    for model in structure:
        for chain in model.get_chains():
            residues = [res for res in chain if PDB.Polypeptide.is_aa(res)]
            pdb_sequence = "".join(seq1(res.get_resname()) for res in residues)

            if cdr3a in pdb_sequence:
                curr_seq = 'CDR3_alpha'
                start_idx = pdb_sequence.find(cdr3a)
                sequence = cdr3a
            elif cdr3b in pdb_sequence:
                curr_seq = 'CDR3_beta'
                start_idx = pdb_sequence.find(cdr3b)
                sequence = cdr3b
            elif peptide in pdb_sequence:
                curr_seq = 'peptide'
                start_idx = pdb_sequence.find(peptide)
                sequence = peptide
            else:
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
                    ca_key = None

                for atom in res.get_atoms():
                    atom_name = atom.get_name()
                    atom_key = (curr_seq, res_name, res_num, atom_name)
                    coords_all[atom_key] = atom.coord.tolist()
                    if ca_key:
                        atom_to_ca_map[atom_key] = ca_key

    return coords_ca, coords_all, atom_to_ca_map


def apply_pca(data: Dict[Any, List[float]], pca) -> Dict[Any, List[float]]:
    """
    Apply PCA transformation to the data.

    Args:
        data: Dictionary of coordinates.
        pca: PCA object.

    Returns:
        Updated data with PCA components.
    """
    keys = list(data.keys())
    points = np.array([data[key] for key in keys])

    projected = pca.transform(points)

    for i, key in enumerate(keys):
        data[key].extend(projected[i].tolist())

    return data


def calculate_3d_distance(coord1: List[float], coord2: List[float]) -> float:
    """
    Calculate the 3D Euclidean distance between two coordinates.

    Args:
        coord1: First coordinate.
        coord2: Second coordinate.

    Returns:
        Distance between the coordinates.
    """
    return np.linalg.norm(np.array(coord1) - np.array(coord2))
