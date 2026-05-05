from typing import Dict, List, Tuple, Optional
import os
import re
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from coordinates import calculate_3d_distance

AMINO_ACID_DICT = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

AMINO_ACID_DICT_REVERSE = {v: k for k, v in AMINO_ACID_DICT.items()}


def process_contacts_file(contacts_path: str) -> pd.DataFrame:
    """
    Process a contacts file and return a cleaned DataFrame.

    Args:
        contacts_path: Path to the contacts file.

    Returns:
        Processed DataFrame with amino acid contacts.
    """
    data_contacts = pd.read_csv(contacts_path, header=None)

    data_contacts[0] = data_contacts[0].apply(lambda x: x[2:-1])
    data_contacts[1] = data_contacts[1].apply(lambda x: x[2:-1])
    data_contacts[3] = data_contacts[3].apply(lambda x: x.split("') - ('"))
    data_contacts[4] = data_contacts[4].apply(lambda x: x[1:])
    data_contacts[4] = data_contacts[4].apply(lambda x: x[1:-1])
    data_contacts[6] = data_contacts[6].apply(lambda x: x[2:-2])

    data_contacts['atom_from'] = data_contacts[3].apply(lambda x: x[0][2:])
    data_contacts['chain_to'] = data_contacts[3].apply(lambda x: x[1][:-1])

    data_contacts.rename(
        columns={
            0: 'chain_from',
            1: 'aa_from',
            2: 'res_num_from',
            4: 'aa_to',
            5: 'res_num_to',
            6: 'atom_to'
        },
        inplace=True
    )
    data_contacts.drop(3, inplace=True, axis=1)

    data_contacts = data_contacts.drop_duplicates(
        ['chain_from', 'aa_from', 'chain_to', 'aa_to']
    )

    data_contacts['aa_from'] = data_contacts['aa_from'].map(AMINO_ACID_DICT)
    data_contacts['aa_to'] = data_contacts['aa_to'].map(AMINO_ACID_DICT)

    data_contacts.drop(['atom_to', 'atom_from'], axis=1, inplace=True)

    return data_contacts


def process_contacts_file(contacts_path):
    data_contacts = pd.read_csv(contacts_path, header=None)

    data_contacts[0] = data_contacts[0].apply(lambda x: x[2:-1])
    data_contacts[1] = data_contacts[1].apply(lambda x: x[2:-1])
    data_contacts[3] = data_contacts[3].apply(lambda x: x.split("') - ('"))
    data_contacts[4] = data_contacts[4].apply(lambda x: x[1:])
    data_contacts[4] = data_contacts[4].apply(lambda x: x[1:-1])
    data_contacts[6] = data_contacts[6].apply(lambda x: x[2:-2])

    data_contacts['atom_from'] = data_contacts[3].apply(lambda x: x[0][2:])
    data_contacts['chain_to'] = data_contacts[3].apply(lambda x: x[1][:-1])

    data_contacts.rename(
        columns={0: 'chain_from', 1: 'aa_from', 2: 'res_num_from', 4: 'aa_to', 5: 'res_num_to', 6: 'atom_to'},
        inplace=True)
    data_contacts.drop(3, inplace=True, axis=1)

    data_contacts = data_contacts.drop_duplicates(['chain_from', 'aa_from', 'chain_to', 'aa_to'])

    data_contacts['aa_from'] = data_contacts['aa_from'].map(AMINO_ACID_DICT)
    data_contacts['aa_to'] = data_contacts['aa_to'].map(AMINO_ACID_DICT)

    data_contacts.drop(['atom_to', 'atom_from'], axis=1, inplace=True)

    return data_contacts


def plot_combined_residue_graph_pca(
    coords_ca: Dict[Tuple[str, str, int], List[float]],
    coords_all: Dict[Tuple[str, str, int, str], List[float]],
    atom_to_ca_map: Dict[Tuple[str, str, int, str], Tuple[str, str, int]],
    max_distance: float = 5.0,
    pdb_filename: Optional[str] = None,
    save_dir: Optional[str] = None,
    simplified: bool = False
) -> None:
    """
    Plot combined residue graph in PCA space.

    Args:
        coords_ca: CA coordinates.
        coords_all: All atom coordinates.
        atom_to_ca_map: Mapping from atoms to CA.
        max_distance: Maximum distance for connections.
        pdb_filename: PDB filename for saving.
        save_dir: Directory to save files.
        simplified: If True, use simplified markers.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)

    chains = sorted(set(key[0] for key in coords_ca))
    chain_to_color = {'CDR3_alpha': 'green', 'CDR3_beta': 'red', 'peptide': 'blue'}
    legend_elements = []

    for chain_id in chains:
        chain_residues = [
            (res_key, coords_ca[res_key])
            for res_key in coords_ca
            if res_key[0] == chain_id
        ]
        chain_residues.sort(key=lambda x: x[0][2])

        legend_elements.append(Line2D(
            [0], [0],
            marker='s',
            color=chain_to_color[chain_id],
            label=chain_id,
            markersize=10,
            linestyle='-'
        ))

        pc1 = [res[1][3] for res in chain_residues]
        pc2 = [res[1][4] for res in chain_residues]

        marker = '.' if simplified else 's'
        markersize = 8 if simplified else 20

        for i, (res_key, coord) in enumerate(chain_residues):
            ax.plot(
                [coord[3]], [coord[4]],
                marker=marker,
                markersize=markersize,
                markerfacecolor=chain_to_color[chain_id],
                markeredgecolor=chain_to_color[chain_id],
                markeredgewidth=1.2,
                alpha=1,
                linestyle='none',
                zorder=3
            )

            if not simplified:
                one_letter = res_key[1]
                ax.text(
                    coord[3], coord[4],
                    f'$\\mathrm{{{one_letter}}}^{{{res_key[2]}}}$',
                    color='white',
                    fontsize=10,
                    ha='center',
                    va='center',
                    fontweight='bold',
                    zorder=10
                )

        for i in range(len(chain_residues) - 1):
            ax.plot(
                [pc1[i], pc1[i + 1]], [pc2[i], pc2[i + 1]],
                c=chain_to_color[chain_id],
                alpha=1,
                linestyle='-'
            )

    drawn_connections = set()
    ca_connections_drawn = set()
    atomic_contacts = []
    keys_all = list(coords_all.keys())

    for i, (atom_key_i, coord_i) in enumerate(coords_all.items()):
        ca_i = atom_to_ca_map.get(atom_key_i)
        if ca_i is None:
            continue
        for j in range(i + 1, len(coords_all)):
            atom_key_j = keys_all[j]
            coord_j = coords_all[atom_key_j]
            ca_j = atom_to_ca_map.get(atom_key_j)
            if ca_j is None or ca_i == ca_j:
                continue
            if ca_i[0] == ca_j[0]:
                continue
            dist = calculate_3d_distance(coord_i, coord_j)
            if dist <= max_distance:
                ca_pair = frozenset((ca_i, ca_j))
                atom_pair = frozenset((atom_key_i, atom_key_j))

                if ca_pair not in ca_connections_drawn:
                    if ca_i[0] == 'peptide' or ca_j[0] == 'peptide':
                        connection_color = 'black'
                        connection_style = 'dotted'
                        linewidth = 1.5
                    else:
                        connection_color = 'black'
                        connection_style = '--'
                        linewidth = 0.2

                    c1 = coords_ca[ca_i]
                    c2 = coords_ca[ca_j]
                    ax.plot(
                        [c1[3], c2[3]], [c1[4], c2[4]],
                        c=connection_color,
                        alpha=1,
                        linestyle=connection_style,
                        linewidth=linewidth
                    )
                    ca_connections_drawn.add(ca_pair)

                if atom_pair not in drawn_connections and atom_pair not in ca_connections_drawn:
                    drawn_connections.add(atom_pair)

                    resname_i = AMINO_ACID_DICT_REVERSE.get(ca_i[1], ca_i[1])
                    resname_j = AMINO_ACID_DICT_REVERSE.get(ca_j[1], ca_j[1])
                    atom_info_i = (atom_key_i[0], resname_i, atom_key_i[2], atom_key_i[3])
                    atom_info_j = (atom_key_j[0], resname_j, atom_key_j[2], atom_key_j[3])
                    atomic_contacts.append((atom_info_i, atom_info_j))

    ax.set_axis_off()
    plt.tight_layout()

    if pdb_filename:
        base_name = os.path.splitext(os.path.basename(pdb_filename))[0]
        suffix = '_simplified' if simplified else ''
        image_filename = f"{base_name}{suffix}.svg"
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            image_filename = os.path.join(save_dir, image_filename)
        plt.savefig(image_filename, dpi=300, bbox_inches='tight', transparent=True)

        contacts_filename = f"{base_name}_contacts.txt"
        if save_dir:
            contacts_filename = os.path.join(save_dir, contacts_filename)
        with open(contacts_filename, 'w') as f:
            for a1, a2 in atomic_contacts:
                f.write(f"{a1} - {a2}\n")

        aa_contacts = process_contacts_file(contacts_filename)
        aa_contacts_filename = f"{base_name}_aa_contacts.tsv"
        if save_dir:
            aa_contacts_filename = os.path.join(save_dir, aa_contacts_filename)
        aa_contacts.to_csv(aa_contacts_filename, sep='\t')

    plt.show()


def svg_to_html(svg_input_path: str, html_output_path: str) -> None:
    """
    Converts an SVG file to an HTML file by embedding the SVG content inline.

    Args:
        svg_input_path: Path to the input SVG file.
        html_output_path: Path to the output HTML file.
    """
    with open(svg_input_path, 'r', encoding='utf-8') as svg_file:
        svg_content = svg_file.read()

    svg_content = re.sub(r'<\?xml.*?\?\>', '', svg_content).strip()

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SVG Embedded in HTML</title>
</head>
<body>
    {svg_content}
</body>
</html>"""

    with open(html_output_path, 'w', encoding='utf-8') as html_file:
        html_file.write(html_template)

    print(f"SVG converted to HTML: {html_output_path}")
    """
    Converts an SVG file to an HTML file by embedding the SVG content inline.
    
    Args:
    svg_input_path (str): Path to the input SVG file.
    html_output_path (str): Path to the output HTML file.
    """
    with open(svg_input_path, 'r', encoding='utf-8') as svg_file:
        svg_content = svg_file.read()

    svg_content = re.sub(r'<\?xml.*?\?\>', '', svg_content).strip()

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SVG Embedded in HTML</title>
</head>
<body>
    {svg_content}
</body>
</html>"""

    with open(html_output_path, 'w', encoding='utf-8') as html_file:
        html_file.write(html_template)

    print(f"SVG converted to HTML: {html_output_path}")
