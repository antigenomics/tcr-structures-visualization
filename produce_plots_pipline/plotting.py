from coordinates import *
import re


def plot_combined_residue_graph_pca(coords_ca, coords_all, atom_to_ca_map, max_distance=5.0,
                                 pdb_filename=None, save_dir=None):
 
    amino_acid_dict = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111,)
    
    
    chains = sorted(set(key[0] for key in coords_ca))
    
    chain_to_color = {'CDR3_alpha': 'green', 'CDR3_beta': 'red', 'peptide': 'blue'}
    legend_elements = []
    j= 0
    for chain_id in chains:
        chain_residues = [(res_key, coords_ca[res_key]) for res_key in coords_ca if res_key[0] == chain_id]
        chain_residues.sort(key=lambda x: x[0][2])

        legend_elements.append(Line2D([0], [0],
                                      marker='s',
                                      color=chain_to_color[chain_id],
                                      label=chain_id,
                                      markersize=10,
                                      linestyle='-'))

        
        pc1 = [res[1][3] for res in chain_residues]
        pc2 = [res[1][4] for res in chain_residues]
        
        
        for i, (res_key, coord) in enumerate(chain_residues):
            ax.plot([coord[3]], [coord[4]],
                    marker='s',
                    markersize=20,
                    markerfacecolor=chain_to_color[chain_id],
                    markeredgecolor=chain_to_color[chain_id],
                    markeredgewidth=1.2,
                    alpha=1,
                    linestyle='none',
                    zorder=3)

            one_letter = res_key[1]
            ax.text(coord[3], coord[4],
                    f'$\\mathrm{{{one_letter}}}^{{{res_key[2]}}}$',
                    color='white',
                    fontsize=10,
                    ha='center', va='center',
                    fontweight='bold',
                    zorder=10)

        for i in range(len(chain_residues) - 1):
            ax.plot([pc1[i], pc1[i + 1]], [pc2[i], pc2[i + 1]],
                    c=chain_to_color[chain_id], alpha=1, linestyle='-')
        j+=1
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
                        connction_style = 'dotted'
                        linewidth = 1.5
                    else:
                        connection_color = 'black'
                        connction_style = '--'
                        linewidth = 0.2
                        
                        
                    c1 = coords_ca[ca_i]
                    c2 = coords_ca[ca_j]
                    ax.plot([c1[3], c2[3]], [c1[4], c2[4]],
                            c=connection_color, alpha=1, linestyle=connction_style, linewidth=linewidth)
                    ca_connections_drawn.add(ca_pair)

              
                if atom_pair not in drawn_connections and atom_pair not in ca_connections_drawn:
                    drawn_connections.add(atom_pair)

                    resname_i = amino_acid_dict.get(ca_i[1], ca_i[1])
                    resname_j = amino_acid_dict.get(ca_j[1], ca_j[1])
                    atom_info_i = (atom_key_i[0], resname_i, atom_key_i[2], atom_key_i[3])
                    atom_info_j = (atom_key_j[0], resname_j, atom_key_j[2], atom_key_j[3])
                    atomic_contacts.append((atom_info_i, atom_info_j))

    num_drawn_connections = len(ca_connections_drawn)         
    num_nonredundant_contacts = len(drawn_connections)        
    
    ax.set_axis_off()
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')

    #ax.legend(handles=legend_elements, loc='best')
    plt.tight_layout()

    if pdb_filename:
        base_name = os.path.splitext(os.path.basename(pdb_filename))[0]
        image_filename = base_name + '.svg'
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            image_filename = os.path.join(save_dir, image_filename)
        plt.savefig(image_filename, dpi=300, bbox_inches='tight', transparent=True)
        print(f'Сохранено изображение: {image_filename}')

        contacts_filename = base_name + '_contacts.txt'
        if save_dir:
            contacts_filename = os.path.join(save_dir, contacts_filename)
        with open(contacts_filename, 'w') as f:
            for a1, a2 in atomic_contacts:
                f.write(f"{a1} - {a2}\n")
        print(f'Сохранён список атомных контактов: {contacts_filename}')

    plt.show()

    return atomic_contacts, num_drawn_connections, num_nonredundant_contacts


def plot_combined_residue_graph_pca_simple(coords_ca, coords_all, atom_to_ca_map, max_distance=5.0,
                                 pdb_filename=None, save_dir=None):
 
    amino_acid_dict = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111,)
    
    
    chains = sorted(set(key[0] for key in coords_ca))
    
    chain_to_color = {'CDR3_alpha': 'green', 'CDR3_beta': 'red', 'peptide': 'blue'}
    legend_elements = []
    j= 0
    for chain_id in chains:
        chain_residues = [(res_key, coords_ca[res_key]) for res_key in coords_ca if res_key[0] == chain_id]
        chain_residues.sort(key=lambda x: x[0][2])

        legend_elements.append(Line2D([0], [0],
                                      marker='s',
                                      color=chain_to_color[chain_id],
                                      label=chain_id,
                                      markersize=10,
                                      linestyle='-'))

        
        pc1 = [res[1][3] for res in chain_residues]
        pc2 = [res[1][4] for res in chain_residues]
        
        
        for i, (res_key, coord) in enumerate(chain_residues):
            ax.plot([coord[3]], [coord[4]],
                    marker='.',
                    markersize=8,
                    markerfacecolor=chain_to_color[chain_id],
                    markeredgecolor=chain_to_color[chain_id],
                    markeredgewidth=1.2,
                    alpha=1,
                    linestyle='none',
                    zorder=3)

        for i in range(len(chain_residues) - 1):
            ax.plot([pc1[i], pc1[i + 1]], [pc2[i], pc2[i + 1]],
                    c=chain_to_color[chain_id], alpha=1, linestyle='-')
        j+=1
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
                        connction_style = 'dotted'
                        linewidth = 1.5
                    else:
                        connection_color = 'black'
                        connction_style = '--'
                        linewidth = 0.2
                        
                        
                    c1 = coords_ca[ca_i]
                    c2 = coords_ca[ca_j]
                    ax.plot([c1[3], c2[3]], [c1[4], c2[4]],
                            c=connection_color, alpha=1, linestyle=connction_style, linewidth=linewidth)
                    ca_connections_drawn.add(ca_pair)

              
                if atom_pair not in drawn_connections and atom_pair not in ca_connections_drawn:
                    drawn_connections.add(atom_pair)

                    resname_i = amino_acid_dict.get(ca_i[1], ca_i[1])
                    resname_j = amino_acid_dict.get(ca_j[1], ca_j[1])
                    atom_info_i = (atom_key_i[0], resname_i, atom_key_i[2], atom_key_i[3])
                    atom_info_j = (atom_key_j[0], resname_j, atom_key_j[2], atom_key_j[3])
                    atomic_contacts.append((atom_info_i, atom_info_j))
        
    
    ax.set_axis_off()

    plt.tight_layout()

    if pdb_filename:
        base_name = os.path.splitext(os.path.basename(pdb_filename))[0]
        image_filename = base_name +'_simplified'+ '.svg'
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            image_filename = os.path.join(save_dir, image_filename)
        plt.savefig(image_filename, dpi=300, bbox_inches='tight', transparent=True)

    plt.show()


def svg_to_html(svg_input_path, html_output_path):
    """
    Converts an SVG file to an HTML file by embedding the SVG content inline.
    
    Args:
    svg_input_path (str): Path to the input SVG file.
    html_output_path (str): Path to the output HTML file.
    """
    # Read the SVG content
    with open(svg_input_path, 'r', encoding='utf-8') as svg_file:
        svg_content = svg_file.read()

    # Remove XML declaration if present (e.g., <?xml ... ?>)
    svg_content = re.sub(r'<\?xml.*?\?\>', '', svg_content).strip()

    # Create the HTML template with the SVG embedded
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

    # Write to the output HTML file
    with open(html_output_path, 'w', encoding='utf-8') as html_file:
        html_file.write(html_template)

    print(f"SVG converted to HTML: {html_output_path}")