import argparse

argument_parser = argparse.ArgumentParser(description = '''Build topology files for IMP.''')
argument_parser.add_argument('--id', type=str, help = 'pdb id')
argument_parser.add_argument('--input', type=str, help = 'input')
args = argument_parser.parse_args()

# init
title = [
    "molecule_name", 
    "color", 
    "fasta_fn", 
    "fasta_id", 
    "pdb_fn", 
    "chain", 
    "residue_range", 
    "pdb_offset", 
    "bead_size", 
    "em_residues_per_gaussian", 
    "rigid_body", 
    "super_rigid_body", 
    "chain_of_super_rigid_bodies"]
chain_name = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
top60_color = [
    "tan",
    "sienna",
    "brown",
    "dark red",
    "firebrick",
    "salmon",
    "red",
    "coral",
    "sandy brown",
    "orange red",
    "orange",
    "goldenrod",
    "gold",
    "yellow",
    "khaki",
    "dark khaki",
    "dark olive green",
    "olive drab",
    "chartreuse",
    "green",
    "dark green",
    "forest green",
    "lime green",
    "light green",
    "sea green",
    "spring green",
    "dark cyan",
    "light sea green",
    "turquoise",
    "aquamarine",
    "cyan",
    "deep sky blue",
    "dodger blue",
    "steel blue",
    "sky blue",
    "light blue",
    "blue",
    "medium blue",
    "cornflower blue",
    "navy blue",
    "dark slate blue",
    "medium purple",
    "purple",
    "plum",
    "orchid",
    "magenta",
    "dark magenta",
    "violet red",
    "hot pink",
    "pink",
    "deep pink",
    "rosy brown",
    "slate gray",
    "dark slate gray",
    "white",
    "light gray",
    "gray",
    "dark gray",
    "dim gray",
    "black"]


with open(args.input) as f:
    lines = f.read().splitlines()

pdb_id = args.id

pdb_path = lines[0]
chains = lines[1].split(' ')
chains_id = lines[2].split(' ')

length = len(chains)

output = []

output.append("|molecule_name|color|fasta_fn|fasta_id |pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|")
for idx in range(length):
    output.append(f'|chain_{chains[idx]}|{top60_color[idx]:17}|{pdb_id}.fasta|{pdb_id}_{chains_id[idx]}|{pdb_path}|{chains[idx]}|1,END|0|20|0|{idx+1:2}||')

with open(f'output/{args.id}/imp/{args.id}_topology.txt', 'w') as f:
    for line in output:
        f.write(f"{line}\n")
