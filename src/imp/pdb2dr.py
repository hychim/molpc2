import sys
import os
import argparse
import itertools

import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser

parser = argparse.ArgumentParser(description = 'Converting protein interaction pdb file to distant restraint csv.')
parser.add_argument('--id', type=str, help = 'job ID')
parser.add_argument('--close_end', type=str, help = 'Path to close end txt from mcts')
parser.add_argument('--mcts_pdb', type=str, help = 'Path to final complex from mcts')
parser.add_argument('--dist', type=int, help = 'distance for setting distant restraint')
args = parser.parse_args()

pdb_id = args.id
close_end = args.close_end
dist_thres = args.dist
mcts_final = args.mcts_pdb

pdb_parser = PDBParser(PERMISSIVE=1)
structure = pdb_parser.get_structure('mcts', mcts_final)

if not os.path.exists(f'output/{pdb_id}/imp/'):
        os.makedirs(f'output/{pdb_id}/imp/')
f = open(f'output/{pdb_id}/imp/{pdb_id}_dr.csv','w')

chain_lst = []
for chain in structure[0]:
    chain_lst.append(chain.get_id())
c = list(itertools.combinations(chain_lst, 2))

for pair in c:
    chainA = structure[0][pair[0]]
    chainB = structure[0][pair[1]]

    A_id = chainA.get_id()
    B_id = chainB.get_id()

    if (A_id==close_end[0] and B_id==close_end[1]) or (A_id==close_end[1] and B_id==close_end[0]):
        None
    else:
        for residue in chainA:
            A_first = residue.get_id()[1]
            break
        for residue in chainB:
            B_first = residue.get_id()[1]
            break

        for resA in chainA:
            for resB in chainB:
                dist = resA['CA'] - resB['CA']
                if dist < dist_thres:
                    line = f"chain_{A_id},{resA.get_id()[1]-A_first+1},chain_{B_id},{resB.get_id()[1]-B_first+1},{dist} \n"
                    f.write(line)

# add close end distant restraint
if os.path.exists(close_end):
    with open(close_end) as f:
        lines = [line.rstrip() for line in f]
    close_end = str(lines[0]).split(',')
    structure = pdb_parser.get_structure('close_end', close_end[2])
    chain_lst = []
    for chain in structure[0]:
        chain_lst.append(chain)
    chainA = chain_lst[0]
    chainB = chain_lst[1]

    for residue in chainA:
        A_first = residue.get_id()[1]
        break
    for residue in chainB:
        B_first = residue.get_id()[1]
        break
    for resA in chainA:
        for resB in chainB:
            dist = resA['CA'] - resB['CA']
            if dist < dist_thres:
                line = f"chain_{close_end[0]},{resA.get_id()[1]-A_first+1},chain_{close_end[1]},{resB.get_id()[1]-B_first+1},{dist} \n"
                f.write(line)

f.close()
