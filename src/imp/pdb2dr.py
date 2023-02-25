import sys
import os
import argparse

import Bio.PDB
import Bio.PDB.StructureBuilder
from Bio.PDB.Residue import Residue
from Bio.PDB import PDBParser

### setting
pdb_parser = PDBParser(PERMISSIVE=1)

parser = argparse.ArgumentParser(description = 'Converting protein interaction pdb file to distant restraint csv.')
parser.add_argument('--dimer_dir', type=str, help = 'Where is the dimer dir')
parser.add_argument('--dimer_lst', type=str, help = 'dimer list')
parser.add_argument('--dist', type=int, help = 'distance for setting distant restraint')
parser.add_argument('--outdir', type=str, help = 'Where to write the result')
args = parser.parse_args()

dimer_dir = args.dimer_dir
dimer_lst = args.dimer_lst
dist = args.dist
outdir = args.outdir

with open(dimer_lst) as f:
    lines = f.read().splitlines()

# loop through dimer prediction pdb
directory = os.fsencode(dimer_dir)

# Read structure from files
f = open(outdir,'w')



for pdb in lines[3:]:
    line = pdb.split(',')
    print(line)
    structure_id = line[2]
    filename = f"{dimer_dir}{structure_id}"
    structure = pdb_parser.get_structure(structure_id, filename)
    model = structure[0]
    a,b = line[0], line[1]

    chain_lst = []
    for chain in model.get_chains():
        chain_lst.append(chain.get_id())
    chainA = model[chain_lst[0]]
    chainB = model[chain_lst[1]]


    for residue in chainA:
        A_first = residue.get_id()[1]
        break
    for residue in chainB:
        B_first = residue.get_id()[1]
        break

    # loop through all residue in chain
    # ONLY WORKS FOR DIMER NOW
    for residue1 in chainA:
        for residue2 in chainB:
            try:
                # compute distance between CA atoms
                # should use NeighborSearch modeule in biopython
                distance = residue1['CA'] - residue2['CA'] 
            except KeyError:
                ## no CA atom, e.g. for H_NAG
                continue
            if distance < dist:
                line_out = f"chain_{a},{residue1.get_id()[1]-A_first+1},chain_{b},{residue2.get_id()[1]-B_first+1},{distance} \n"
                f.write(line_out)

f.close()
