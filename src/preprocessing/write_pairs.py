from Bio.PDB import *
import glob
import numpy as np
import os

import argparse

argument_parser = argparse.ArgumentParser(description = '''Search for assemble path''')
argument_parser.add_argument('--trimer_dir', type=str, help = 'Path to pdb files.')
argument_parser.add_argument('--output_dir', type=str, help = 'Path to pdb files.')
args = argument_parser.parse_args()

parser = PDBParser()

class ChainSelectBC(Select):
    def accept_chain(self, chain):
        if chain.get_id() == "B" or chain.get_id() == "C":
            return 1
        else:
            return 0
class ChainSelectBD(Select):
    def accept_chain(self, chain):
        if chain.get_id() == "B" or chain.get_id() == "D":
            return 1
        else:
            return 0
class ChainSelectCD(Select):
    def accept_chain(self, chain):
        if chain.get_id() == "C" or chain.get_id() == "D":
            return 1
        else:
            return 0

def contact_map_btw(chainA, chainB):
    contact_map = np.zeros((1,len(chainB)))
    
    for r_A in chainA:
        dists = []
        for r_B in chainB:
            dists.append(r_A["CA"] - r_B["CA"])
        contact_map = np.vstack([contact_map, dists])
    return np.delete(contact_map, 0, 0)

def interacting(chainA, chainB):
    contact_map = contact_map_btw(chainA, chainB)
    if np.min(contact_map)>6:
        return False
    else:
        return True

def write_pairs(directory):
    for pdb in glob.glob(directory):
        pdb_id = pdb[-12:-4]
        trimer_name = pdb[-7:-4]
        structure = parser.get_structure("trimer", pdb)

        io = PDBIO()
        io.set_structure(structure)
        if interacting(structure[0]["B"], structure[0]["C"]):
            io.save(args.output_dir+pdb_id+"_"+trimer_name[0]+trimer_name[1]+"_BC.pdb", ChainSelectBC())
        if interacting(structure[0]["B"], structure[0]["D"]):
            io.save(args.output_dir+pdb_id+"_"+trimer_name[0]+trimer_name[2]+"_BD.pdb", ChainSelectBD())
        if interacting(structure[0]["C"], structure[0]["D"]):
            io.save(args.output_dir+pdb_id+"_"+trimer_name[1]+trimer_name[2]+"_CD.pdb", ChainSelectCD())

def main():
    os.makedirs(args.output_dir) if not os.path.exists(args.output_dir) else None
    write_pairs(args.trimer_dir+"*")

if __name__ == '__main__':
    main()