from Bio.PDB import *
import glob
import itertools
import numpy as np
import argparse

argument_parser = argparse.ArgumentParser(description = '''Search for assemble path''')
argument_parser.add_argument('--trimer_dir', type=str, help = 'Path to pdb files.')
argument_parser.add_argument('--output_dir', type=str, help = 'Path to pdb files.')
args = argument_parser.parse_args()

parser = PDBParser()

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
    if np.min(contact_map)>8:
        return False
    else:
        return True

class ChainSelect(Select):
    def __init__(self, pairs):
        self.pairs = pairs
    def accept_chain(self, chain):
        if chain.get_id() == self.pairs[0] or chain.get_id() == self.pairs[1]:
            return 1
        else:
            return 0

def write_pairs(directory):
    for pdb in glob.glob(directory):
        pdb_id = pdb[-12:-4]
        trimer_name = pdb[-7:-4]
        structure = parser.get_structure("trimer", pdb)

        chain_lst = []
        for chain in structure[0]:
            chain_lst.append(chain.get_id())
        perm = [p for p in itertools.combinations(chain_lst, 2)]        # for chain name, e.g. B, C, D
        perm_ind = [p for p in itertools.combinations(trimer_name, 2)]  # for individual chain id, e.g. 0, 1, 2, 3
        io = PDBIO()
        io.set_structure(structure)
        for i in range(len(perm)):
            chainA = perm[i][0]
            chainB = perm[i][1]
            if interacting(structure[0][chainA], structure[0][chainB]):
                io.save(f'{args.output_dir}{pdb_id}_{perm_ind[i][0]}{perm_ind[i][1]}_{chainA}{chainB}.pdb', ChainSelect(perm[i]))
            
def main():
    write_pairs(args.trimer_dir+"*")

if __name__ == '__main__':
    main()
