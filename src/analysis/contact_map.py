import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import argparse

argument_parser = argparse.ArgumentParser(description = '''Plotting protein contact map from pdb files''')
argument_parser.add_argument('--pdb', type=str, help = 'Path to pdb files.')
argument_parser.add_argument('--output', type=str, help = 'Path for output.')
args = argument_parser.parse_args()

def cal_dist(r1, r2):
        # r1 means residue 1, r2 means residue 2
        return r1["CA"] - r2["CA"]
    
def clean_pdb(structure):
    new_structure = []

    for chain in structure[0]:
        new_chain = []
        for r in chain:
            if r.id[0] == " ":
                new_chain.append(r)
        new_structure.append(new_chain)
    return new_structure
    
def contact_map_btw(chainA, chainB):
    contact_map = np.zeros((1,len(chainB)))

    for r_A in chainA:
        dists = []
        for r_B in chainB:
            dists.append(r_A["CA"] - r_B["CA"])
        contact_map = np.vstack([contact_map, dists])
    return np.delete(contact_map, 0, 0)

def contact_map_total(structure):
    #borderline = 10
    total_res = 0
    num_chain = 0
    
    structure = clean_pdb(structure)
    
    for chain in structure:
        total_res += len(chain)
        num_chain += 1
    
    contact_map = np.zeros((1, total_res+num_chain*10+1))
    for chainA in structure:
        contact_map_row = np.zeros((len(chainA), 1))
        for chainB in structure:
            contact_map_sub = contact_map_btw(chainA, chainB)
            contact_map_row = np.hstack([contact_map_row, contact_map_sub])
            contact_map_row = np.hstack([contact_map_row, np.zeros((len(chainA), 10))])
        #contact_map_row = np.delete(contact_map_row, 0, 0)
        contact_map = np.vstack([contact_map, contact_map_row])
        contact_map = np.vstack([contact_map, np.zeros((10, total_res+num_chain*10+1))])
    return contact_map

    
def contact_map_plot(contact_map):
    fig, axs = plt.subplots(1,2)
    axs[0].imshow(contact_map, cmap='hot', interpolation='nearest')
    axs[0].invert_yaxis()
    axs[1].hist(contact_map)
    plt.savefig(args.output)
    
def main():
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure("protein", args.pdb)
    contact_map = contact_map_total(structure)
    contact_map_plot(contact_map)

if __name__ == '__main__':
    main()