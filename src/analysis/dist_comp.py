import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import argparse

argument_parser = argparse.ArgumentParser(description = '''Plotting protein contact map from pdb files''')
argument_parser.add_argument('--pdb1', type=str, help = 'Path to the experimetal pdb files.')
argument_parser.add_argument('--pdb2', type=str, help = 'Path to the your model pdb files.')
argument_parser.add_argument('--mmalign', type=str, help = 'Path to mmalign result(outfmt) txt file.')
argument_parser.add_argument('--output', type=str, help = 'Path for output.')
args = argument_parser.parse_args()

plt.rcParams['figure.figsize'] = [8, 8]

def cal_dist(r1, r2):
        # r1 means residue 1, r2 means residue 2
        return r1["CA"] - r2["CA"]

def clean_pdb(structure):
    new_structure = {}

    for chain in structure[0]:
        new_chain = []
        for r in chain:
            if r.id[0] == " ":
                new_chain.append(r)
        new_structure[chain.get_id()] = new_chain
    return new_structure
    
def clean_pdb_fixed(structure):
    # for fixed pdb files from perl script only(only protein model needed this step
    # to fix the pdb numbering)
    new_structure = {}

    for chain in structure[0]:
        new_chain = []
        for r in chain:
            if r.id[2] == " ":
                new_chain.append(r)

        new_structure[chain.get_id()] = new_chain
    return new_structure

def contact_map_btw(chainA, chainB):
    contact_map = np.zeros((1,len(chainB)))
    
    for r_A in chainA:
        dists = []
        for r_B in chainB:
            dists.append(r_A["CA"] - r_B["CA"])
        contact_map = np.vstack([contact_map, dists])
    return np.delete(contact_map, 0, 0)

def contact_map_total_raw(structure, chain_lst, native):
    #borderline = 10
    total_res = 0
    num_chain = 0

    if native:
        structure = clean_pdb(structure)
    else:
        structure = clean_pdb_fixed(structure)
    
    for chain_id in chain_lst:
        total_res += len(structure[chain_id])
        num_chain += 1

    print("total res:", total_res)
    print("num chain:", num_chain)

    contact_map = np.zeros((1, total_res+1))
    for chainA_id in chain_lst:
        chainA = structure[chainA_id]
        contact_map_row = np.zeros((len(chainA), 1))
        for chainB_id in chain_lst:
            chainB = structure[chainB_id]
            contact_map_sub = contact_map_btw(chainA, chainB)
            contact_map_row = np.hstack([contact_map_row, contact_map_sub])
        contact_map = np.vstack([contact_map, contact_map_row])
    return np.delete(np.delete(contact_map, 0, 0),0,1)

def dist_comp(stru1, lst1, stru2, lst2):
    # stru1 should be the experimental structure
    # stru2 should be your model
    if len(lst1)>len(lst2):
        lst1 = lst1[:len(lst2)]
    elif len(lst2)>len(lst1):   
        lst2 = lst2[:len(lst1)]

    print(lst1)
    print(lst2)

    contact_map1 = contact_map_total_raw(stru1, lst1, True)
    contact_map2 = contact_map_total_raw(stru2, lst2, False)

    print(contact_map1.shape)
    print(contact_map2.shape)

    plt.rcParams['figure.figsize'] = [7, 7]

    plt.plot(contact_map1.flatten(), contact_map2.flatten(), ',')
    plt.xlabel("experimental")
    plt.ylabel("model")

    return None

def main():
    with open(args.mmalign) as f:
        lines = f.readlines()

    parser = Bio.PDB.PDBParser()
    structure1 = parser.get_structure("experimental", args.pdb1)
    structure2 = parser.get_structure("model", args.pdb2)

    exper_chain_lst = lines[1].split("\t")[0].split(":")[1:]
    model_chain_lst = lines[1].split("\t")[1].split(":")[1:]
    # remove empty str in the list(due to unmatch number of chain in pdbs)
    exper_chain_lst = [x for x in exper_chain_lst if x]
    model_chain_lst = [x for x in model_chain_lst if x]

    dist_comp(structure1, exper_chain_lst, structure2, model_chain_lst)
    plt.savefig(args.output)

if __name__ == "__main__":
    main()