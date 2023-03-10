import argparse
import os
import numpy as np
import glob
import copy
from collections import Counter, defaultdict
from Bio.PDB import PDBParser, Superimposer, PDBIO
import string

import argparse

argument_parser = argparse.ArgumentParser(description = '''Search for assemble path''')
argument_parser.add_argument('--id', type=str, help = 'Path to pdb files.')
argument_parser.add_argument('--pairs_dir', type=str, help = 'Path to pdb files.')
argument_parser.add_argument('--steps', type=int, help = 'No. of simluation steps for each moves.')
argument_parser.add_argument('--moves', type=int, help = 'No. of moves.')
argument_parser.add_argument('--output', type=str, help = 'Output')
args = argument_parser.parse_args()

# define pdb parser and superimposer for biopython pdb
parser = PDBParser()
super_imposer = Superimposer()
# define all sources and edges
output = args.output
pairs_dir = args.pairs_dir + "*"
pair_paths = glob.glob(pairs_dir)

edges = []
edges_pdb = []
sources = pair_paths
for file_path in pair_paths:
    name = os.path.basename(file_path)[:-4]
    edges_lst = (name.split('_')[3]).split('-')
    edges_pdb_lst = (name.split('_')[4]).split('-')
    edges.append(edges_lst)
    edges_pdb.append(edges_pdb_lst)

    # bidriection
    edges_lst_rev = edges_lst[::-1]
    edges.append(edges_lst_rev)
    edges_pdb_lst_rev = edges_pdb_lst[::-1]
    edges_pdb.append(edges_pdb_lst_rev)

edges = np.array(edges)
edges_pdb = np.array(edges_pdb)
sources = np.array(sources)
sources = np.repeat(sources, 2)

def superimpose(structureA, structureB, shared_chain_inA, shared_chain_inB, added_chain):
    ref_atoms = []
    alt_atoms = []

    for res in structureA[0][shared_chain_inA]:
        ref_atoms.append(res["CA"])
    for res in structureB[0][shared_chain_inB]:
        alt_atoms.append(res["CA"])

    super_imposer.set_atoms(ref_atoms, alt_atoms)
    super_imposer.apply(structureB.get_atoms())

    structure_merged  = structureA.copy()
    structure_merged[0].add(structureB[0][added_chain])
    
    return structure_merged

def dist_comp(chainA, chainB):
    contact_map = []
    
    for r_A in chainA:
        dists = []
        for r_B in chainB:
            dists.append(r_A["CA"] - r_B["CA"])
        contact_map += dists
    return contact_map

def check_overlaps(structure):
    # check only the last chain
    threshold = 3
    overlap_percent_threshold = 0.1
    close_end = None
    
    for chain in structure[0]:
        last_chain_id = chain.get_id()

    for chain in structure[0]:
        if structure[0][last_chain_id].get_id() != chain.get_id():
            dist = dist_comp(structure[0][last_chain_id], chain)
            overlap_percent = len([i for i in dist if i<threshold])/(len(chain))
            if overlap_percent > overlap_percent_threshold:
                if 0.9> overlap_percent > 0.5:
                    #print(overlap_percent)
                    #print(f'close end: {structure[0][last_chain_id].get_id()} and {chain.get_id()}')
                    close_end = chain.get_id()
                return True, close_end
    return False, close_end

def count_interface(structure):
    num_interface = 0
    for chainA in structure[0]:
        for chainB in structure[0]:
            if chainA.get_id() != chainB.get_id():
                contact_map = dist_comp(chainA, chainB)
                if 8 > min(contact_map): # 8 is the threshold for interacting protein
                    num_interface += 1
    return num_interface/2

def count_interface_chain(structure):
    num_interface_lst = []
    for chainA in structure[0]:
        num_interface = 0
        for chainB in structure[0]:
            if chainA.get_id() != chainB.get_id():
                contact_map = dist_comp(chainA, chainB)
                if 12 > min(contact_map): # 8 is the threshold for interacting protein
                    num_interface += 1
        num_interface_lst.append(num_interface)
    return num_interface_lst

def re_name_chain(structure):
    new_structure = structure.copy()
    i = 0

    for chain in new_structure[0]:
        chain.id = string.ascii_lowercase[i]
        i += 1
    return new_structure

def get_plddt(structure):
    plddt = []
    for chain in structure[0]:
        for res in chain:
            plddt.append(res["CA"].get_bfactor())
    return plddt
        
def score_complex(structure):
    '''Score all interfaces in the current complex
    '''
    plddt = get_plddt(structure)
    complex_score =  np.log10(count_interface(structure))*(sum(plddt)/len(plddt))
    return complex_score

def save_pdb(structure, outpath):
    io = PDBIO()
    io.set_structure(structure)
    io.save(outpath)

class MonteCarloTreeSearchNode():
    def __init__(self, chain, edge_chain, source=None, structure=None, complex_scores=[0], parent=None, 
                parent_path=[], parent_pdb_path=[]):
        self.chain = chain                  # chain that decided to be added in this decision
        self.edge_chain = edge_chain        # chain that new chain "added" on
        self.source = source                # where the chain comes from
        self.structure = structure          # the structure exist as a biopython pdb class
        self.complex_scores = complex_scores

        self.parent = parent #Parent node
        self.path = copy.deepcopy(parent_path) #All nodes up to (and including) the parent
        self.path.append(chain)

        self.pdb_path = copy.deepcopy(parent_pdb_path)
        self.pdb_path.append(string.ascii_lowercase[len(self.path)-1])

        self.children = [] #All nodes branching out from the current
        self._number_of_visits = 0
        
        #self._untried_edges, self._untried_edges_pdb, self._untried_sources, self._untried_edgesA = self.get_possible_edges_all()

        if self.structure == None:
            self._untried_edges, self._untried_edges_pdb, self._untried_sources, self._untried_edgesA = self.get_possible_edges_all()
        else:
            self._untried_edges, self._untried_edges_pdb, self._untried_sources, self._untried_edgesA = self.get_possible_edges()

        self.early_stop = False
        self.close_end = None
        return

    def get_possible_edges_all(self):
        untried_edges = []
        untried_edges_pdb = []  # for pdb chain, eg. ["B", "C"], ["C", "D"]
        untried_sources = []
        untried_edgesA = []     # for re named chain path, eg. ["0" ,"0"], ["0" ,"1"]
        for j in range(len(self.path)):
            #Get all edges to the current node
            cedges = edges[np.argwhere(edges[:,0]==self.path[j])[:,0]]
            cedges_pdb = edges_pdb[np.argwhere(edges[:,0]==self.path[j])[:,0]]
            csources = sources[np.argwhere(edges[:,0]==self.path[j])[:,0]]
            for i in range(len(cedges)):
                untried_edges.append(cedges[i])
                untried_edges_pdb.append(cedges_pdb[i])
                untried_sources.append(csources[i])
                untried_edgesA.append(self.pdb_path[j])

        return untried_edges, untried_edges_pdb, untried_sources, untried_edgesA

    def get_possible_edges(self):
        untried_edges = []
        untried_edges_pdb = []  # for pdb chain, eg. ["B", "C"], ["C", "D"]
        untried_sources = []
        untried_edgesA = []     # for re named chain path, eg. ["0" ,"0"], ["0" ,"1"]
        # need to fix none type
        interface_lst = count_interface_chain(self.structure)
        shortlisted = []
        shortlisted_pdb_path = []
        for i in range(len(self.path)):
            if interface_lst[i] <= min(interface_lst)+1:
                shortlisted.append(self.path[i])
                shortlisted_pdb_path.append(self.pdb_path[i])

        for j in range(len(shortlisted)):
            #Get all edges to the current node
            cedges = edges[np.argwhere(edges[:,0]==shortlisted[j])[:,0]]
            cedges_pdb = edges_pdb[np.argwhere(edges[:,0]==shortlisted[j])[:,0]]
            csources = sources[np.argwhere(edges[:,0]==shortlisted[j])[:,0]]
            for i in range(len(cedges)):
                untried_edges.append(cedges[i])
                untried_edges_pdb.append(cedges_pdb[i])
                untried_sources.append(csources[i])
                untried_edgesA.append(shortlisted_pdb_path[j])

        return untried_edges, untried_edges_pdb, untried_sources, untried_edgesA

    def expand(self):
        new_edge = self._untried_edges.pop()
        new_edge_pdb = self._untried_edges_pdb.pop()
        new_source = self._untried_sources.pop()
        new_edgeA = self._untried_edgesA.pop()

        chainA = new_edge_pdb[0]    # chain name in pairs pdb, e.g. B, C, D
        chainA_renamed = new_edgeA  # renamed chain name, e.g. a, b, c, d
        chainB = new_edge_pdb[1]    # chain name in pairs pdb, e.g. B, C, D

        chainA_ind = new_edge[0]    # edge chain as individual chain name, e.g. 0, 1, 2
        chainB_ind = new_edge[1]    # added chain as individual chain name, e.g. 0, 1, 2


        if self.structure == None:
            #first node
            child_structure = parser.get_structure("child", new_source)
        elif self.structure != None:
            source_structure = parser.get_structure("child", new_source)
            child_structure = superimpose(self.structure, source_structure, chainA_renamed, chainA, chainB)

        child_structure = re_name_chain(child_structure)
        overlap, close_end = check_overlaps(child_structure)

        if not overlap:
            complex_score = score_complex(child_structure)
            child_node = MonteCarloTreeSearchNode(chainB_ind, chainA_ind, source=new_source, 
                                                    structure=child_structure, complex_scores=[complex_score], 
                                                    parent=self, parent_path=self.path, parent_pdb_path=self.pdb_path)
        elif overlap:
            if close_end != None:
                self.close_end = f'{chainA_renamed}, {close_end}'
                self.close_end_source = new_source
            return self

        self.children.append(child_node)
        return child_node

    def rollout(self):
            '''Simulate an assembly path until
            1. all chains are in complex
            2. an overlap is found
            '''
            overlap = False

            path_node = copy.deepcopy(self)

            while len(path_node.path)<26 and overlap==False:
                if len(path_node._untried_edges)>0:
                    edge_ind = np.random.randint(len(path_node._untried_edges))
                else:
                    overlap=True
                    break

                new_edge = path_node._untried_edges[edge_ind]
                new_edge_pdb = path_node._untried_edges_pdb[edge_ind]
                new_source = path_node._untried_sources[edge_ind]
                new_edgeA = path_node._untried_edgesA[edge_ind]

                chainA = new_edge_pdb[0]
                chainA_renamed = new_edgeA
                chainB = new_edge_pdb[1]

                chainA_ind = new_edge[0]    # edge chain as individual chain name, e.g. 0, 1, 2
                chainB_ind = new_edge[1]    # added chain as individual chain name, e.g. 0, 1, 2

                if path_node.structure != None:
                    source_structure = parser.get_structure("child", new_source)
                    child_structure = superimpose(path_node.structure, source_structure, chainA_renamed, chainA, chainB)
                elif path_node.structure == None: #root node
                    child_structure = parser.get_structure("child", new_source)
                    
                child_structure = re_name_chain(child_structure)
                
                #Check overlaps
                overlap, close_end = check_overlaps(child_structure)

                #If no overlap - score and create a child node
                if overlap==False:
                    path_node = MonteCarloTreeSearchNode(chainB_ind, chainA_ind, source=new_source, structure=child_structure, 
                                                            complex_scores=[0], parent=path_node, parent_path=path_node.path, 
                                                            parent_pdb_path=path_node.pdb_path)
                elif overlap==True:
                    break
            #Score rollout
            if path_node.structure != None:
                rollout_score = score_complex(path_node.structure)
            elif path_node.structure == None:
                rollout_score = 0
            return rollout_score

    def back_prop(self, rollout_score):
        '''Update the previous nodes in the path
        '''
        self._number_of_visits += 1
        self.complex_scores.append(rollout_score)
        #This is recursive and will back_prop to all parents
        if self.parent:
            self.parent.back_prop(rollout_score)

    def best_child(self):
        '''Calculate the UCB

        Vi is the average reward/value of all nodes beneath this node (sum of interface scores)
        N is the number of times the parent node has been visited, and
        ni is the number of times the child node i has been visited

        The first component of the formula above corresponds to exploitation;
        it is high for moves with high average win ratio.
        The second component corresponds to exploration; it is high for moves with few simulations.
        '''
        choices_weights = [(np.average(c.complex_scores) + 2 * np.sqrt(np.log(c.parent._number_of_visits+1e-3) / (c._number_of_visits+1e-12))) for c in self.children]
        return self.children[np.argmax(choices_weights)]

    def tree_policy(self):
        current_node = self
        fully_expanded = (len(current_node._untried_edges) == 0)
        while len(current_node.path) < 26:
            if fully_expanded:
                if len(current_node.children) != 0:
                    current_node = current_node.best_child()
                elif len(current_node.children)==0:
                    current_node.early_stop = True
                return current_node
            elif not fully_expanded:
                return current_node.expand()

    def best_action(self):
        simulation_no = args.steps
        
        for i in range(simulation_no):
            v = self.tree_policy()
            if v != None: # if not None
                reward = v.rollout()
                v.back_prop(reward)
                print(v.path)
        if len(self.children) > 0:
            return self.best_child()
        else:
            return self

def main():
    root = MonteCarloTreeSearchNode('0', '', source=None, structure=None, complex_scores=[0], parent=None, parent_path=[], parent_pdb_path=[])

    v = root

    move_count = 1
    
    if not os.path.exists(f'{output}'):
        os.makedirs(f'{output}')
    
    for _ in range(args.moves):
        v = v.best_action()
        save_pdb(v.structure, f'{output}{args.id}_step{str(move_count)}.pdb')
        move_count += 1
        if v.early_stop:
            print("early stop")
            break

    save_pdb(v.structure, f'{output}{args.id}_final.pdb')

    # writing output csv
    pairs_lst = []
    source_lst = []

    if v.close_end != None:
        print('close end', v.close_end, v.close_end_source)
        with open(f'{output}{args.id}_close.txt', "w") as text_file:
            text_file.write(f'{v.close_end},{v.close_end_source}')

    for i in range(len(v.pdb_path)-1):
        pairs_lst.append(f'{v.pdb_path[i]},{v.pdb_path[i+1]}')
    if v.close_end != None:
        pairs_lst.append(v.close_end)
        source_lst.append(v.close_end_source)

    w = v
    
    while w.source:
        source_lst.append(w.source)
        w = w.parent
    source_lst.reverse()

    out_lst = []
    out_lst.append(f'output/{args.id}/mcts/{args.id}_final.pdb')
    out_lst.append(" ".join(v.pdb_path))
    out_lst.append(" ".join(v.path))
    out_lst = out_lst + [f'{pairs_lst[i]},{source_lst[i]}' for i in range(len(source_lst))]

    with open(f'{output}{args.id}_path.txt', 'w') as f:
        for line in out_lst:
            f.write(f"{line}\n")

if __name__ == "__main__":
    main()
