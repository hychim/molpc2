import itertools
import argparse
import os

parser = argparse.ArgumentParser(description = '''Write the fasta file in to all permutation.''')
parser.add_argument('--fasta', type=str, help = 'Path to fasta file.')
parser.add_argument('--ID', type=str, help = 'job ID')
parser.add_argument('--output', type=str, help = 'Output directory')
args = parser.parse_args()

fasta = args.fasta
ID = args.ID
output = args.output

def gen_perm(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    lst = range(len(lines[::2]))
    #perm = [p for p in itertools.product(lst, repeat=3)]
    perm = [p for p in itertools.combinations_with_replacement(lst, 3)]
    return perm
    
def gen_perm_fasta(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    for i in gen_perm(filename):
        fasta = []
        for j in i:
            fasta.append(lines[j*2])
            fasta.append(lines[j*2+1])
        name = ''.join(str(x) for x in i)
        os.makedirs(output) if not os.path.exists(output) else None
        with open(f'{output}/{ID}_{name}.fasta', 'w') as f:
            for line in fasta:
                f.write(f"{line}\n")

def main():
    gen_perm_fasta(fasta)

if __name__ == "__main__":
    main()
