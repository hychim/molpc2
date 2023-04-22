import argparse
import os

parser = argparse.ArgumentParser(description = '''Write the fasta file in to all permutation.''')
parser.add_argument('--input', type=str, help = 'Path to mcts output txt.')
parser.add_argument('--fasta', type=str, help = 'fasta file with all individual chains.')
parser.add_argument('--ID', type=str, help = 'job ID')
parser.add_argument('--output', type=str, help = 'Output directory')
args = parser.parse_args()

with open(args.input, 'r') as f:
    lines = f.read().splitlines()
with open(args.fasta, 'r') as f:
    fasta = f.read().splitlines()

stoi = lines[2].split()

with open(f'{args.output}/{args.ID}_stoichiometry.fasta', 'w') as f:
            for indiv in stoi:
                f.write(f"{fasta[int(indiv)*2]}\n")
                f.write(f"{fasta[int(indiv)*2+1]}\n")