import glob
import itertools
import os

import argparse

argument_parser = argparse.ArgumentParser(description = '''Remove highly similar trimer''')
argument_parser.add_argument('--trimer_dir', type=str, help = 'Path to trimer directory.')
args = argument_parser.parse_args()

perm = [p for p in itertools.combinations(glob.glob(args.trimer_dir+'*'), 2)]

remove_lst = []
for pairs in perm:
    complexA, complexB = pairs[0], pairs[1]
    command = f'src/preprocessing/MMalign {complexA} {complexB} -outfmt 2'
    result = os.popen(command).read()
    tmscore = float(result.split('\n')[1].split('\t')[2])
    if tmscore>0.9:
        remove_lst.append(complexB)

for complex_extra in remove_lst:
    os.system(f'rm {complex_extra}')