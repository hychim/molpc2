# plot pair distances comparison between 2 different pdb
# $1 is the experimental
# $2 is the model

./fix_numbering.pl $2 $1
./MMalign $1 $2'.fixed' -outfmt 2 > mmalign.txt
python dist_comp.py --pdb1 $1 --pdb2 $2'.fixed' --mmalign mmalign.txt --output $3
