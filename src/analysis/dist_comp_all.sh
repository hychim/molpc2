for file in ../../data/1IJG/pairs/*
do
    echo $file
    ./fix_numbering.pl $file data/1IJG.pdb
    ./MMalign data/1IJG.pdb $file'.fixed' -outfmt 2 > mmalign.txt
    python dist_comp.py --pdb1 data/1IJG.pdb --pdb2 $file'.fixed' --mmalign mmalign.txt --output ../../result/analysis/dist_comp_dimer_1IJG/$(basename "${file%.*}")'.png'
    rm $file'.fixed'
done
#./MMalign $1 $2 -outfmt 2 > mmalign.txt
#python dist_comp.py --pdb1 $1 --pdb2 $2 --mmalign mmalign.txt --output $3