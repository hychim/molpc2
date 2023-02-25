cd ../../data/MoLPC_source/
for file in *
do
    echo ${file%%.*}
    python ../../src/imp/pdb2dr.py --dimer_dir ../../data/${file%%.*}/pairs/ --dimer_lst ../../data/MoLPC_source/$file --dist 6 --outdir ../dr/${file%%.*}_dr.txt
done
