FASTA=$1

FILENAME=$(basename -- "$FASTA")
NAME="${FILENAME%.*}"
echo $NAME

#mkdir output/$NAME

### Fasta Preprocessing
python src/preprocessing/af_fasta_pre.py --ID $NAME --fasta $FASTA --output output/$NAME/fasta_trimer/

### AlphaFold
echo 'Start modeling trimer with alphafold'
for T_FASTA in output/$NAME/fasta_trimer/*
do
	#mkdir output/$NAME/trimers/
	echo $T_FASTA
	#bash run_alphafold.sh
done
echo 'Trimer modeling done'

### Dimer2Pairs
echo 'Converting trimers to protein pairs'
python src/preprocessing/write_pairs.py --trimer_dir output/$NAME/trimer/ --output_dir output/$NAME/pairs/
echo 'Pairs converting done'

### MCTS
echo 'Start MCTS'
python src/mcts/mcts.py --id $NAME --pairs_dir output/$NAME/pairs/ --output output/$NAME/mcts/
echo 'MCTS done'

### IMP
echo 'Start IMP'
# convert the mcts final structure and the close end pairs to distant restraint
python src/imp/pdb2dr.py --id $NAME --close_end output/$NAME/mcts/$NAME'_close.txt' --mcts_pdb output/$NAME/mcts/$NAME'_final.pdb' --dist 6
# generate the topology file for IMP
python src/imp/gen_topo.py --id $NAME --input output/$NAME/mcts/$NAME'_path.txt'
# Running IMP
python src/imp/modeling.py --topology output/$NAME/imp/$NAME'_topology.txt' --fasta data/ --pdbdir . --dr output/$NAME/imp/$NAME'_dr.csv' --steps 20000 --outdir output/$NAME/imp/imp

echo 'IMP done'


cp output/$NAME/mcts/$NAME'_final.pdb' output/$NAME/$NAME'_mcts.pdb'
cp output/$NAME/imp/imp/pdbs/model.0.pdb output/$NAME/$NAME'_imp.pdb'
rm -r output/$NAME/imp/imp

echo 'MoLPC done'
