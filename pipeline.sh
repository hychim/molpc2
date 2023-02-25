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
