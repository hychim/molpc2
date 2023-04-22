#!/bin/bash
############################################################
# Help 
help() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-f <fasta_path>       Path to fasta file with all unique chain(s), e.g. data/2BL2.fasta"
        echo "Optional Parameters:"
		echo "-o <output>           Output path(default: ./output)"
        echo "-m <mer>              Number of chains of the sub-unit predicted from the AlphaFold. (default: 3)"
		echo "-d <alphafold_data>   Path to directory of AlphaFold supporting data. (default: ../alphafold_data_v2.3)"
		echo "-c <moves>            Maximum moves in monte carlo tree search, if your complexes have more than 30 chains, please increase the no. of moves. (default: 30)"
		echo "-s <steps>            Number of simulations in each moves in mcts, more the steps, more accurate the modeling will be. (default: 50)"
		echo "-r <remodel>          Remodel the final structure with AlphaFold (AF), IMP (IMP) or no re-modeling(False) (default: 'False')"
        echo ""
        exit 1
}

############################################################
while getopts h:f:o:m:d:c:s:r: flag
do
    case "${flag}" in
		h) usage;;
        f) FASTA=${OPTARG};;
		o) OUTPUT=${OPTARG};;
        m) MER=${OPTARG};;
        d) AF_DATA=${OPTARG};;
        c) MOVES=${OPTARG};;
        s) STEPS=${OPTARG};;
        r) REMODEL=${OPTARG};;
    esac
done

if [[ "$FASTA" == "" ]] ; then
    help
fi

if [[ "$OUTPUT" == "" ]] ; then
    OUTPUT='./output'
fi

if [[ "$MER" == "" ]] ; then
    MER='3'
fi

if [[ "$AF_DATA" == "" ]] ; then
    AF_DATA='../alphafold_data_v2.3'
fi

if [[ "$MOVES" == "" ]] ; then
    MOVES='30'
fi

if [[ "$STEPS" == "" ]] ; then
    STEPS='50'
fi

if [[ "$REMODEL" == "" ]] ; then
    IMP=FALSE
fi

############################################################
### MAIN PROGRAM
FILENAME=$(basename -- "$FASTA")
NAME="${FILENAME%.*}"
echo $NAME

mkdir $OUTPUT/$NAME

### Fasta Preprocessing
python src/preprocessing/af_fasta_pre.py --ID $NAME --fasta $FASTA --mer $MER --output $OUTPUT/$NAME/fasta_trimer/

### AlphaFold
echo 'Start modeling trimer with alphafold'
for T_FASTA in $OUTPUT/$NAME/fasta_trimer/*
do
	echo $T_FASTA
	bash src/alphafold/run_alphafold.sh -d $AF_DATA -o $OUTPUT/$NAME/alphafold/ -f $T_FASTA -t 2021-11-01 -r false -m multimer -l 1
done
mkdir $OUTPUT/$NAME/trimer/
for FILE in $OUTPUT/$NAME/alphafold/*
do
	FILENAME=$(basename -- "$FILE")
	FILENAME="${FILENAME%.*}"
	cp $FILE/ranked_0.pdb $OUTPUT/$NAME/trimer/$FILENAME'_0.pdb'
	cp $FILE/ranked_1.pdb $OUTPUT/$NAME/trimer/$FILENAME'_1.pdb'
	cp $FILE/ranked_2.pdb $OUTPUT/$NAME/trimer/$FILENAME'_2.pdb'
	cp $FILE/ranked_3.pdb $OUTPUT/$NAME/trimer/$FILENAME'_3.pdb'
	cp $FILE/ranked_4.pdb $OUTPUT/$NAME/trimer/$FILENAME'_4.pdb'
done

rm -r $OUTPUT/$NAME/alphafold/
python src/preprocessing/pairs_filter.py --trimer_dir $OUTPUT/$NAME/trimer/
echo 'Trimer modeling done'

### Dimer2Pairs
echo 'Converting trimers to protein pairs'
mkdir $OUTPUT/$NAME/pairs/
python src/preprocessing/write2pairs.py --trimer_dir $OUTPUT/$NAME/trimer/ --output_dir $OUTPUT/$NAME/pairs/
echo 'Pairs converting done'

### MCTS
echo 'Start MCTS'
python src/mcts/mcts.py --id $NAME --pairs_dir $OUTPUT/$NAME/pairs/ --output $OUTPUT/$NAME/mcts/ --moves $MOVES --steps $STEPS
cp $OUTPUT/$NAME/mcts/$NAME'_final.pdb' $OUTPUT/$NAME/$NAME'_mcts.pdb'
python src/mcts/write_stoi.py --input $OUTPUT/$NAME/mcts/$NAME'_path.txt' --fasta $FASTA --ID $NAME --output $OUTPUT
echo 'MCTS done'

### IMP
if [[$IMP]] ; then
	echo 'Start IMP'
	# convert the mcts final structure and the close end pairs to distant restraint
	python src/imp/pdb2dr.py --id $NAME --close_end $OUTPUT/$NAME/mcts/$NAME'_close.txt' --mcts_pdb $OUTPUT/$NAME/mcts/$NAME'_final.pdb' --dist 6
	# generate the topology file for IMP
	python src/imp/gen_topo.py --id $NAME --input $OUTPUT/$NAME/mcts/$NAME'_path.txt'
	# Running IMP
	python src/imp/modeling.py --topology $OUTPUT/$NAME/imp/$NAME'_topology.txt' --fasta data/ --pdbdir . --dr $OUTPUT/$NAME/imp/$NAME'_dr.csv' --steps 20000 --outdir $OUTPUT/$NAME/imp/imp

	echo 'IMP done'
	cp $OUTPUT/$NAME/imp/imp/pdbs/model.0.pdb $OUTPUT/$NAME/$NAME'_imp.pdb'
	rm -r $OUTPUT/$NAME/imp/imp
fi

echo 'MoLPC done'
