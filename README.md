# MoLPC-IMP
## Introduction
**M**odelling **o**f **L**arge **P**rotein **C**omplexes with **I**ntegrative **M**odeling **P**latform

This directory contains the pipeline for modeling large protein complexes without knowing the stochiometry using AlphaFold2, Monte Carlo Tree Search(MCTS), and IMP. MoLPC-IMP can be run using predictions of subcomponents from any method and is thus not directly dependent on AlphaFold2.

Given a set of unique protein sequences, this pipeline predicts the structure of an entire complex composed of the supplied sequences. The pipeline is developed for protein complexes with 7-26 chains, but is also functional for smaller protein complexes.

### | [1AVO](https://www.rcsb.org/structure/1avo) | C7 | Hetero 14-mer - A7B7 |
**Native complex in grey, prediction colored by chain**

<img src="analysis/1AVO.gif" width="75%" height="75%" />

## Computational Requirement
Before beginning the process of setting up this pipeline on your local system, make sure you have adequate computational resources.
The main bottleneck here is the structure prediction of trimeric subcomponents with AlphaFold2, which can require >40Gb of GPU RAM
depending on the number of residues in the subcomponent that is being predicted. Make sure you have available GPUs suitable for this
type of structure prediction as predicting with CPU will take an unreasonable amount of time. This pipeline assumes you have NVIDIA GPUs
on your system, readily available.

## Setup
All needed packages(including AlphaFold and IMP) are supplied through Conda. The only requirement for running MoLPC is therefore Conda, which can be installed by following: https://docs.conda.io/en/latest/miniconda.html
To setup this pipeline, clone this gitlab repository:
```bash
git clone https://github.com/hychim/molpc-imp.git
```
Run the following to install a conda environment with the necessary dependencies.
```bash
conda env create --name molpc-imp -f environment.yml
```
The AlphaFold database are also required which can be downloaded following: https://github.com/deepmind/alphafold. The database is assumed named as alphafold_data and placed outside the molpc-imp. You can also change it in the pipeline.sh script,
```
molpc-imp
├── data
├── environment.yml
├── output
├── pipeline.sh
├── README.md
├── script
└── src
alphafold_data
├── bfd
├── mgnify
├── params
├── pdb70
├── pdb_mmcif
├── pdb_seqres
├── uniprot
├── uniref30
└── uniref90
```

## Pipeline
To run the molpc-imp pipeline. The only inut required is a fasta file containing the sequence of all individual chains in the complexes.

For example, to fold a homomer. The input fasta should be:
```
>ID_0
<SEQUENCE>
```

For example, to fold a heteromer with 2 individual chains. The input fasta should be:
```
>ID_0
<SEQUENCE A>
>ID_1
<SEQUENCE B>
>ID_2
<SEQUENCE C>
```
Then to run molpc-pipeline, simply do:
```
bash pipeline YOUR_FASTA.fasta
```

### Help
```
Please make sure all required parameters are given
Usage: pipeline.sh <OPTIONS>
Required Parameters:
-f <fasta_path>       Path to fasta file with all unique chain(s), e.g. data/2BL2.fasta
Optional Parameters:
-o <output>           Output path(default: ./output)
-m <mer>              Number of chains of the sub-unit predicted from the AlphaFold. (default: 3)
-d <alphafold_data>   Path to directory of AlphaFold supporting data. (default: ../alphafold_data_v2.3)
-c <moves>            Maximum moves in monte carlo tree search, if your complexes have more than 30 chains, please increase the no. of moves. (default: 30)
-s <steps>            Number of simulations in each moves in mcts, more the steps, more accurate the modeling will be. (default: 50)
-r <remodel>          Remodel the final structure with AlphaFold (AF), IMP (IMP) or no re-modeling(False) (default: 'False')

```

### MoLPC-IMP output
The outputs will be saved in output directory(/molpc/output/<NAME>/) of the directory. The outputs include the computed MSAs, unrelaxed structures, relaxed structures, ranked structures, raw model outputs, prediction metadata, and section timings. The output directory will have the following structure:

```
<NAME>
├── <NAME>_imp.pdb
├── <NAME>_mcts.pdb
├── <NAME>_stiochiometry.fasta
├── fasta_trimer
│   └── <NAME>_{XXX}.fasta
├── imp
│   ├── <NAME>_dr.csv
│   └── <NAME>_topology.txt
├── mcts
│   ├── <NAME>_close.txt
│   ├── <NAME>_final.pdb
│   ├── <NAME>_path.txt
│   └── <NAME>_step{1,2,...,30}.pdb
├── pairs
│   └── <NAME>_{XXX}_{YY}_{ZZ}.pdb
└── trimer
    └── <NAME>_{XXX}.pdb
```
The contents of each output file are as follows:
 - `<NAME>_imp.pdb` - Final model from IMP
 - `<NAME>_mcts.pdb` - Final model from the Monte Carlo Tree Search
 - `<NAME>_stoichiometry.fasta` - Predicted stoichiometry from MCTS in fasta format
 - `fasta_trimer/<NAME>_{X-X-X}.fasta` - Combination of the input fasta. The <NAME> is the input fasta file name. The {X-X-X is the combination of unique chain. For example, for trimer that have 1 unique chain 0 and 2 unique chain 1, {X-X-X} will be {0-1-1}.
 - `trimer/<NAME>_{X-X-X}_{rank}.pdb` - Trimer structure predicted from AlphaFold. The {X-X-X} is the combination of unique chain. The {rank} is the rank order from AlphaFold.
 - `pairs<NAME>_{X-X-X}_{rank}_{Y-Y}_{Z-Z}.pdb` - Protein pairs that extracted from trimer. {Y-Y} is the interacting unique chain pairs, for example {0-0}. {Z-Z} is the interacting chain name, for example {B-C}.
 - `mcts/<NAME>_move{1,2,...,30}.pdb` - Every moves' result from MCTS.
 - `mcts/<NAME>_close.txt` - 'Closing protein pairs' that detected from MCTS. Needed for better IMP modeling expecially for cyclic protein.
 - `mcts/<NAME>_path.txt` - Asseble path from MCTS.
 - `mcts/<NAME>_final.pdb` - Final pdb structure from MCTS. Same as `<NAME>_mcts.pdb`
 - `imp/<NAME>_dr.csv` - Distant restraint information extracted from MCTS final structure and closing protein pairs.
 - `imp/<NAME>_topology.txt` - Toplogy file for IMP. For details, please refer to [IMP manual](https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html)
 
## Procedure
The pipeline consists of four steps:
1. Compute all combination without replacement
    Combination without replacement of trimer sequences will be computed and generate the fasta file for trimer modeling in AlphaFold2.
1. AlphaFold-Multimer
    Model the trimer combinations with AlphaFold2-Multimer
1. Converting trimer to protein pairs
    Convert the trimer from AlphaFold2 to protein pairs. Pairs with distance larger than the threshold(default 6Å) will be filtered out.
1. MCTS
    From the interactions in the predicted subcomponents, we add chains sequentially following a predetermined path through the interaction network (graph). If two pairwise interactions are A-B and B-C, we assemble the complex A-B-C by superposing chain B from A-B and B-C using BioPython’s SVD and rotating the missing chain to its correct relative position. To find the optimal assembly route for a complex, we search for an optimal path using Monte Carlo Tree Search
1. Extract distant restrain from MCTS final complexes
    
1. Model complexes with IMP
    
