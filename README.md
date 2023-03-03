# MoLPC-IMP
## Introduction
**M**odelling **o**f **L**arge **P**rotein **C**omplexes with **I**ntegrative **M**odeling **P**latform

This directory contains the pipeline for modeling large protein complexes without knowing the stochiometry using AlphaFold2, Monte Carlo Tree Search(MCTS), and IMP. MoLPC-IMP can be run using predictions of subcomponents from any method and is thus not directly dependent on AlphaFold2.

Given a set of unique protein sequences, this pipeline predicts the structure of an entire complex composed of the supplied sequences. The pipeline is developed for protein complexes with 7-26 chains, but is also functional for smaller protein complexes.

### |1AVO | TMscores | A7B7 |
**Native complex in grey, prediction colored by chain**
<img src="1AVO.gif" width="40" height="40" />


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

### molpc-imp output
The outputs will be saved in output directory(/molpc/output/<NAME>/) of the directory. The outputs include the computed MSAs, unrelaxed structures, relaxed structures, ranked structures, raw model outputs, prediction metadata, and section timings. The output directory will have the following structure:

```
<NAME>
├── <NAME>_imp.pdb
├── <NAME>_mcts.pdb
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
 - `<NAME>_imp.pdb`  - Final model from IMP
 - `<NAME>_mcts.pdb` - Final model from the Monte Carlo Tree Search
 
 
## Procedure
The pipeline consists of four steps:
1. Compute all combination without replacement

2. AlphaFold-Multimer

3. Converting trimer to protein pairs

4. MCTS

5. Extract distant restrain from MCTS final complexes

6. Model complexes with IMP
