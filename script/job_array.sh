#!/bin/bash
#SBATCH --gpus 1
#SBATCH -t 3-00:00:00
#SBATCH -A Berzelius-2022-230
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --mail-user ho-yeung.chim.2766@student.uu.se
#SBATCH --array=1-151%75

module load Anaconda/2021.05-nsc1
conda activate /proj/berzelius-2021-29/users/x_hoych/conda_molpc-imp

cd /proj/berzelius-2021-29/users/x_hoych/molpc-imp

config=molpc_array_sub.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "Array Task ID: ${ArrayTaskID}, pdb_id: ${sample}"
bash pipeline.sh -f data/molpc_dataset/$sample'.fasta' -m 4
