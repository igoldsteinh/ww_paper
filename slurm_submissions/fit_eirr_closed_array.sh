#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 4:00:00   ## 4 hr run time limit
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=1-100

module purge
module load julia-1_8_5
cd //dfs6/pub/igoldst1/wastewater2
sim_num=101


julia --project --threads 4 scripts/fit_models/fit_eirrc_closed.jl $sim_num $SLURM_ARRAY_TASK_ID

