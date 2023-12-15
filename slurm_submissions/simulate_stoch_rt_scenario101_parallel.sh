#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 8:00:00   ## 4 hr run time limit
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=1-100

module purge
module load R
cd //dfs6/pub/igoldst1/ww_paper

Rscript scripts/simulated_data/simulate_stochastic_rt_scenario101.R $SLURM_ARRAY_TASK_ID

