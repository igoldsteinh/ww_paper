#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 4:00:00   ## 4 hr run time limit
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu


module purge
module load R
cd //dfs6/pub/igoldst1/ww_paper

sim_num=1
Rscript scripts/process_results/summarise_eir_cases.R $sim_num

