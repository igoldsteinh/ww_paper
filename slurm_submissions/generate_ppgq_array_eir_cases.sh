#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 4:00:00   ## 4 hr run time limit
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu
#SBATCH --array=1-2


module purge
module load julia-1_8_5
cd //dfs6/pub/igoldst1/ww_paper
sim_num=1

if [ $SLURM_ARRAY_TASK_ID == 1 ]; then
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID slurm_submissions/process_results_array_eir_cases.sh
fi

julia --project --threads 4 scripts/generate_quantities/eir_cases_generate_pp_and_gq.jl $sim_num $SLURM_ARRAY_TASK_ID

