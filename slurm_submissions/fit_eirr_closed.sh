#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 4:00:00   ## 4 hr run time limit
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=igoldst1@uci.edu


module purge
module load julia-1_8_5
cd //dfs6/pub/igoldst1/ww_paper
sim_num=1
julia --project --threads 4 scripts/fit_models/fit_eirr_closed.jl $sim_num

sbatch --depend=afterany:$SLURM_JOB_ID slurm_submissions/generate_pp_and_gq_eirr_closed.sh
