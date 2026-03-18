#!/bin/bash

#SBATCH --job-name=subsample_10_1000_shotgun_5clus_spearman_average
#your HPC account 
#SBATCH --account=dongyunsheng
#RAM
#SBATCH --mem=10GB
# num Task for job
#SBATCH --ntasks=3
# NOT CHANGE
#SBATCH --nodes=1
# you can choose ricerca or gpu
#SBATCH --partition=ricerca
# your log file to tracking the pipeline
#SBATCH --error="/share/project4/home/dongyunsheng/sbatchlogs/subsample_10_1000_shotgun_5clus_spearman_average.err"
#SBATCH --output="/share/project4/home/dongyunsheng/sbatchlogs/subsample_10_1000_shotgun_5clus_spearman_average.out"
# insert your mail
#SBATCH --mail-user=yunsheng.dong@istitutotumori.mi.it
#SBATCH --mail-type=ALL



# specify conda env 
# eval "$(/share/project4/home/dongyunsheng/myprograms/micromamba shell hook -s bash)"
# micromamba activate r_env
# specify the script 

bash myscript_compDisMatrix_folderUnix.sh
# bash myscript_compDisMatrix_single.sh results/spearman_average_shotgun
# ~/myprograms/micromamba run -n r_env Rscript createdbMain.R