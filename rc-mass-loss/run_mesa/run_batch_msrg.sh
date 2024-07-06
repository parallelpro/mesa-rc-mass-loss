#!/bin/bash
#SBATCH --job-name=msrg
#SBATCH --array=1-100
#SBATCH --partition=shared
#SBATCH --time=03-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    
#SBATCH --error=msrg-%A_%a.err ## %A - filled with jobid
#SBATCH --output=msrg-%A_%a.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=yaguangl@hawaii.edu

## All options and environment variables found on schedMD site: http://slurm.schedmd.com/sbatch.html

# record time
date
hostname

# change to zsh
# module purge
source /home/yaguangl/custom_setup.sh
source /home/yaguangl/.zshrc

# navigate to the mesa directory
cd /home/yaguangl/koa_scratch/modelbase/rc-mass-loss/hpc/

# activate astro
micromamba activate astro

# think $SLURM_ARRAY_TASK_ID as a row number
# extract the number in that row of batch_index.lst
# this number will be the track index 
fid=$(sed -n "$SLURM_ARRAY_TASK_ID"p batch_index.lst)
cp -r work_ms_rg "run_ms_rg_${fid}"

# run mesa
cd "run_ms_rg_${fid}"
sh clean
sh mk
python driver.py ${fid}

# upon completion
# do not use "" when use *
mv mesa_terminal_output* ../log/
mv *.mod ../final_model/
mv *.h5 ../history/
cd ../

rm -r "run_ms_rg_${fid}"

# record time
date
