#!/bin/bash
#SBATCH --job-name=cheb
#SBATCH --array=1-8192
#SBATCH --partition=shared
#SBATCH --time=01-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G ## max amount of memory per node you require
##SBATCH --core-spec=0 ## Uncomment to allow jobs to request all cores on a node    
#SBATCH --error=cheb-%A_%a.err ## %A - filled with jobid
#SBATCH --output=cheb-%A_%a.out ## %A - filled with jobid
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
cd /home/yaguangl/koa_scratch/modelbase/reimer-mass-loss/hpc/

# activate astro
micromamba activate astro

# think $SLURM_ARRAY_TASK_ID as a row number
# extract the number in that row of batch_index.lst
# this number will be the track index 
fid=$(sed -n "$SLURM_ARRAY_TASK_ID"p batch_index.lst)
cp -r work_cheb "run_cheb_${fid}"

# run mesa
cd "run_cheb_${fid}"
sh clean
sh mk
index=$(printf "%06d" $fid)
cp ../final_model/index${index}_zacheb.mod start.mod
python driver.py ${fid}

# upon completion
# do not use "" when use *
mv mesa_terminal_output* ../log/
mv *.mod ../final_model/
mv *.h5 ../history/
cd ../

rm -r "run_cheb_${fid}"

# record time
date
