#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH -A phys56xx
#SBATCH --partition=standard
#SBATCH --job-name=simpixAtoB
#SBATCH --output=simpixAtoB_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

cd $SLURM_SUBMIT_DIR
make

/usr/bin/time -p ./simpix A_1920x1080.png B_1920x1080.png out_AtoB.png --batch
echo "done"
