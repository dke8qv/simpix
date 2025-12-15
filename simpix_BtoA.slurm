#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH -A phys56xx
#SBATCH --partition=standard
#SBATCH --job-name=simpixBtoA
#SBATCH --output=simpixBtoA_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

cd $SLURM_SUBMIT_DIR
make

/usr/bin/time -p ./simpix B_1920x1080.png A_1920x1080.png out_BtoA.png --batch
echo "done"
