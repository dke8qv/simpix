#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH -A phys56xx
#SBATCH --partition=standard
#SBATCH --job-name=simpix1080
#SBATCH --output=simpix1080_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

source ~/.bashrc

make

# run both directions
./simpix A_1920x1080.png B_1920x1080.png out_AtoB.png --batch
./simpix B_1920x1080.png A_1920x1080.png out_BtoA.png --batch

echo "done"


