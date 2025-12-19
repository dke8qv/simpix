# simpix

C++ starter code
* simpix_start.cpp
use make to build this example

Usage: simapix_start image1 image2 <output=out.png>

Python starter code
* simpix_start.py

Usage: simapix_start image1 image2 <output=out.png>



## Simpix (Simulated Annealing Pixel Swap)

This program takes a source image A and rearranges (swaps) its pixels using simulated annealing to match a target image B.
Each source pixel is used exactly once.

## Images + Resizing

I selected two images A and B and resized them to 1920x1080 using ImageMagick:

convert A.jpg -resize 1920x1080\! A_1920x1080.png
convert B.jpg -resize 1920x1080\! B_1920x1080.png

## Build

make

## Run (local)

./simpix A_1920x1080.png B_1920x1080.png out_AtoB.png --batch
./simpix B_1920x1080.png A_1920x1080.png out_BtoA.png --batch

## Run (Rivanna Slurm)

sbatch simpix_AtoB.slurm
sbatch simpix_BtoA.slurm

## Output files

collage_A_1920x1080_to_B_1920x1080.png
collage_B_1920x1080_to_A_1920x1080.png

(out images written via the 3rd command-line argument)

## Runtimes

A -> B: runtime_seconds = 37.574716076
B -> A: runtime_seconds = 34.647599371
