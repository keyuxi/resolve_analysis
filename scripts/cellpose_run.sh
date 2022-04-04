#!/bin/bash
#SBATCH --job-nam=run_cellpose
#SBATCH --output=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf,normal
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00

source ~/.bashrc
conda activate cellpose
python resolve/run_cellpose.py \
    --diameter 65 --mask_threshold -0.2 \
    -i /home/users/kyx/gs/resolve/data/32801-slideC/32801-slideC_A1_DAPI.tiff \
    -o /home/users/kyx/gs/resolve/data/cellpose_32801-slideC/cellpose_mask_A1.tiff