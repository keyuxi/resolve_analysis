#!/bin/bash
#SBATCH --job-nam=run_baysor
#SBATCH --output=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf,normal
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00

/oak/stanford/groups/wjg/kyx/software/Baysor/bin/Baysor run -s 45.0 --scale-std=50% -c config.toml -p --exclude-genes='Adra2a' --force-2d -o ./A1_out 32801-slideC_A1_results.txt
