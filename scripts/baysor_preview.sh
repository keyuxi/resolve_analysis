#!/bin/bash
#SBATCH --job-nam=preview_baysor
#SBATCH --output=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.out
#SBATCH --error=/scratch/groups/wjg/kyx/resolve/out/R-%x.%j.err
#SBATCH -n 1
#SBATCH --partition=wjg,biochem,sfgf,normal
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

/oak/stanford/groups/wjg/kyx/software/Baysor/bin/Baysor preview -c config.toml -o A1_out --exclude-genes='Adra2a' 32801-slideC_A1_results.txt
