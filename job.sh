#!/bin/bash
#SBATCH --job-name=CAD_main
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10gb
#SBATCH --time=28-00:00:00
#SBATCH --partition=batch_30d
#SBATCH --mail-type=ALL
#SBATCH -o logs/%x_%j.out

# Load necessary modules
module load R/4.4.1-foss-2022b

# Start targets pipeline
R -e "targets::tar_make()"
