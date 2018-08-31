#!/usr/bin/env bash
#SBATCH -J r_analysis
#SBATCH -p hns,dpetrov
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH -o /scratch/PI/dpetrov/grace/SlurmFiles/slurm-%A_%a.out
#SBATCH --mem=50000
#SBATCH --requeue
#SBATCH --mail-user=gracelam95122@gmail.com
#SBATCH --mail-type=END

module load R
Rscript rna_seq_sherlock.r
