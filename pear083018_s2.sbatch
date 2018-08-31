#!/usr/bin/env bash
#SBATCH -J assemble 20, S2
#SBATCH -p hns,dpetrov
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH -o /scratch/PI/dpetrov/grace/SlurmFiles/slurm-%A_%a.out
#SBATCH --mem=50000
#SBATCH --requeue
#SBATCH --mail-user=gracelam95122@gmail.com
#SBATCH --mail-type=END

/scratch/PI/dpetrov/pear -f /scratch/PI/dpetrov/grace/Grace_Experiment2/GL_081318_20_S98_L001_R1_001.fastq.gz -r /scratch/PI/dpetrov/grace/Grace_Experiment2/GL_081318_20_S98_L001_R2_001.fastq.gz -o /scratch/PI/dpetrov/grace/output

