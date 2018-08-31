#!/usr/bin/env bash
#SBATCH -J assemble 1, S1
#SBATCH -p hns,dpetrov
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH -o /scratch/PI/dpetrov/grace/SlurmFiles/slurm-%A_%a.out
#SBATCH --mem=50000
#SBATCH --requeue
#SBATCH --mail-user=gracelam95122@gmail.com
#SBATCH --mail-type=END

/home/groups/dpetrov/SOFTWARE/bin/pear-0.9.6-bin-64 -f /scratch/PI/dpetrov/grace/Grace_Experiment2/GL_081318_1_S97_L001_R1_001.fastq.gz -r /scratch/PI/dpetrov/grace/Grace_Experiment2/GL_081318_1_S97_L001_R2_001.fastq.gz -o /scratch/PI/dpetrov/grace/output


