#!/usr/bin/env Rscript
#run once to install packages
#source("http://bioconductor.org/biocLite.R")
#biocLite(pkgs=c("Rsubread","limma","edgeR"))

#load libraries
library(Rsubread)
library(edgeR)

#read in fastq files in specified directory
#fastq.files <- list.files(path = "/scratch/PI/dpetrov/grace", pattern = ".fastq$", full.names = TRUE)
#fastq.files

#build index using downloaded reference fasta files from ftp://ftp.ensembl.org/pub/release-93/fasta/saccharomyces_cerevisiae/
#buildindex(basename="fullchr",reference="/scratch/PI/dpetrov/grace/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa")

#map indices
#should be from 33-34 from QC files

#align(index="fullchr",readfile1=fastq.files, phredOffset = 33)

bam.files <- list.files(path = "/scratch/PI/dpetrov/grace", pattern = ".BAM$", full.names = TRUE)

#calculate counts using external annotation gtf file
fc <- featureCounts(bam.files, annot.inbuilt=NULL, annot.ext="/scratch/PI/dpetrov/grace/Saccharomyces_cerevisiae.R64-1-1.84.gtf", isGTFAnnotationFile = TRUE, GTF.featureType="exon",
GTF.attrType="gene_id")


