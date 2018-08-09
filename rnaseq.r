#run once to install packages
source("http://bioconductor.org/biocLite.R")
biocLite(pkgs=c("Rsubread","limma","edgeR"))
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Scerevisiae.UCSC.sacCer3")
install.packages("statmod")
biocLite("topGO")
biocLite("goseq")


#for gene annotations
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")

install.packages(gplots)

#load library
library(Rsubread)
library(edgeR)
library(gplots)
library(goseq)
library(statmod)

#read in fastq files in specified directory
fastq.files <- list.files(path = "./datarna", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files


index.files <- list.files(path = "./datarna", pattern = ".fastq.gz$", full.names = TRUE)

#build index using downloaded reference fasta files from ftp://ftp.ensembl.org/pub/release-93/fasta/saccharomyces_cerevisiae/
buildindex(basename="chr1",reference="./datarna/chrI.fa")

#map indices
#figure out what phredOffset is!
align(index="chr1",readfile1=fastq.files, phredOffset = 64)

bam.files <- list.files(path = "./datarna", pattern = ".BAM$", full.names = TRUE)

#calculate counts using external annotation gtf file
fc <- featureCounts(bam.files, annot.inbuilt=NULL, annot.ext="./RNAseq/ref/Saccharomyces_cerevisiae.R64-1-1.84.gtf", isGTFAnnotationFile = TRUE)

fc <- featureCounts(bam.files, annot.inbuilt=NULL, annot.ext="./RNAseq/ref/Saccharomyces_cerevisiae.R64-1-1.84.gtf", isGTFAnnotationFile = TRUE, GTF.featureType="exon",
GTF.attrType="gene_id")

#create a DGEList object

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID", "Length")])

#Only keep in analysis genes that have >10 reads per million mapped reads in at least two libraries
isexpr <- rowSums(cpm(x) > 10) >= 2

x <- x[isexpr,]


#normalizing, if <1 then small number of high count genes monopolzing sequencing adn causing the counts for other genes to be lower than what would be usual, and library size should be scaled down for that sample.

y <- calcNormFactors(x)
colors <- c("darkgreen", "red", "blue")
pch <- c(0, 1, 2)

plotMDS(x, col=colors, pch=pch)
plotMD(x, column = 1, main="YPD_rep1")
abline(h=0, col="red", lty=2, lwd=2)
plotMD(x, column = 2, main="YPD_rep2")
plotMD(x, column = 3, main="Glycerol")

#heat map clustering
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
col.pan <- colorpanel(100, "blue", "red", "white")

colnames(y) <- c("YPD_rep1", "YPD_rep2", "Glycerol")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
 margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

#checking library sizes
barplot(y$samples$lib.size,main="Barplot of library sizes", names=colnames(y),las=2)

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
#horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#dispersion estimation
design <- model.matrix(~0+y$samples$lib.size)
colnames(design) <- levels(y$samples$lib.size)
y <- estimateDisp(y, design, robust=TRUE)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)


