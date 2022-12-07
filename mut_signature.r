#setwd
#load maf tools
library(maftools)
#download cesc maf.gz file

my_maf = "path to the file"
laml = read.maf(maf = my_maf)

# to see the object of class maf
laml

#Shows sample summry.
getSampleSummary(laml)

#Shows gene summary.
getGeneSummary(laml)

#shows clinical data associated with samples
getClinicalData(laml)

#Shows all fields in MAF
getFields(laml)

#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'laml')

#visualization
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10)

#transitions and transversions
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1))

#mutational signatures
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)

laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

#Signature analysis
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)

plotCophenetic(res = laml.sign)

laml.sig = extractSignatures(mat = laml.tnm, n = 3)

#Compate against original 30 signatures 
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")

#Compate against updated version3 60 signatures 
laml.v3.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "SBS")

install.packages('pheatmap')
library(pheatmap)

pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

#Finally plot signatures
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")
