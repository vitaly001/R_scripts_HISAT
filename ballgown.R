library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
setwd('/Volumes/HD3/G_RNA/')
gf = "/Volumes/HD3/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"

pheno_data = read.csv("samples.csv")
data_directory = ('ballgown')

file.path(path = dataDir, pattern = samplePattern)

bg = ballgown(dataDir = "ballgown", samplePattern = "S", pData = pheno_data) #worked, FIRST column of pData should match the names of the folders containing the ballgown data and have common attribute (mus in this case) and be in the same oder as it is in the directory
# bg = ballgown(dataDir = "ballgown", samplePattern = "") # alternative way to load the stringtie data into ballgown (without pheno data)
save(bg, file="bg.rda")
load("bg.rda")
#sample_IDs = c("YOUNG_PLUS","YOUNG_MINUS", "OLD_PLUS", "OLD_MINUS")
#sample_paths = paste0('/Volumes/HD2/HISAT_OUTPUT/ballgown/', sample_IDs)
#sample_paths
#sessionInfo() to check the versions of the R packages
annot = gffReadGR(gf, splitByTranscript = TRUE)
checkAssembledTx(annotated = annot, assembled = structure(bg)$trans, ind=23425) # ind integer index of annotated specifying which annotated transcript to plot. All transcripts (assembled and annotated) overlapping annotated[[ind]] will be plotted.
#data_directory = system.file('ballgown', package='ballgown')

#data_directory

gene_expression = gexpr(bg)
gene_expression = data.frame(geneNames=ballgown::geneNames(bg), geneIDs=ballgown::geneIDs(bg), gene_expression)

write.csv(gene_expression, "gene_expression.csv", row.names=TRUE)
full_table = texpr(bg, 'all')
write.csv(full_table, "full_table.csv", row.names=TRUE)
#bg_samples = ballgown(dataDir = data_directory, samplePattern = "S", pData=pheno_data, meas = 'all')

exon_data_frame = eexpr(bg, 'all')
write.csv(exon_data_frame, "exon_data_frame.csv", row.names=TRUE)





bg_filt = subset(bg,"rowVars(texpr(bg)) >5",genomesubset=TRUE)

results_transcripts = stattest(bg_filt, feature="transcript", covariate="Treatment", adjustvars = c("Age"), getFC=TRUE, meas="FPKM")

timecourse_results = stattest(bg_filt, feature='transcript',getFC=TRUE, meas='FPKM',covariate='Treatment', timecourse=TRUE)
timecourse_results = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), timecourse_results)
timecourse_results = arrange(timecourse_results,fc)
write.csv(timecourse_results, "timecourse_results.csv", row.names=FALSE)
results_genes = stattest(bg_filt, feature="gene", covariate="Treatment", adjustvars = c("Age"), getFC=TRUE, meas="FPKM")

#results_genes = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_genes)

results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt), geneIDs=ballgown::geneIDs(bg_filt), results_transcripts)
results_transcripts = arrange(results_transcripts,fc)
results_genes = arrange(results_genes,fc)
write.csv(results_transcripts, "transcript_results.csv", row.names=FALSE)
write.csv(results_genes, "gene_results.csv", row.names=FALSE)
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

tropical= c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
fpkm = texpr(bg, meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$conditions),las=2,ylab='log2(FPKM+1)')
ballgown::transcriptNames(bg)[54295]
ballgown::geneNames(bg)[54295]
plot(fpkm[36295,] ~ pheno_data$conditions, border=c(1,2), main=paste(ballgown::geneNames(bg)[36295],' : ', ballgown::transcriptNames(bg)[36295]),pch=19, xlab="conditions", ylab='log2(FPKM+1)')
points(fpkm[36294,] ~ jitter(as.numeric(pheno_data$Treatment)), col=as.numeric(pheno_data$Treatment))
over100 = exprfilter(bg, cutoff=100, meas = 'FPKM') 
plotTranscripts(ballgown::geneIDs(bg)[36294], bg, main=c('Gene Cdkn2a in sample OLD_PLUS'), samples = c('OLD_PLUS', 'OLD_MINUS', 'YOUNG_PLUS', 'YOUNG_MINUS'))
plotTranscripts(gene="Cdkn2a", gown= bg, samples =c('OLD_PLUS', 'OLD_MINUS', 'YOUNG_PLUS', 'YOUNG_MINUS'), colorby = "transcript",
                meas = "FPKM", legend = TRUE, labelTranscripts = TRUE, main = NULL,
                blackBorders = TRUE, log = TRUE, logbase = 2, customCol = NULL,
                customOrder = NULL)

bg_ob = select(full_table, t_id, gene_name) #require dplyr library

list_of_genes = read.csv("up1.csv", header = FALSE)$V1
as.character(list_of_genes)
list_of_genes = readLines("up.csv") # short way to create list from csv file
pdf_path = '/Volumes/HD3/NGS/G_RNASeq/HISAT_OUTPUT/plot_file.pdf'
pdf(file = pdf_path)

for (i in list_of_genes) {for (n in dat) {plotTranscripts(gene= i, gown= bg, samples =c('OLD_PLUS', 'OLD_MINUS', 'YOUNG_PLUS', 'YOUNG_MINUS'), colorby = "transcript",
                                          meas = "FPKM", legend = TRUE, labelTranscripts = TRUE, main = n,
                                          blackBorders = TRUE, log = FALSE, logbase = 2, customCol = NULL,
                                          customOrder = NULL) }}
dev.off()

lapply(list_of_genes, plotTranscripts, gown=bg, samples = c('OLD_PLUS', 'OLD_MINUS','YOUNG_PLUS', 'YOUNG_MINUS'), labelTranscripts = TRUE)


