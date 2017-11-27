#Align the reads (using HISAT2) to reference genome
#Using R string manipulation, construct the Unix commands to call HISAT2
getwd()
setwd('/Volumes/HD3/NGS/AML')
samples = read.csv("samples.csv", stringsAsFactors=FALSE)

#gf = "/Volumes/HD2/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
hisatind = "/Volumes/HD3/NGS/Bowtie_HISAT_index_and_genome/hg19/genome"
# these arguments (ignore quala and score-min are not optimal amd very loosely, do not use)
cmd = with(samples, print(paste0("hisat2  --ignore-quals --score-min L,0,-1.5 -p 6 --dta -x ", hisatind, " -U ", samples$conditions,'.fastq',' -S ', samples$conditions,'.sam')))

for (i in cmd) {system(i)
        }
