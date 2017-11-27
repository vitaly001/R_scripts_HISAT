#Align the reads (using tophat2) to reference genome
#Using R string manipulation, construct the Unix commands to call tophat2
getwd()
setwd("/Volumes/HD2/G_RNASeq/")
samples = read.csv("samples.csv", stringsAsFactors=FALSE)
rm(hisat_com)
hisat_com = ''
#gf = "/Volumes/HD2/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
hisatind = "/Volumes/HD2/mm10/genome"

for(i in seq_len(nrow(samples))) {
        cmd = paste0("hisat2 -p 6 --dta -x ", hisatind, " -U ", samples$fastq1, ' -S ', samples$conditions,'.sam')
                hisat_com = paste(hisat_com, cmd, collapse = ' && ' )

}

print(hisat_com)

print(cmd)

#align = do.call(paste, c(as.list(cmd), sep =' && ')) #concatenate all commands in vector separated with shell control operator

#align # print the consecutive commands && : Used to build AND lists, it allows you to run one command only if another exited successfully. 
#system(align)  # invoke commands
# system(cmd) # invoke command
#for (i in cmd) {system(i)}  # invoke commands using loop (different method)
