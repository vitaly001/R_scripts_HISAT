setwd('/Volumes/HD3/NGS/AML')
samples = read.csv("samples.csv", stringsAsFactors=FALSE)

for(i in seq_len(nrow(samples))) {
#        lib = samples$conditions[i]
#        samFile = paste0(lib, ".sam")
        # sort by position and index for IGV and create bam files
#        system(paste0("samtools sort -o ",lib,".bam ", samFile))
        print(i)
}

cmd = with(samples, print(paste0("samtools sort -o ",samples$LibraryName,".bam ", samples$LibraryName,'.sam')))

for (i in cmd) {system(i)
}
