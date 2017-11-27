
#Using R string manipulation, construct the Unix commands to call stringtie
getwd()

samples = read.csv("samples.csv", stringsAsFactors=FALSE)
gf = "/Volumes/HD3/NGS/Bowtie_HISAT_index_and_genome/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf"
hisatind = "/Volumes/HD2/mm10/genome"

# 1st step
for(i in seq_len(nrow(samples))) {
        lib = samples$LibraryName[i]
        stringtie_cmd = (paste0("stringtie -p 6 -G " , gf, " -o ", lib, '.gtf', ' -l ', lib, ' ', lib,'.bam', ' -A ', lib,'.tab'))
        print(stringtie_cmd)
        system(stringtie_cmd) # invoke command from R directly
        mergelist = paste0(lib,".gtf", collapse = '')
        cat(mergelist, file = "mergelist.txt", collapse = '\n', append = TRUE)
}
# IMPORTANT do not use samples$conditions instead of lib variable (otherwise vector will be created with list of commands)


#2nd step

system(paste0('stringtie --merge -p 6 -G ', gf, ' -o stringtie_merged.gtf mergelist.txt'))

#3rd step
for(i in seq_len(nrow(samples))) {
        lib = samples$LibraryName[i]
        stringtie_ballgown = paste0("stringtie -e -B -p 6 -G stringtie_merged.gtf -o ballgown/", lib, '/',lib,'.gtf ', lib,'.bam')
        system(stringtie_ballgown)
}

print(stringtie_ballgown)
