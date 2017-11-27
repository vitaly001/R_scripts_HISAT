
#Using R string manipulation, construct the Unix commands to call stringtie
getwd()
setwd("/Volumes/HD3/G_RNA/")
samples = read.csv("samples.csv", stringsAsFactors=FALSE)
gf = "/Volumes/HD3/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
#hisatind = "/Volumes/HD2/mm10/genome"

# 1st step
for(i in seq_len(nrow(samples))) {
        lib = samples$conditions[i]
        stringtie_cmd = (paste0("stringtie -p 6 -G " , gf, " -o ", lib, '.gtf', ' ', lib,'/accepted_hits.bam'))
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
        lib = samples$conditions[i]
        stringtie_ballgown = paste0("stringtie -e -B -p 6 -G stringtie_merged.gtf -o ballgown/", lib, '/',lib,'.gtf ', lib,'/accepted_hits.bam -A ', lib,'.tab')
        system(stringtie_ballgown)
}

print(stringtie_ballgown)
