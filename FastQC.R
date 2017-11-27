# set working directory to the directoy with fastq files
setwd('/Volumes/HD3/NGS/AML')
library("ShortRead")

fqQC = qa(dirPath=".", pattern=".fastq$", type="fastq")
report(fqQC, type="html", dest="fastqQAreport")
