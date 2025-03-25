#!/usr/bin/env Rscript
require(tidyverse)
require(reshape2)
require(dplyr)


main_dir = commandArgs(trailingOnly=TRUE)

path0 <- paste0(main_dir,"/multiqc/")

file_list <- list.files(path0, recursive=T)
files<-file_list[grep("./multiqc_general_stats.txt", file_list)]

datasets <- list()

for(i in 1:length(files)){
  if (file.exists(paste0(path0, files[i])) && str_detect(files[i], "star|FCounts|raw\\b", negate = T)){
    tmp<-data.table::fread(paste0(path0, files[i]), header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
    step <- sub("^(\\w+).*", "\\1", files[i])
    colnames(tmp)<- c("Sample", step)
    tmp <- tmp[str_detect(tmp$Sample, "_R2"),]
    tmp$Sample <- sub("_L001.*$|_subs.*$|_demux.*$|_i7Trim.*$|_trim2.*$", "", tmp$Sample)
    tmp$Sample <- sub("_R2.*$", "", tmp$Sample)
    datasets[[step]]<-tmp
  } 
}

df_reads <- Reduce(function(x, y) merge(x, y, all=TRUE), datasets)
#saving files
print("writing csv")
# write out this table to a csv file
write.table(df_reads, "./ReadsNumber.csv", sep=",", quote=F, col.names=NA)