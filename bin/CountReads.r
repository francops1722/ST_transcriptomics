#!/usr/bin/env Rscript
require(tidyverse)
require(reshape2)
require(dplyr)

color_panel2<-c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","#ffbf00")

main_dir = commandArgs(trailingOnly=TRUE)
setwd(main_dir)

path0 <- "./multiqc/"

file_list <- list.files(path0, recursive=T)
files<-file_list[grep("./multiqc_general_stats.txt", file_list)]

datasets <- list()

for(i in 1:length(files)){
  if (file.exists(paste0(path0, files[i])) && str_detect(files[i], "star|FCounts", negate = T)){
    tmp<-data.table::fread(paste0(path0, files[i]), header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
    step <- sub("^(\\w+).*", "\\1", files[i])
    colnames(tmp)<- c("Sample", step)
    tmp <- tmp[str_detect(tmp$Sample, "_R2"),]
    tmp$Sample <- sub("_L001.*$|_subs.*$|_demux.*$|_i7Trim.*$|_Trim4.*$", "", tmp$Sample)
    tmp$Sample <- sub("_R2.*$", "", tmp$Sample)
    datasets[[step]]<-tmp
  } 
}

df_reads <- Reduce(function(x, y) merge(x, y, all=TRUE), datasets)
#saving files
print("writing csv")
# write out this table to a csv file
write.table(df_reads, "./ReadsNumber.csv", sep=",", quote=F, col.names=NA)



reads_melt <- df_reads %>% pivot_longer(!Sample, names_to="status", values_to="reads")

p1 <- ggplot(reads_melt, aes(x=reorder(status, desc(reads)), y=reads, group = Sample, color=Sample))+
    geom_line()+
      geom_point() +
      scale_y_continuous(trans = "log10") +
      labs(y = "# reads (log10)", x = "stage in pipeline", title = "Reads over the course of pipeline") +
      theme(axis.text.x=element_text(angle=90, hjust=1))+scale_color_manual(values=color_panel2) 

ggsave("./Reads.pdf", p1, dpi = 300, width = 6, height = 4)

