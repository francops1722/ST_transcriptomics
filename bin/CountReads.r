#!/usr/bin/env Rscript
require(tidyverse)
require(reshape2)
require(dplyr)

color_panel2<-c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","#ffbf00")

main_dir = commandArgs(trailingOnly=TRUE)
setwd(main_dir)


raw <- "multiqc/raw/raw_multiqc_report_data/multiqc_general_stats.txt"
raw_sub <- "multiqc/raw_sub/raw_sub_multiqc_report_data/multiqc_general_stats.txt"
trim1 <- "multiqc/trim/trim_multiqc_report_data/multiqc_general_stats.txt"
demux <- "multiqc/demux/demux_multiqc_report_data/multiqc_general_stats.txt"
trim2 <- "multiqc/trim2/trim2_multiqc_report_data/multiqc_general_stats.txt"

paths <- c(raw, raw_sub, trim1, demux, trim2)


#Reading and formating Raw reads
tmp_raw<-data.table::fread(paths[1], header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
tmp_raw$Sample <- str_replace(tmp_raw$Sample, "_L001", "")
tmp_raw$Sample <- str_replace(tmp_raw$Sample, "_001", "")
colnames(tmp_raw)<- c("Sample", "Raw")
tmp_raw <- tmp_raw[str_detect(tmp_raw$Sample, "_R2"),]

#Reading and formating Subsampled raw reads
if (file.exists(paths[2])){
    status = "subs"
    tmp_sub<-data.table::fread(paths[2], header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
    tmp_sub$Sample <- str_replace(tmp_sub$Sample, "_subs", "")
    colnames(tmp_sub)<- c("Sample", "Raw_Sub")
    tmp_sub <- tmp_sub[str_detect(tmp_sub$Sample, "_R2"),]
    } else {status = "No_subs"}

#Reading and formating trimmed reads
tmp_trim1<-data.table::fread(paths[3], header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
tmp_trim1$Sample <- str_replace(tmp_trim1$Sample, "_i7Trim", "")
colnames(tmp_trim1)<- c("Sample", "Trim1")
tmp_trim1 <- tmp_trim1[str_detect(tmp_trim1$Sample, "_R2"),]

#Reading and formating read after Barcode demultiplexing
tmp_demux<-data.table::fread(paths[4], header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
tmp_demux$Sample <- str_replace(tmp_demux$Sample, "_demux", "")
colnames(tmp_demux)<- c("Sample", "Demux")
tmp_demux <- tmp_demux[str_detect(tmp_demux$Sample, "_R2"),]

#Reading and formating reads after second trimming
tmp_trim2<-data.table::fread(paths[5], header=T, sep="\t", data.table=FALSE)[c("Sample", "FastQC_mqc-generalstats-fastqc-total_sequences")]
tmp_trim2$Sample <- str_replace(tmp_trim2$Sample, "_Trim4", "")
colnames(tmp_trim2)<- c("Sample", "Trim2")

#Merging reads
if (status == "subs") {
    df_reads <- tmp_raw %>% merge(tmp_sub, by='Sample') %>% merge(tmp_trim1, by='Sample') %>% merge(tmp_demux, by='Sample') %>% merge(tmp_trim2, by='Sample')
} else {
    df_reads <- tmp_sub %>% merge(tmp_trim1, by='Sample') %>% merge(tmp_demux, by='Sample') %>% merge(tmp_trim2, by='Sample')
}

#saving files
print("writing csv")
# write out this table to a csv file
write.table(df_reads, "ReadsNumber.csv", sep=",", quote=F, col.names=NA)
#write.csv(df_reads, "./ReadsNumber.csv", quote=FALSE, row.names=FALSE)

#making plot
if (status=="subs"){
    reads_melt <- df_reads %>% pivot_longer(cols=c("Raw", "Raw_Sub", "Trim1", "Demux", "Trim2"), names_to="status", values_to="reads")
} else {
    reads_melt <- df_reads %>% pivot_longer(cols=c("Raw", "Trim1", "Demux", "Trim2"), names_to="status", values_to="reads")
}

p1 <- ggplot(reads_melt,aes(x=reorder(status, desc(reads)), y=reads, group = Sample, color=Sample))+
    geom_line()+
      geom_point() +
      scale_y_continuous(trans = "log10") +
      labs(y = "# reads (log10)", x = "stage in pipeline", title = "Reads over the course of pipeline") +
      theme(axis.text.x=element_text(angle=90, hjust=1))+scale_color_manual(values=color_panel2) 
      
ggsave("Reads.pdf", p1, dpi = 300, width = 6, height = 4)