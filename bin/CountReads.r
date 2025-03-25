#!/usr/bin/env Rscript
require(tidyverse)
require(reshape2)
require(dplyr)



main_dir = commandArgs(trailingOnly=TRUE)
path_csv <- list.files(main_dir, recursive=TRUE)%>%str_subset(., "ReadsNumber.csv")

df_reads <- read.csv(paste(main_dir, path_csv, sep="/"), header = T, row.names = 1)

reads_melt <- df_reads %>% pivot_longer(!Sample, names_to="status", values_to="reads")

p1 <- ggplot(reads_melt, aes(x=reorder(status, desc(reads)), y=reads, group = Sample, color=Sample))+
    geom_line()+
      geom_point() +
      scale_y_continuous(trans = "log10") +
      labs(y = "# reads (log10)", x = "stage in pipeline", title = "Reads over the course of pipeline") +
      theme(axis.text.x=element_text(angle=90, hjust=1)) 

ggsave("Reads.pdf", p1, dpi = 300, width = 6, height = 4)

