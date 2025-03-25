#!/usr/bin/env Rscript

require(tidyverse)
require(reshape2)
require(dplyr)



#setting paths
main_dir = commandArgs(trailingOnly=TRUE)
#setwd(main_dir)
#reading report file
path_dedup = paste0(main_dir, "/6_UMI_Counts/logs_umi/UMI_dedup.txt")

if (file.exists(path_dedup)){
    dedup <- read.table(path_dedup)
    colnames(dedup)<-c("Sample", "input_reads", "output_reads")
    dedup <- dedup %>% mutate("percent_passing_dedup" = round(100*output_reads/input_reads, 2))
    dedup$Sample <- sub("_umiCounts.log", "", dedup$Sample)
} else { stop("Gene counts do not exists")}

#making plot

dedup_plot<- dedup %>%
  as_tibble() %>%
  ggplot(aes(x = Sample)) +
  geom_col(aes(y = input_reads, color = "black"), position = "identity", fill = "transparent", width = 0.8)  +
  geom_col(aes(y = output_reads , color = "#428bca"), position = "identity", width = 0.8) + 
  geom_text(aes(label = paste0(percent_passing_dedup, "%"), y = output_reads), size=3)+
  theme_bw() + coord_flip()+
  labs(title="Reads passing deduplication", y = "# Reads", x = "Sample")+
  scale_color_identity(name = "Status",
                       aesthetics = c("color"),
                       breaks = c("black", "#428bca"),
                       labels = c("Input reads", "Output reads"), guide = "legend") 

ggsave("Duplication.pdf", plot=dedup_plot, width = 6, height = 4, dpi = 300, useDingbats=F)
