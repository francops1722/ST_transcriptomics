#!/usr/bin/env Rscript

require(tidyverse)
require(reshape2)
require(dplyr)

color_panel2<-c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","#ffbf00")

main_dir = commandArgs(trailingOnly=TRUE)
setwd(main_dir)

star_path <- "multiqc/star/star_multiqc_report_data/multiqc_star.txt"
star_df <- read_delim(star_path)

reads_df <- read.csv("ReadsNumber.csv", header = T, row.names = 1)
reads_df$Sample <- str_replace(reads_df$Sample, "_R2", "")


Nreads <- merge(reads_df, star_df[, c("Sample", "uniquely_mapped")], by='Sample')
Nreads <- Nreads %>% mutate("map_prct" =round((uniquely_mapped/Raw_Sub)*100,1))


Nreads2_plot<- Nreads %>%
  as_tibble() %>%
  ggplot(aes(x = Sample)) +
  geom_col(aes(y = Raw_Sub, color = "black"), position = "identity", fill = "transparent", width = 0.8)  +
  geom_text(aes(label = "100%", y = Raw_Sub), position = position_stack(vjust = 1), size=2.5)+
  geom_col(aes(y = uniquely_mapped , color = "#428bca"), position = "identity", width = 0.8) + 
  geom_text(aes(label = paste0(map_prct, "%"), y = uniquely_mapped), size=3)+
  theme_bw() + coord_flip()+
  labs(title="Uniquely mapped reads vs Raw reads", y = "# Raw Reads", x = "Sample", subtitle = "The 100% represents the # of subsampled raw reads")+
  scale_color_identity(name = "Status",
                       aesthetics = c("color"),
                       breaks = c("black", "#428bca"),
                       labels = c("Subsampled reads", "Uniquely mapped"), guide = "legend") + scale_fill_manual(values = color_panel2)

ggsave( "UniqvsRaw_subs.pdf", Nreads2_plot, dpi = 300, width = 8, height = 6)



star_data <- star_df %>%  merge(reads_df, by='Sample')

melt_data2 <- star_data %>% 
  pivot_longer(cols=c("uniquely_mapped_percent", "multimapped_percent", "multimapped_toomany_percent", "unmapped_mismatches_percent", "unmapped_tooshort_percent", "unmapped_other_percent"), names_to="status", values_to="reads")


alignment_plot<-melt_data2 %>% ggplot(aes(x=Sample, y=reads, color=status)) + geom_bar(stat='identity', position = 'fill', width = 0.8, fill="transparent") + 
  geom_text(aes(label = ifelse(reads < 2, NA, paste0(round(reads,1), "%")), y = reads), position = position_fill(vjust = 0.5), size=2)+
  theme_bw()+ coord_flip()+
  scale_color_manual(values=color_panel2) +
  labs(title="Alignment of read2 - Subsampling", y = 'Fraction from trimmed reads', x = "Sample", subtitle = "The 100% represent the #reads after 2nd trimming")

ggsave( "Alignment_subs.pdf", alignment_plot, dpi = 300, width = 8, height = 6)