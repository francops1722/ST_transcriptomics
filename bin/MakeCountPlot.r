#!/usr/bin/env Rscript

require(tidyverse)
require(reshape2)
require(dplyr)
require(biomaRt)

color_panel2<-c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","#ffbf00")
theme_point<-theme_classic()+theme(strip.background = element_blank())

main_dir = commandArgs(trailingOnly=TRUE)
setwd(main_dir)


FC_counts_id<- read.csv("./8_Counts_summary/gene_counts.csv")
colnames(FC_counts_id)[1] <- 'ensembl_gene_id'
ensembl <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", version = 109)
genes_ens_2 <- getBM(attributes=c('ensembl_gene_id','gene_biotype'),mart=ensembl)

genes_long_id <-  FC_counts_id %>% 
  gather(., key='Sample', value="est_counts", -c(ensembl_gene_id)) %>% 
  left_join(., genes_ens_2, by="ensembl_gene_id")

genes_cutoff_id_0 <- genes_long_id %>% filter(est_counts > 0) %>% filter(gene_biotype=="protein_coding")
genes_cutoff_id_10 <- genes_long_id %>% filter(est_counts > 3) %>% filter(gene_biotype=="protein_coding")

number_genes_detected <- genes_cutoff_id_0 %>% group_by(Sample) %>% dplyr::summarize(genes_above0=n())
number_genes_detected_10 <- genes_cutoff_id_10 %>% group_by(Sample) %>% dplyr::summarize(genes_above10=n())

p1 <- ggplot(number_genes_detected,aes(x=Sample,y=genes_above0, col=Sample))+
  geom_point() +
  theme_point+
  theme(panel.grid.major.x=element_line(linetype = "dashed",color="gray88")) +
  labs(x="", y="# protein coding genes",title="Protein coding genes (unfiltered)") +
  scale_y_continuous(limits = c(0, NA))+
  scale_color_manual(values=color_panel2) +
  coord_flip()

p3 <- ggplot(number_genes_detected_10, aes(x=Sample,y=genes_above10, col=Sample))+
  geom_point() +
  theme_point+
  theme(panel.grid.major.x=element_line(linetype = "dashed",color="gray88")) +
  labs(x="", y="# protein coding genes",title="Protein coding genes (>3)") +
  scale_y_continuous(limits = c(0, NA))+
  scale_color_manual(values=color_panel2) +
  coord_flip()

ggsave("ProtCodGenes_noFilt.pdf", plot=p1, height=4, width=6, dpi = 300, useDingbats=F)
ggsave("ProtCodGenes_Filt.pdf", plot=p3, height=4, width=6, dpi = 300, useDingbats=F)