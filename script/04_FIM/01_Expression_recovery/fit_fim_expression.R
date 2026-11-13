#!/usr/bin/env Rscript
# Chicken Productive Aging Atlas | 04_FIM
# Purpose: FIM expression changes and recovery.
# Inputs/outputs: see the module README and config/runtime_layout.json.
# Parameters: retain the study settings; supply tissue/sample names through arguments or manifests.

###########################################

##              DEseq2 for Aging Tissues                   ##

##              created by chzhong 3.22/2024             ##

##              version.1.0                                          ##

###########################################

setwd("/path/to/chicken-aging/ont/shell/15.Transcripts_DESeq/DETs") # Configure paths.

if (length(commandArgs(trailingOnly = TRUE)) != 1) {

  stop("Usage: Rscript fit_fim_expression.R tissue")

}

file1 <- commandArgs(trailingOnly = TRUE)[1]

# Load packages.

library(DESeq2)

library(ggplot2)

# Read the expression matrix.

count_data <- read.table(paste(file1,"_counts_final.txt",sep=""),header=TRUE,stringsAsFactors = TRUE,row.names = 1)

batch <- read.table("group.txt",header=TRUE,stringsAsFactors = TRUE)

#all(rownames(batch) %in% colnames(count_data))

#all(rownames(batch) == colnames(count_data))

ddsHTSeq=DESeqDataSetFromMatrix(countData=count_data, colData=batch, design=~group)

dds<-DESeq(ddsHTSeq,modelMatrixType="standard")

###### 100vsM

res_100vsM<-results(dds,contrast=c("group","100wk","M"))

res_100vsM <- res_100vsM[order(res_100vsM$padj),]

res_100vsM <- na.omit(res_100vsM)

write.table(res_100vsM,paste(file1,"_100vsM.txt",sep=""),col.names = T,row.names = T,quote = F,sep = "\t")

# Select features with adjusted P < 0.05 and absolute log2 fold change > 1.

diff_gene_deseq2 <-subset(res_100vsM,padj < 0.05 & (log2FoldChange>1 | log2FoldChange < -1))

cut_off_pvalue = 0.05  # Statistical significance.

cut_off_logFC = 1      # Expression fold change.

# Assign 100wk_high, CR_high or Stable according to the thresholds; retain legacy field labels.

DEG100vsM <- as.data.frame(diff_gene_deseq2)

DEG100vsM$change = ifelse(DEG100vsM$padj < cut_off_pvalue & abs(DEG100vsM$log2FoldChange) >= cut_off_logFC, 

                           ifelse(DEG100vsM$log2FoldChange > cut_off_logFC ,'100wk_high','CR_high'),

                           'Stable')

write.table(DEG100vsM,paste(file1,"_100vsM_DEG.txt",sep=""),col.names = T,row.names = T,quote = F,sep = "\t")

p <- ggplot(

  # Data, mappings and colors.

  DEG100vsM, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +

  geom_point(alpha=0.4, size=3.5) +

  scale_color_manual(values=c("#546de5","#ff4757"))+

  # Reference lines.

  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +

  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +

  # Axes.

  labs(x="log2(fold change)",

       y="-log10 (p-value)")+

  theme_bw()+

  # Legend.

  theme(plot.title = element_text(hjust = 0.5), 

        legend.position="right", 

        legend.title = element_blank())

#ggsave("output_graph.png", plot=last_plot(), width=10, height=5)

savename = paste("/path/to/chicken-aging/ont/shell/15.Transcripts_DESeq/DETs/",file1,"/100vsM.pdf",sep="")

ggsave(savename, p, width = 120, height = 100, units = "mm")

