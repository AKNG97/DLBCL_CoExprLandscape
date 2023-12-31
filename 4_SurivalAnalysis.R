#NOTES: The files row_data.RDS and clinical_info.RDS can be obtained through the DataProcessing.R scripts. The 88ConstantlyChangingGenes.RDS file 
# file can be substitued by any expression matrix of interest.

#Load required data

library(dplyr)
library(survival)
library(ggsurvfit)

#Load annot file from TCGA to get gene names
annot <- readRDS("Results/row_data.RDS")

#Load normalized expression data from the set of genes of interest, in this case we 
# analyzed 88 genes with a constant log2 fold change across the sequential order of the COO classification

Norm_data <- readRDS("Results/Norm_data.RDS")
UnclassvsGCB <- read.delim("Results/DEG/resLFC_Group_Unclass_vs_GCB.tsv")
UnclassvsABC <- read.delim("Results/DEG/resLFC_Group_Unclass_vs_ABC.tsv")

up_genes <- unique(intersect(resLFC_UnclassvsGCB[resLFC_UnclassvsGCB$log2FoldChange > 1 & resLFC_UnclassvsGCB$padj < 0.05,]$X,
            resLFC_UnclassvsABC[resLFC_UnclassvsABC$log2FoldChange < -1 & resLFC_UnclassvsABC$padj < 0.05,]$X))

down_genes <- unique(intersect(resLFC_UnclassvsGCB[resLFC_UnclassvsGCB$log2FoldChange < -1 & resLFC_UnclassvsGCB$padj < 0.05,]$X,
                             resLFC_UnclassvsABC[resLFC_UnclassvsABC$log2FoldChange > 1 & resLFC_UnclassvsABC$padj < 0.05,]$X))
genes <- c(up_genes, down_genes)

ProgressExpressGenes <- Norm_data[rownames(Norm_data) %in% genes, ]

saveRDS(Norm_data,"Results/DEG/ConstantlyChangingGenes.RDS")

ProgressExpressGenes <- as.data.frame(ProgressExpressGenes)

#Get gene names
rownames(ProgressExpressGenes) <- annot$gene_name[match(rownames(ProgressExpressGenes), annot$gene_id)]

dir.create("Results/SurvivalAnalysis_results")

#Classify samples as High or Low according with the median of each gene
for(i in 1:nrow(ProgressExpressGenes)){
  
  gene <- rownames(ProgressExpressGenes)[i]
  
  vector_i <- as.data.frame(t(ProgressExpressGenes[i,]))
  
  median_i <- median(vector_i[,1])
  
  low<- rownames(vector_i)[which(vector_i < median_i)]
  high<- rownames(vector_i)[which(vector_i >= median_i)]
  
  factors <- tibble(Group = paste0("High_", gene),
                    Sample = colnames(Lymph_raw))
  factors$Group[which(factors$Sample %in% low)] <- paste0("Low_", gene)
  
  if(i == 1){
    factors_global <- factors
    colnames(factors_global)[1] <- gene
  } else {
    colnames(factors)[1] <- gene
    factors_global <- factors_global %>% inner_join(factors, by = c("Sample" = "Sample"))
  }
  } 

write.csv(factors_global, "Results/SurvivalAnalysis_results/GroupsByMedianExpression.csv", quote = FALSE)

#Get survival data from clinical information
clinical_data <- readRDS("Results/clinical_info.RDS")
clinical_data[is.na(clinical_data$days_to_last_follow_up), ]$days_to_last_follow_up <- 0

factors_global <- factors_global %>% inner_join(clinical_data[,c("sample_submitter_id", "days_to_last_follow_up", "vital_status")], 
                                                           by = c("Sample" = "sample_submitter_id"))

factors_global <- factors_global %>% relocate(c("Sample", "days_to_last_follow_up", "vital_status"))

factors_global <- as.data.frame(factors_global)

#Perform surival analysis
dir.create("Results/SurvivalAnalysis_results/Plots")
for(i in 4:ncol(factors_global)) {
  
  gene.name <- colnames(factors_global)[i]
  gene.name
  factors_global[,i] <- gsub(paste0("Low_", gene.name), 0, factors_global[,i])
  factors_global[,i] <- gsub(paste0("High_", gene.name), 1, factors_global[,i])
  
  surv.matrix <- factors_global[,c(1, 2, 3, i)]
  colnames(surv.matrix)[4] <- "Gene"
  
  surv.plot <- survfit2(Surv(days_to_last_follow_up, vital_status) ~ Gene, data = surv.matrix) %>% 
    ggsurvfit() +
    labs(
      x = "Days to last follow up",
      y = "Overall survival probability",
      title = gene.name
    ) +
    add_risktable()
  
  ggsave(paste0("SurvivalAnalysis_results/Plots/KM_", gene.name, ".png"), surv.plot, height = 10, width = 10, dpi = 300)
  
  suvdiff <- survdiff(Surv(days_to_last_follow_up, vital_status) ~ Gene, data = surv.matrix)
  
  if(suvdiff$pvalue < 0.05) {
    print(paste0(gene.name, " less 0.05"))
    if(i == 4){
      pvalues <- data.frame(Genes = gene.name,
                            KM_pvalue = suvdiff$pvalue)
    } else if(exists("pvalues") == FALSE) {
      pvalues <- data.frame(Genes = gene.name,
                            KM_pvalue = suvdiff$pvalue) 
    } else {
      x <- data.frame(Genes = gene.name,
                      KM_pvalue = suvdiff$pvalue)
      
      pvalues <- rbind(pvalues, x)
    }
  } else {
    print(paste0(gene.name, " over 0.05"))
  }
  
}

write.csv(pvalues, "SurvivalAnalysis_results/SignificantGenesForSurvival.csv", quote = FALSE)

