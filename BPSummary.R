library(dplyr)

communities <-  read.csv("communities_ABC_4.csv", row.names = 1)
communities[communities == "" | communities == " "] <- NA

resLFC_GCBvsABC <- read.delim("resLFC_Group_GCB_vs_ABC.tsv")
resLFC_UnclassvsABC <- read.delim("resLFC_Group_Unclass_vs_ABC.tsv")

resLFC_ABCvsGCB <- resLFC_GCBvsABC %>% select(baseMean, log2FoldChange) %>% mutate(log2FoldChange = -log2FoldChange)
resLFC_ABCvsUnclass <- resLFC_UnclassvsABC %>% select(baseMean, log2FoldChange) %>% mutate(log2FoldChange = -log2FoldChange)

for(i in 1:nrow(communities)) {
  
  query_GCB <- resLFC_ABCvsGCB[rownames(resLFC_ABCvsGCB) %in% communities[i,], ]
  query_ABC <- resLFC_ABCvsUnclass[rownames(resLFC_ABCvsUnclass) %in% communities[i,], ]
  
  CommNumGenes <- length(which(!is.na(communities[i,])))
  
  Num_Diff_up_GCB <- sum(query_GCB$log2FoldChange > 0)
  Num_Diff_down_GCB <- sum(query_GCB$log2FoldChange < 0)

  Num_Diff_up_Unclass <- sum(query_Unclass$log2FoldChange > 0)
  Num_Diff_down_Unclass <- sum(query_Unclass$log2FoldChange < 0)
  
  DE_Comm_GCB <- data.frame(Community_name = rownames(communities[i,]), 
                        Diff_up = Num_Diff_up_GCB/CommNumGenes, 
                        Diff_down = Num_Diff_down_GCB/CommNumGenes,
                        Num_Diff_up = Num_Diff_up_GCB,
                        Num_Diff_down = Num_Diff_down_GCB,
                        NA_DEG = (1 - ((Num_Diff_up_GCB + Num_Diff_down_GCB)/CommNumGenes)),
                        CommNumGenes = CommNumGenes,
                        DENumGenes = nrow(query))

    DE_Comm_Unclass <- data.frame(Community_name = rownames(communities[i,]), 
                        Diff_up = Num_Diff_up_Unclass/CommNumGenes, 
                        Diff_down = Num_Diff_down_Unclass/CommNumGenes,
                        Num_Diff_up = Num_Diff_up_Unclass,
                        Num_Diff_down = Num_Diff_down_Unclass,
                        NA_DEG = (1 - ((Num_Diff_up_Unclass + Num_Diff_down_Unclass)/CommNumGenes)),
                        CommNumGenes = CommNumGenes,
                        DENumGenes = nrow(query))
  
  if(i==1){
    x <- DE_Comm_GCB
    y <- DE_Comm_Unclass
  } else {
    x <- rbind(x, DE_Comm_GCB)
    y <- rbind(x, DE_Comm_Unclass)
  } 
}

Comm_DE_ABCvsGCB <- x
Comm_DE_ABCvsUnclass <- y

ABC_enriched <- read.delim("ABC_Communities_Enriched_4.sif", header = F)
colnames(ABC_enriched) <- c("Community", "Interaction", "GO_BP", "ID", "GeneRatio",  
                "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")

ABC_enriched <- ABC_enriched %>% inner_join(Comm_DE_ABCvsGCB, by = c("Community" = "Community_name"))
ABC_enriched <- ABC_enriched %>% inner_join(Comm_DE_ABCvsUnclass, by = c("Community" = "Community_name"), suffix = c("_GCB", "_Unclass"))
                                      
ABC_enriched <- ABC_enriched %>% group_by(GO_BP) %>% summarise(mean_up_GCB = mean(Diff_up_GCB), mean_down_GCB =  mean(Diff_down_GCB),
                                                 Global_sign_GCB = mean_up_GCB - mean_down_GCB,
                                                 mean_up_Unclass = mean(Diff_up_Unclass), mean_down_Unclass =  mean(Diff_down_Unclass),
                                                 Global_sign_Unclass = mean_up_Unclass - mean_down_Unclass
                                                 )

HT <- as.matrix(ABC_enriched[,c(4,7)])
rownames(HT) <- ABC_enriched$GO_BP         
ABC_HT <- ComplexHeatmap::Heatmap(HT, row_km = 6, row_km_repeats = 100,border = TRUE, heatmap_width = unit(9, "cm"),
                        show_column_names = T,
                        show_row_names = T, row_names_gp = grid::gpar(fontsize = 3),
                        col = colorRamp2(c(-1, 0, 1), 
                                         c("blue", "white", "red")))

draw(ABC_HT)
png("ABC_BP_DEG.png", units="in", width=10, height=10, res=300)
draw(ABC_HT, heatmap_legend_side = "left")
dev.off()

